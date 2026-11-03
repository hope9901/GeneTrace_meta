"""Cross-strain ACT1 recovery from complete, unannotated NCBI contig assemblies.

No annotation is invented. Perturbations test the hit-seeded locus separately
from homologs elsewhere in the assembly. This is not an orthology benchmark.
"""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.engine import EvidenceStore, FIELDS, THRESHOLDS, find_binary, parse_hits
from genetrace.io import read_fasta, write_fasta, write_json, sha256

ROOT = Path(__file__).resolve().parents[1]
STRAINS = {
    'CLIB215': ('GCA_019394775.1', 'ASM1939477v1', '019/394/775'),
    'KACC30049': ('GCA_060154965.1', 'ASM6015496v1', '060/154/965'),
}


def assembly_metrics(sequences):
    lengths = sorted((len(s) for s in sequences.values()), reverse=True)
    total, cumulative, n50 = sum(lengths), 0, 0
    for length in lengths:
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return dict(sequence_count=len(lengths), total_bp=total, n50_bp=n50,
                n_fraction=sum(s.count('N') for s in sequences.values()) / total)


def overlap(hit, locus):
    return (hit['sseqid'] == locus['seqid'] and hit['strand'] == locus['strand']
            and hit['start'] <= locus['end'] and hit['end'] >= locus['start'])


def search(store, genome, query):
    binary = find_binary('tblastn')
    version = store.command('tblastn_version', [binary, '-version'])
    record = store.command('tblastn', [binary, '-query', str(query.resolve()),
                                     '-subject', str(genome.resolve()), '-evalue', '1e-5',
                                     '-outfmt', '6 ' + FIELDS, '-max_target_seqs', '100', '-seg', 'yes'])
    if record['status'] != 'ok':
        raise RuntimeError('Translated search failed; see ' + str(store.directory))
    hits = parse_hits(record['result']['stdout'])
    store.add('parse_alignments', dict(hits=hits, thresholds=THRESHOLDS),
              input_evidence=record['id'], version_evidence=version['id'])
    return hits


def run_strain(strain, root):
    accession, name, path = STRAINS[strain]
    source = ROOT / 'data/strains' / strain
    report_text = (source / 'assembly_report.txt').read_text()
    if '# Assembly level: Contig' not in report_text:
        raise ValueError('NCBI report does not confirm Contig assembly level')
    if accession not in report_text or 'Saccharomyces cerevisiae' not in report_text:
        raise ValueError('NCBI report organism or accession mismatch')
    sequences = read_fasta(source / 'genomic.fna.gz')
    metrics = assembly_metrics(sequences)
    strain_root = root / strain
    store = EvidenceStore(strain_root / 'original')
    query = ROOT / 'data/demo/intact/query.faa'
    write_fasta(store.directory / 'genome.fna', sequences)
    write_fasta(store.directory / 'query.faa', read_fasta(query))
    base = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + path + '/' + accession + '_' + name + '/'
    manifest = dict(strain=strain, accession=accession, assembly_level='Contig',
                    annotation_status='unavailable_in_ncbi_ftp_listing', query='S288C ACT1 NP_116614.1',
                    query_sha256=sha256(query), metrics=metrics,
                    files=[dict(file=p.name, sha256=sha256(p),
                                url=base + accession + '_' + name + '_' + p.name)
                           for p in (source / 'genomic.fna.gz', source / 'assembly_report.txt', source / 'feature_count.txt')],
                    query_selection='Same preselected ACT1 protein as earlier baseline',
                    seed_method='Best observed qualifying TBLASTN HSP, not independent synteny localization',
                    created_utc=datetime.now(timezone.utc).isoformat())
    write_json(strain_root / 'manifest.json', manifest)
    store.add('source_metadata', manifest)
    hits = search(store, store.directory / 'genome.fna', store.directory / 'query.faa')
    qualifying = [h for h in hits if h['passes_thresholds']]
    if not qualifying:
        raise RuntimeError('No qualifying candidate; do not manufacture a seed')
    best = max(qualifying, key=lambda h: h['bitscore'])
    seq = sequences[best['sseqid']]
    # Cover the main HSP plus flanks, including the short exon/intron region.
    locus = dict(seqid=best['sseqid'], start=max(1, best['start'] - 500),
                 end=min(len(seq), best['end'] + 500), strand=best['strand'])
    summary = dict(strain=strain, accession=accession, metrics=metrics,
                   best_hit=best, locus=locus, annotation_status=manifest['annotation_status'],
                   scenarios=[])
    for scenario in ('original', 'locus_N_masked', 'locus_A_replaced'):
        run_store = store if scenario == 'original' else EvidenceStore(strain_root / scenario)
        if scenario == 'original':
            mutated, results = sequences, hits
        else:
            mutated = dict(sequences)
            letter = 'N' if scenario == 'locus_N_masked' else 'A'
            mutated[locus['seqid']] = (seq[:locus['start'] - 1] + letter * (locus['end'] - locus['start'] + 1)
                                       + seq[locus['end']:])
            write_fasta(run_store.directory / 'genome.fna', mutated)
            write_fasta(run_store.directory / 'query.faa', read_fasta(query))
            results = search(run_store, run_store.directory / 'genome.fna', run_store.directory / 'query.faa')
        at_locus = [h for h in results if h['passes_thresholds'] and overlap(h, locus)]
        other = [h for h in results if h['passes_thresholds'] and not overlap(h, locus)]
        locus_seq = mutated[locus['seqid']][locus['start'] - 1:locus['end']]
        edge_distance = min(locus['start'] - 1, len(seq) - locus['end'])
        if at_locus:
            verdict = 'homolog_detected_annotation_unavailable'
        elif locus_seq.count('N') / len(locus_seq) >= .1 or edge_distance < 100:
            verdict = 'assembly_limited'
        else:
            verdict = 'inconclusive'
        run_store.add('locus_assessment', dict(locus=locus, n_fraction=locus_seq.count('N') / len(locus_seq),
                                              distance_to_contig_edge=edge_distance,
                                              qualifying_locus_hits=at_locus, other_qualifying_hits=other,
                                              verdict=verdict))
        scenario_report = dict(scenario=scenario, controlled_perturbation=scenario != 'original',
                               verdict=verdict, locus_hits=len(at_locus), other_qualifying_hits=len(other),
                               best_locus_hit=max(at_locus, key=lambda h: h['bitscore']) if at_locus else None,
                               evidence_ids=[e['id'] for e in run_store.items],
                               genome_sha256=sha256(run_store.directory / 'genome.fna'),
                               limitations=['No original CDS/protein annotation was available; annotation gaps cannot be assessed.',
                                            'Seed is derived from the observed best hit; not a blind absence benchmark.',
                                            'Other actin-related hits are not proof that ACT1 exists at the seeded locus.',
                                            'Homology does not establish orthology or function. No biological loss call.',
                                            'Single HSP thresholds; no spliced alignment or read support. No Meta model call.'])
        write_json(run_store.directory / 'report.json', scenario_report)
        write_json(run_store.directory / 'evidence.json', run_store.items)
        summary['scenarios'].append(scenario_report)
    write_json(strain_root / 'summary.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, help='New output root')
    args = parser.parse_args()
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    summaries = [run_strain(strain, root) for strain in STRAINS]
    write_json(root / 'comparison.json', summaries)
    for summary in summaries:
        hit = summary['best_hit']
        print(summary['strain'], summary['metrics'],
              'identity=', hit['pident'], 'coverage=', hit['query_coverage'],
              [(s['scenario'], s['verdict']) for s in summary['scenarios']])


if __name__ == '__main__':
    main()
