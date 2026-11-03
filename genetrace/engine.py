import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .io import read_fasta, read_gff, sha256, write_json

ROOT = Path(__file__).resolve().parents[1]
FIELDS = 'qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen qseq sseq'
THRESHOLDS = dict(min_identity=40.0, min_single_hsp_query_coverage=0.60, max_evalue=1e-10)


def find_binary(name):
    configured = os.environ.get('BLAST_BIN')
    candidates = ([Path(configured) / name] if configured else [])
    candidates += sorted((ROOT / 'tools').glob('ncbi-blast-*/bin/' + name))
    installed = shutil.which(name)
    if installed:
        candidates.append(Path(installed))
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    raise FileNotFoundError('BLAST+ missing: set BLAST_BIN to a directory containing ' + name)


class EvidenceStore:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=False)
        self.items = []

    def add(self, tool, result, status='ok', **metadata):
        evidence_id = 'E%03d' % (len(self.items) + 1)
        record = dict(id=evidence_id, tool=tool, status=status, result=result, **metadata)
        self.items.append(record)
        write_json(self.directory / (evidence_id + '.json'), record)
        return record

    def command(self, tool, args):
        try:
            process = subprocess.run(args, capture_output=True, text=True, timeout=120, shell=False)
            return self.add(tool, dict(stdout=process.stdout, stderr=process.stderr,
                                      exit_code=process.returncode),
                            status='ok' if process.returncode == 0 else 'error', command=args)
        except (OSError, subprocess.TimeoutExpired) as error:
            return self.add(tool, dict(error=str(error)), status='error', command=args)


def parse_hits(stdout):
    hits = []
    for line in stdout.splitlines():
        values = line.split('\t')
        if len(values) != 13:
            raise ValueError('Unexpected BLAST output column count')
        hit = dict(zip(FIELDS.split(), values))
        for key in ('length', 'qstart', 'qend', 'sstart', 'send', 'qlen'):
            hit[key] = int(hit[key])
        for key in ('pident', 'evalue', 'bitscore'):
            hit[key] = float(hit[key])
        hit['query_coverage'] = (abs(hit['qend'] - hit['qstart']) + 1) / hit['qlen']
        hit['start'], hit['end'] = sorted((hit['sstart'], hit['send']))
        hit['strand'] = '+' if hit['send'] >= hit['sstart'] else '-'
        hit['passes_thresholds'] = (hit['pident'] >= THRESHOLDS['min_identity'] and
                                   hit['query_coverage'] >= THRESHOLDS['min_single_hsp_query_coverage'] and
                                   hit['evalue'] <= THRESHOLDS['max_evalue'])
        hits.append(hit)
    return hits


def classify(hits, features, locus, locus_metrics):
    """Conservative local-locus triage; never infer loss from search absence."""
    strong = [h for h in hits if h['passes_thresholds'] and h['sseqid'] == locus['seqid']
              and h['strand'] == locus['strand'] and h['start'] <= locus['end'] and h['end'] >= locus['start']]
    for hit in strong:
        overlapping = [f for f in features if f['type'] == 'CDS' and f['seqid'] == hit['sseqid']
                       and f['strand'] == hit['strand']
                       and max(0, min(hit['end'], f['end']) - max(hit['start'], f['start']) + 1)
                       / (hit['end'] - hit['start'] + 1) >= 0.5]
        if not overlapping:
            return 'annotation_gap_candidate', 'A strong translated hit overlaps the seeded locus without a same-strand CDS covering at least half of this hit.'
    if strong:
        return 'annotated_homolog', 'Strong sequence evidence at the seeded locus overlaps a same-strand annotated CDS.'
    if locus_metrics['n_fraction'] >= 0.1 or locus_metrics['distance_to_contig_edge'] < 100:
        return 'assembly_limited', 'The seeded locus has substantial ambiguous sequence or is near a contig edge; absence is unresolved.'
    return 'inconclusive', 'No qualifying hit at the seeded locus. This is not evidence sufficient to infer biological gene loss.'


def investigate(dataset, output):
    dataset = Path(dataset).resolve()
    manifest = json.loads((dataset / 'manifest.json').read_text())
    required = ('genome.fna', 'annotation.gff3', 'proteins.faa', 'query.faa')
    if set(manifest['files']) != set(required):
        raise ValueError('Manifest must list exactly the four supported input files')
    for filename in required:
        if sha256(dataset / filename) != manifest['files'][filename]:
            raise ValueError('Input checksum mismatch: ' + filename)
    sequences = read_fasta(dataset / 'genome.fna')
    query = read_fasta(dataset / 'query.faa')
    if len(query) != 1:
        raise ValueError('This MVP supports exactly one query protein')
    features = read_gff(dataset / 'annotation.gff3')
    locus = manifest['expected_locus']
    seq = sequences[locus['seqid']]
    if not (1 <= locus['start'] <= locus['end'] <= len(seq)):
        raise ValueError('Seeded locus outside sequence bounds')
    store = EvidenceStore(output)
    write_json(store.directory / 'input_manifest.json', manifest)
    input_dir = store.directory / 'inputs'
    input_dir.mkdir()
    for filename in required:
        shutil.copyfile(dataset / filename, input_dir / filename)
    annotation = store.add('inspect_annotation', dict(
        gene_count=sum(f['type'] == 'gene' for f in features),
        seeded_locus=locus,
        overlapping_genes=[{k: f[k] for k in ('seqid', 'start', 'end', 'strand', 'attrs')}
                           for f in features if f['type'] == 'gene' and f['seqid'] == locus['seqid']
                           and f['start'] <= locus['end'] and f['end'] >= locus['start']]))
    searches = {}
    for tool, subject in (('blastp', 'proteins.faa'), ('tblastn', 'genome.fna')):
        try:
            binary = find_binary(tool)
        except FileNotFoundError as error:
            searches[tool] = store.add(tool, dict(error=str(error)), status='error')
            continue
        version = store.command(tool + '_version', [binary, '-version'])
        args = [binary, '-query', str((input_dir / 'query.faa').resolve()),
                '-subject', str((input_dir / subject).resolve()), '-evalue', '1e-5',
                '-outfmt', '6 ' + FIELDS, '-max_target_seqs', '100', '-seg', 'yes']
        result = store.command(tool, args)
        if result['status'] == 'ok':
            try:
                result['hits'] = parse_hits(result['result']['stdout'])
                result['version_evidence'] = version['id']
                write_json(store.directory / (result['id'] + '.json'), result)
            except ValueError as error:
                result['status'] = 'error'
                result['parse_error'] = str(error)
                write_json(store.directory / (result['id'] + '.json'), result)
        searches[tool] = result
    locus_seq = seq[locus['start'] - 1:locus['end']]
    metrics = dict(n_fraction=locus_seq.count('N') / len(locus_seq),
                   distance_to_contig_edge=min(locus['start'] - 1, len(seq) - locus['end']),
                   locus_length=len(locus_seq),
                   neighboring_genes=[dict(start=f['start'], end=f['end'], name=f['attrs'].get('Name'))
                                      for f in features if f['type'] == 'gene' and f['seqid'] == locus['seqid']
                                      and f['start'] <= locus['end'] + 5000 and f['end'] >= locus['start'] - 5000])
    local = store.add('inspect_locus', metrics)
    if any(search['status'] != 'ok' for search in searches.values()):
        verdict, explanation = 'inconclusive', 'A required sequence search failed. No biological conclusion is made.'
    else:
        verdict, explanation = classify(searches['tblastn']['hits'], features, locus, metrics)
    report = dict(schema_version=1, created_utc=datetime.now(timezone.utc).isoformat(),
                  execution_mode='deterministic_baseline_no_model', scenario=manifest['scenario'],
                  controlled_perturbation=manifest['controlled_perturbation'], description=manifest['description'],
                  scope=manifest['scope'], verdict=verdict, thresholds=THRESHOLDS,
                  claims=[dict(text=explanation, evidence_ids=[annotation['id'], searches['blastp']['id'],
                                                            searches['tblastn']['id'], local['id']])],
                  limitations=['Single-reference cropped-locus control; not a blind comparative benchmark.',
                               'Expected locus is provided as a seed from source annotation.',
                               'Single-HSP threshold may miss intron-rich or divergent homologs.',
                               'Homology does not establish orthology, intact function, or true gene loss.',
                               'Meta model has not participated in this baseline run.'])
    validate_claims(report, store.items)
    write_json(store.directory / 'evidence.json', store.items)
    write_json(store.directory / 'report.json', report)
    markdown = ['# GeneTrace evidence report', '', '**' + verdict + '**', '', manifest['description'], '',
                'Mode: deterministic baseline (no Meta model call)', '',
                explanation, '', 'Evidence: ' + ', '.join(report['claims'][0]['evidence_ids']), '',
                '## Limits', ''] + ['- ' + item for item in report['limitations']]
    (store.directory / 'report.md').write_text('\n'.join(markdown) + '\n')
    return report


def validate_claims(report, evidence):
    known = {item['id'] for item in evidence}
    for claim in report['claims']:
        if not claim['evidence_ids'] or not set(claim['evidence_ids']) <= known:
            raise ValueError('Claim contains missing or unknown evidence IDs')
