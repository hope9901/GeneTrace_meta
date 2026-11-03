"""Create explicitly labeled local-locus controls from public NCBI files."""
import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import read_fasta, read_gff, write_fasta, write_json, sha256

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/'


def main():
    source, output = ROOT / 'data/source', ROOT / 'data/demo'
    output.mkdir(parents=True, exist_ok=True)
    features = read_gff(source / 'genomic.gff.gz')
    genes = [f for f in features if f['type'] == 'gene' and f['attrs'].get('gene') == 'ACT1']
    if len(genes) != 1:
        raise ValueError('Expected exactly one ACT1 gene')
    gene = genes[0]
    locus = gene['attrs']['locus_tag']
    cds = [f for f in features if f['type'] == 'CDS' and f['attrs'].get('locus_tag') == locus]
    proteins = {f['attrs']['protein_id'] for f in cds}
    if len(proteins) != 1:
        raise ValueError('Expected exactly one reference protein')
    protein_id = proteins.pop()
    genome = read_fasta(source / 'genomic.fna.gz')
    proteome = read_fasta(source / 'protein.faa.gz')
    chromosome = genome[gene['seqid']]
    left, right = max(0, gene['start'] - 1 - 15000), min(len(chromosome), gene['end'] + 15000)
    sequence = chromosome[left:right]
    selected = [f for f in features if f['seqid'] == gene['seqid']
                and left < f['start'] <= f['end'] <= right
                and f['type'] not in ('region', 'chromosome')]
    # Keep only fully contained models to avoid dangling Parent references at crop edges.
    complete_genes = {f['attrs'].get('locus_tag') for f in selected if f['type'] == 'gene'}
    selected = [f for f in selected if f['attrs'].get('locus_tag') in complete_genes]
    start, end = gene['start'] - left, gene['end'] - left
    scenarios = {
        'intact': ('Unmodified annotation within a cropped public locus', sequence, False),
        'annotation_removed': ('Controlled perturbation: ACT1 annotation and protein removed; DNA unchanged', sequence, True),
        'assembly_gap': ('Controlled perturbation: ACT1 annotation/protein removed and locus masked with N',
                         sequence[:start - 1] + 'N' * (end - start + 1) + sequence[end:], True),
        'sequence_absent': ('Controlled perturbation: ACT1 locus replaced with A; not a natural biological loss',
                            sequence[:start - 1] + 'A' * (end - start + 1) + sequence[end:], True),
    }
    for scenario, (description, dna, remove) in scenarios.items():
        directory = output / scenario
        directory.mkdir(exist_ok=True)
        write_fasta(directory / 'genome.fna', {gene['seqid']: dna})
        write_fasta(directory / 'query.faa', {protein_id: proteome[protein_id]})
        write_fasta(directory / 'proteins.faa', {k: v for k, v in proteome.items() if not (remove and k == protein_id)})
        with open(directory / 'annotation.gff3', 'w') as handle:
            handle.write('##gff-version 3\n')
            for f in selected:
                if remove and f['attrs'].get('locus_tag') == locus:
                    continue
                fields = list(f['fields'])
                fields[3], fields[4] = str(f['start'] - left), str(f['end'] - left)
                handle.write('\t'.join(fields) + '\n')
        manifest = dict(schema_version=1, dataset='S. cerevisiae S288C ACT1 local-locus controls',
                        assembly='GCF_000146045.2', gene='ACT1', locus_tag=locus,
                        scenario=scenario, description=description,
                        controlled_perturbation=remove, scope='Cropped single-reference locus, not comparative genomes',
                        reference_chromosome=gene['seqid'], reference_offset_0based=left,
                        expected_locus=dict(seqid=gene['seqid'], start=start, end=end, strand=gene['strand']),
                        coordinate_system='1-based inclusive; add reference_offset_0based for original coordinates',
                        sources=[dict(file=name, url=BASE + 'GCF_000146045.2_R64_' + name,
                                      sha256=sha256(source / name))
                                 for name in ('genomic.fna.gz', 'genomic.gff.gz', 'protein.faa.gz')],
                        retrieval_recorded_utc=datetime.now(timezone.utc).isoformat(),
                        data_policy='https://www.ncbi.nlm.nih.gov/home/about/policies/',
                        files={name: sha256(directory / name) for name in
                               ('genome.fna', 'query.faa', 'proteins.faa', 'annotation.gff3')})
        write_json(directory / 'manifest.json', manifest)
    print('Prepared four ACT1 controls in ' + str(output))


if __name__ == '__main__':
    main()
