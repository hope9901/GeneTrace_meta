"""Audit NCBI current assembly files for species_taxid=4932 (not historical versions).

Fetch only small md5 manifests, not thousands of genomes. Persist every response.
"""
import argparse
import collections
import concurrent.futures
import csv
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import sha256, write_json

REQUIRED = ('genomic.gff.gz', 'cds_from_genomic.fna.gz', 'protein.faa.gz')


def load(path, database):
    with open(path) as handle:
        header = next(line for line in handle if line.startswith('#assembly_accession'))
        rows = list(csv.DictReader(handle, fieldnames=header.strip().lstrip('#').split('\t'), delimiter='\t'))
    return [dict(row, database=database) for row in rows if row['species_taxid'] == '4932'
            and row['version_status'] == 'latest']


def check(row, cache):
    accession = row['assembly_accession']
    cached = cache / (accession + '.md5.txt')
    url = row['ftp_path'].rstrip('/') + '/md5checksums.txt'
    last_error = None
    for attempt in range(3):
        try:
            if cached.exists():
                body = cached.read_text()
            else:
                with urllib.request.urlopen(url, timeout=30) as response:
                    body = response.read().decode()
                cached.write_text(body)
            # Match basenames exactly, not a substring such as protein.gpff.
            names = {Path(line.split()[-1]).name for line in body.splitlines() if line.strip()}
            prefix = row['ftp_path'].rstrip('/').split('/')[-1] + '_'
            available = {suffix: prefix + suffix in names for suffix in REQUIRED}
            return dict(row, file_check_status='ok', files=available,
                        all_three=all(available.values()), manifest_url=url,
                        manifest_sha256=sha256(cached))
        except Exception as error:
            last_error = str(error)
            time.sleep(attempt + 1)
    return dict(row, file_check_status='error', all_three=None, manifest_url=url, error=last_error)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank', required=True)
    parser.add_argument('--refseq', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    cache = output / 'file_manifests'
    cache.mkdir(exist_ok=True)
    rows = load(args.genbank, 'GenBank') + load(args.refseq, 'RefSeq')
    write_json(output / 'assemblies.json', rows)
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
        futures = [executor.submit(check, row, cache) for row in rows]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
            if len(results) % 100 == 0:
                print('Checked %d/%d manifests (%d errors)' %
                      (len(results), len(rows), sum(r['file_check_status'] == 'error' for r in results)), flush=True)
    results.sort(key=lambda r: r['assembly_accession'])
    write_json(output / 'file_audit.json', results)
    complete = [r for r in results if r['all_three']]
    # Treat GCF/GCA paired records as one assembly, even when one side lacks files.
    def canonical(row):
        if row['database'] == 'RefSeq' and row['gbrs_paired_asm'] != 'na':
            return row['gbrs_paired_asm']
        return row['assembly_accession']
    unique = {canonical(row) for row in complete}
    unique_rows = {canonical(row): row for row in complete}
    summary = dict(checked_utc=datetime.now(timezone.utc).isoformat(), species_taxid=4932,
                   scope='Latest live NCBI fungi assemblies, species_taxid=4932; includes partial assemblies',
                   required_files=REQUIRED, total_records=len(rows),
                   by_database=dict(collections.Counter(r['database'] for r in rows)),
                   eligible_records=len(complete), eligible_by_database=dict(collections.Counter(r['database'] for r in complete)),
                   eligible_unique_assemblies=len(unique),
                   eligible_unique_full_assemblies=sum(r['genome_rep'] == 'Full' for r in unique_rows.values()),
                   eligible_unique_by_level=dict(collections.Counter(r['assembly_level'] for r in unique_rows.values())),
                   errors=[r['assembly_accession'] for r in results if r['file_check_status'] == 'error'],
                   metadata_annotated_genbank=sum(r['database'] == 'GenBank' and r['annotation_provider'] != 'na' for r in rows),
                   summary_sha256=dict(genbank=sha256(args.genbank), refseq=sha256(args.refseq)))
    write_json(output / 'summary.json', summary)
    fields = ['assembly_accession', 'database', 'organism_name', 'infraspecific_name', 'assembly_level',
              'genome_rep', 'annotation_provider', 'protein_coding_gene_count', 'gbrs_paired_asm',
              'file_check_status', 'all_three', 'manifest_url']
    with (output / 'eligible_assemblies.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(complete)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
