"""Enrich the audited 435-assembly cohort with actual protein FASTA counts.

Compressed protein FASTA is downloaded; --recompute-gc also reads genome FASTA.
NCBI-published MD5 is checked; counts and checksums are cached instead of
retaining hundreds of FASTA files.
"""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import write_json

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / 'data/inventory/s_cerevisiae_2026-10-30'


def count_proteins(row, cache):
    accession = row['assembly_accession']
    basename = row['ftp_path'].rstrip('/').split('/')[-1] + '_protein.faa.gz'
    url = row['ftp_path'].rstrip('/') + '/' + basename
    manifest = (INVENTORY / 'file_manifests' / (accession + '.md5.txt')).read_text()
    expected = next(line.split()[0] for line in manifest.splitlines()
                    if line.strip() and Path(line.split()[-1]).name == basename)
    destination = cache / (accession + '.json')
    if destination.exists():
        cached = json.loads(destination.read_text())
        if cached.get('expected_md5') == expected and cached.get('status') == 'ok':
            return cached
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as response:
                compressed = response.read()
            actual = hashlib.md5(compressed).hexdigest()
            if actual != expected:
                raise ValueError('NCBI protein FASTA snapshot changed: MD5 mismatch')
            text = gzip.decompress(compressed).decode('ascii')
            identifiers, lengths, length = [], [], None
            for line in text.splitlines():
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if length is not None:
                        lengths.append(length)
                    identifiers.append(line[1:].split()[0])
                    length = 0
                elif length is None:
                    raise ValueError('Sequence appears before FASTA header')
                else:
                    length += len(line)
            if length is not None:
                lengths.append(length)
            if not lengths or min(lengths) == 0:
                raise ValueError('No proteins or empty protein record')
            result = dict(assembly_accession=accession, status='ok', number_of_proteins=len(lengths),
                          unique_protein_ids=len(set(identifiers)), total_protein_aa=sum(lengths),
                          compressed_bytes=len(compressed), expected_md5=expected, observed_md5=actual,
                          protein_fasta_sha256=hashlib.sha256(compressed).hexdigest(), protein_fasta_url=url,
                          counted_utc=datetime.now(timezone.utc).isoformat())
            write_json(destination, result)
            return result
        except Exception as error:
            last_error = str(error)
            time.sleep(attempt + 1)
    result = dict(assembly_accession=accession, status='error', error=last_error, protein_fasta_url=url)
    write_json(destination, result)
    return result


def count_genome(row, cache):
    accession = row['assembly_accession']
    basename = row['ftp_path'].rstrip('/').split('/')[-1] + '_genomic.fna.gz'
    url = row['ftp_path'].rstrip('/') + '/' + basename
    manifest = (INVENTORY / 'file_manifests' / (accession + '.md5.txt')).read_text()
    expected = next(line.split()[0] for line in manifest.splitlines()
                    if line.strip() and Path(line.split()[-1]).name == basename)
    destination = cache / (accession + '.json')
    if destination.exists():
        cached = json.loads(destination.read_text())
        if cached.get('expected_md5') == expected and cached.get('status') == 'ok':
            return cached
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                compressed = response.read()
            actual = hashlib.md5(compressed).hexdigest()
            if actual != expected:
                raise ValueError('NCBI genome snapshot changed: MD5 mismatch')
            gc, atgc, total, n_count, sequence_count = 0, 0, 0, 0, 0
            lengths, current_length = [], None
            for line in gzip.decompress(compressed).splitlines():
                if line.startswith(b'>'):
                    if current_length is not None:
                        lengths.append(current_length)
                    current_length = 0
                    sequence_count += 1
                    continue
                line = line.strip().upper()
                if not line:
                    continue
                if current_length is None:
                    raise ValueError('Genome FASTA sequence before header')
                count_gc = line.count(b'G') + line.count(b'C')
                gc += count_gc
                atgc += count_gc + line.count(b'A') + line.count(b'T')
                n_count += line.count(b'N')
                total += len(line)
                current_length += len(line)
            if current_length is not None:
                lengths.append(current_length)
            if not lengths or not atgc or min(lengths) == 0:
                raise ValueError('Empty genome sequence or zero ACGT bases')
            cumulative, n50 = 0, None
            for length in sorted(lengths, reverse=True):
                cumulative += length
                if cumulative >= total / 2:
                    n50 = length
                    break
            result = dict(assembly_accession=accession, status='ok', gc_percent=100 * gc / atgc,
                          gc_bases=gc, atgc_bases=atgc, total_bases=total, n_bases=n_count,
                          number_of_fasta_sequences=sequence_count, fasta_sequence_n50_bp=n50,
                          compressed_bytes=len(compressed), expected_md5=expected, observed_md5=actual,
                          genome_fasta_sha256=hashlib.sha256(compressed).hexdigest(), genome_fasta_url=url,
                          counted_utc=datetime.now(timezone.utc).isoformat())
            write_json(destination, result)
            return result
        except Exception as error:
            last_error = str(error)
            time.sleep(attempt + 1)
    result = dict(assembly_accession=accession, status='error', error=last_error, genome_fasta_url=url)
    write_json(destination, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--recompute-gc', action='store_true', help='Download whole-genome FASTA to compute unrounded GC%')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    cache = output / 'protein_counts'
    cache.mkdir(exist_ok=True)
    cohort = [row for row in json.loads((INVENTORY / 'file_audit.json').read_text())
              if row['all_three'] and row['database'] == 'GenBank']
    if len(cohort) != 435 or len({r['assembly_accession'] for r in cohort}) != 435:
        raise ValueError('Expected the audited, deduplicated 435-assembly snapshot')
    counted = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
        futures = [executor.submit(count_proteins, row, cache) for row in cohort]
        for future in concurrent.futures.as_completed(futures):
            record = future.result()
            counted[record['assembly_accession']] = record
            if len(counted) % 25 == 0:
                print('Protein FASTA checked: %d/%d; errors: %d' %
                      (len(counted), len(cohort), sum(r['status'] != 'ok' for r in counted.values())), flush=True)
    if any(r['status'] != 'ok' for r in counted.values()):
        raise RuntimeError('Failed protein downloads/counts. Inspect cache and rerun; no incomplete table published.')
    genome_counts = {}
    if args.recompute_gc:
        genome_cache = output / 'genome_counts'
        genome_cache.mkdir(exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
            futures = [executor.submit(count_genome, row, genome_cache) for row in cohort]
            for future in concurrent.futures.as_completed(futures):
                record = future.result()
                genome_counts[record['assembly_accession']] = record
                if len(genome_counts) % 25 == 0:
                    print('Genome FASTA checked: %d/%d; errors: %d' %
                          (len(genome_counts), len(cohort), sum(r['status'] != 'ok' for r in genome_counts.values())), flush=True)
        if any(r['status'] != 'ok' for r in genome_counts.values()):
            raise RuntimeError('Genome download/count failed; no incomplete table published. Inspect cache and rerun.')
    table = []
    for row in cohort:
        count = counted[row['assembly_accession']]
        genome_count = genome_counts.get(row['assembly_accession'])
        strain = row['infraspecific_name']
        if strain.startswith('strain='):
            strain = strain[7:]
        def number(key, converter=int):
            return None if row[key] in ('na', '') else converter(row[key])
        result = dict(assembly_accession=row['assembly_accession'], strain=None if strain == 'na' else strain,
                      isolate=None if row['isolate'] == 'na' else row['isolate'], organism=row['organism_name'],
                      assembly_level=row['assembly_level'], genome_representation=row['genome_rep'],
                      included_in_primary_plots=row['genome_rep'] == 'Full',
                      genome_size_bp=number('genome_size'), genome_size_mb=number('genome_size') / 1e6,
                      ungapped_genome_size_bp=number('genome_size_ungapped'), gc_percent=number('gc_percent', float),
                      number_of_proteins=count['number_of_proteins'],
                      protein_coding_gene_count=number('protein_coding_gene_count'),
                      total_gene_count=number('total_gene_count'),
                      contig_count=number('contig_count'), scaffold_count=number('scaffold_count'),
                      assembly_type=None if row['assembly_type'] == 'na' else row['assembly_type'],
                      release_date=row['seq_rel_date'], annotation_date=row['annotation_date'],
                      annotation_provider=row['annotation_provider'], submitter=row['asm_submitter'],
                      biosample=None if row['biosample'] == 'na' else row['biosample'],
                      bioproject=None if row['bioproject'] == 'na' else row['bioproject'],
                      paired_refseq=None if row['gbrs_paired_asm'] == 'na' else row['gbrs_paired_asm'],
                      ncbi_url='https://www.ncbi.nlm.nih.gov/datasets/genome/' + row['assembly_accession'] + '/',
                      metadata_source='https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt',
                      protein_fasta_url=count['protein_fasta_url'], protein_fasta_sha256=count['protein_fasta_sha256'],
                      protein_count_method='FASTA record count; MD5 verified against NCBI manifest')
        if genome_count:
            # NCBI summary size is primary-assembly span, whereas the FASTA can
            # also contain non-nuclear/other assembly units. Preserve both scopes.
            result.update(ncbi_primary_genome_size_bp=result['genome_size_bp'],
                          fasta_size_minus_ncbi_primary_bp=genome_count['total_bases'] - result['genome_size_bp'],
                          genome_size_bp=genome_count['total_bases'], genome_size_mb=genome_count['total_bases'] / 1e6,
                          ncbi_gc_percent=result['gc_percent'], gc_percent=genome_count['gc_percent'],
                          gc_method='100*(G+C)/(A+C+G+T) from full genomic FASTA; N and other ambiguity codes excluded',
                          n_fraction=genome_count['n_bases'] / genome_count['total_bases'],
                          genome_fasta_sequence_count=genome_count['number_of_fasta_sequences'],
                          fasta_sequence_n50_bp=genome_count['fasta_sequence_n50_bp'],
                          genome_fasta_sha256=genome_count['genome_fasta_sha256'],
                          genome_fasta_url=genome_count['genome_fasta_url'])
        table.append(result)
    write_json(output / 'genome_metadata.json', table)
    with (output / 'genome_metadata.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    print('Saved metadata for %d assemblies (%d full representation)' %
          (len(table), sum(r['included_in_primary_plots'] for r in table)), flush=True)


if __name__ == '__main__':
    main()
