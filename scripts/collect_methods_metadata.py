"""Collect explicit sequencing/assembly and structural annotation provenance.

Assembly reports and complete GFF source columns are audited. A bounded first
GenBank-record header is also read for structured annotation metadata; it is
not represented as a complete GenBank-file audit. No literature or SRA inference.
"""
import concurrent.futures
import gzip
import hashlib
import io
import json
import re
import sys
import time
import urllib.request
from urllib.parse import unquote
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import write_json

ROOT = Path(__file__).resolve().parents[1]
INV = ROOT / 'data/inventory/s_cerevisiae_2026-10-30'
OUT = ROOT / 'data/inventory/s_cerevisiae_methods_2026-10-30'


def get_bytes(url):
    with urllib.request.urlopen(url, timeout=90) as response:
        return response.read()


def collect(row):
    accession = row['assembly_accession']
    dest = OUT / accession
    dest.mkdir(parents=True, exist_ok=True)
    result_file = dest / 'raw.json'
    if result_file.exists():
        record = json.loads(result_file.read_text())
        if record.get('status') == 'ok' and 'gff_structural_method_evidence' in record:
            return record
    base = row['ftp_path'].rstrip('/')
    prefix = base.split('/')[-1]
    md5_text = (INV / 'file_manifests' / (accession + '.md5.txt')).read_text()
    expected = {Path(line.split()[-1]).name: line.split()[0] for line in md5_text.splitlines() if line.strip()}

    def verified(suffix):
        name = prefix + '_' + suffix
        url = base + '/' + name
        body = get_bytes(url)
        if hashlib.md5(body).hexdigest() != expected[name]:
            raise ValueError('Source snapshot changed: ' + name)
        return body, url

    for attempt in range(3):
        try:
            report_path = dest / 'assembly_report.txt'
            if report_path.exists():
                report = report_path.read_text()
            else:
                raw, report_url = verified('assembly_report.txt')
                report = raw.decode()
                report_path.write_text(report)
            report_fields = dict(re.findall(r'^# ([^:\n]+):\s*(.*?)\s*$', report, flags=re.M))
            raw_gff, gff_url = verified('genomic.gff.gz')
            comments, sources, method_evidence = [], {}, set()
            for line in gzip.decompress(raw_gff).decode().splitlines():
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#'):
                    if not line.startswith('##sequence-region'):
                        comments.append(line)
                    continue
                parts = line.split('\t')
                if len(parts) == 9 and parts[2] in ('gene', 'mRNA', 'transcript', 'CDS', 'exon'):
                    sources[parts[1]] = sources.get(parts[1], 0) + 1
                    attrs = dict(part.split('=', 1) for part in parts[8].split(';') if '=' in part)
                    inference = unquote(attrs.get('inference', ''))
                    if 'ab initio prediction' in inference.lower():
                        method_evidence.add(inference)
            gbff_url = base + '/' + prefix + '_genomic.gbff.gz'
            # Header only: exclude feature products/functions and record sequence.
            header_path = dest / 'genbank_first_record_header.txt'
            if header_path.exists():
                header = header_path.read_text()
            else:
                with urllib.request.urlopen(gbff_url, timeout=90) as response:
                    with gzip.GzipFile(fileobj=response) as zipped:
                        lines = []
                        for _ in range(4000):
                            line = zipped.readline().decode('utf-8', errors='replace')
                            if not line or line.startswith('FEATURES') or line.startswith('ORIGIN'):
                                break
                            lines.append(line)
                        else:
                            raise ValueError('GenBank header limit reached')
                header = ''.join(lines)
            (dest / 'genbank_first_record_header.txt').write_text(header)
            (dest / 'gff_header.txt').write_text('\n'.join(comments) + '\n')
            record = dict(assembly_accession=accession, status='ok', report_fields=report_fields,
                          gff_structural_feature_sources=sources, gff_header=comments,
                          gff_structural_method_evidence=sorted(method_evidence),
                          genbank_first_record_header=header,
                          report_url=base + '/' + prefix + '_assembly_report.txt', gff_url=gff_url, gbff_url=gbff_url,
                          gff_compressed_sha256=hashlib.sha256(raw_gff).hexdigest(),
                          genbank_header_scope='First record header only, prior to FEATURES; not whole GBFF')
            write_json(result_file, record)
            return record
        except Exception as error:
            last_error = str(error)
            time.sleep(attempt + 1)
    record = dict(assembly_accession=accession, status='error', error=last_error)
    write_json(result_file, record)
    return record


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [r for r in json.loads((INV / 'file_audit.json').read_text())
            if r['all_three'] and r['database'] == 'GenBank']
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
        for result in executor.map(collect, rows):
            results.append(result)
            if len(results) % 25 == 0:
                print('Method sources checked: %d/435; errors: %d' %
                      (len(results), sum(r['status'] != 'ok' for r in results)), flush=True)
    write_json(OUT / 'raw_methods.json', results)
    errors = [r for r in results if r['status'] != 'ok']
    if errors:
        raise RuntimeError('Incomplete method audit; inspect error records and rerun: ' + str(errors))
    print('Collected all 435 method-source records', flush=True)


if __name__ == '__main__':
    main()
