"""Resolve strain/sample duplicates without collapsing unknown strain labels."""
import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import write_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/orthology_20261030'
CACHE = ROOT / 'data/orthology/biosamples'
MISSING = {'', 'na', 'n/a', 'unknown', 'not provided', 'not applicable', 'missing', 'unspecified'}


def label(value):
    value = re.sub(r'\s+', ' ', value or '').strip()
    return '' if value.casefold() in MISSING else value


def normalized(value):
    # Deliberately retain punctuation: BY4741, BY4742, CEN.PK strains and derivatives
    # must not become identical merely through aggressive string normalization.
    return re.sub(r'\s+', '', label(value)).casefold()


def biosamples(ids, fetch):
    CACHE.mkdir(parents=True, exist_ok=True)
    records = {}
    for path in CACHE.glob('batch_*.xml'):
        for node in ET.parse(path).getroot().iter('BioSample'):
            accession = node.attrib['accession']
            attributes = {a.attrib.get('harmonized_name', a.attrib.get('attribute_name', '')): a.text or ''
                          for a in node.findall('./Attributes/Attribute')}
            records[accession] = dict(accession=accession, attributes=attributes,
                                      title=node.findtext('./Description/Title') or '')
    missing = sorted(set(ids) - set(records))
    if missing and not fetch:
        raise ValueError('Missing BioSamples; run with --fetch-biosamples')
    for index in range(0, len(missing), 30):
        batch = missing[index:index + 30]
        term = ' OR '.join(accession + '[Accession]' for accession in batch)
        search_query = urllib.parse.urlencode(dict(db='biosample', term=term, retmax=100, retmode='json', tool='GeneTrace'))
        with urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + search_query, timeout=90) as response:
            uids = json.loads(response.read())['esearchresult']['idlist']
        if not uids:
            raise ValueError('No NCBI BioSample UIDs resolved for accession batch')
        time.sleep(.4)
        query = urllib.parse.urlencode(dict(db='biosample', id=','.join(uids), retmode='xml', tool='GeneTrace'))
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + query
        with urllib.request.urlopen(url, timeout=90) as response:
            body = response.read()
        tree = ET.fromstring(body)
        found = {n.attrib['accession'] for n in tree.iter('BioSample')}
        if not set(batch) <= found:
            raise ValueError('BioSample response incomplete: ' + str(set(batch) - found))
        # Batch name is content-derived to avoid overwriting a previous snapshot.
        import hashlib
        (CACHE / ('batch_' + hashlib.sha256(','.join(batch).encode()).hexdigest()[:12] + '.xml')).write_bytes(body)
        print('BioSamples retrieved: %d/%d' % (min(index + 30, len(missing)), len(missing)), flush=True)
        time.sleep(.4)
    if missing:
        return biosamples(ids, False)
    return records


def rank(row):
    # Never reward a large predicted-gene count or large genome size: both may
    # be artifacts or combined haplotypes. Keep all retained strains regardless QC.
    levels = {'Complete Genome': 4, 'Chromosome': 3, 'Scaffold': 2, 'Contig': 1}
    reference = row.get('paired_refseq') == 'GCF_000146045.2'
    return (reference, row['genome_representation'] == 'Full',
            row.get('number_of_proteins', 0) >= 1000,
            levels.get(row['assembly_level'], 0), -row.get('n_fraction', 1),
            row.get('fasta_sequence_n50_bp', 0), row.get('release_date', ''), row['assembly_accession'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fetch-biosamples', action='store_true')
    args = parser.parse_args()
    rows = json.loads((ROOT / 'outputs/s_cerevisiae_methods_20261030/methods_analysis.json').read_text())['rows']
    samples = biosamples([r['biosample'] for r in rows if label(r.get('biosample'))], args.fetch_biosamples)
    OUT.mkdir(parents=True, exist_ok=True)
    parent = {r['assembly_accession']: r['assembly_accession'] for r in rows}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a,b):
        parent[find(a)] = find(b)
    sample_seen, strain_seen = {}, {}
    enriched = []
    edges = []
    for row in rows:
        row = dict(row)
        accession = row['assembly_accession']
        sample = label(row.get('biosample'))
        attrs = samples.get(sample, {}).get('attributes', {})
        strain = label(row.get('strain'))
        bs_strain = label(attrs.get('strain'))
        isolate = label(attrs.get('isolate'))
        row.update(biosample_strain=bs_strain or 'N/A', biosample_isolate=isolate or 'N/A',
                   biosample_title=samples.get(sample, {}).get('title', 'N/A'),
                   biosample_url='https://www.ncbi.nlm.nih.gov/biosample/' + sample if sample else 'N/A')
        aliases = [label(x) for x in strain.split(';') if label(x)] if strain else []
        row['strain_labels_conflict'] = bool(strain and bs_strain and normalized(bs_strain) not in {normalized(x) for x in aliases})
        row['resolved_strain'] = strain or bs_strain or isolate or 'N/A'
        row['strain_resolution_source'] = ('assembly strain' if strain else 'BioSample strain' if bs_strain
                                           else 'BioSample isolate' if isolate else 'accession fallback')
        if sample:
            if sample in sample_seen:
                union(accession, sample_seen[sample]); edges.append(dict(a=accession,b=sample_seen[sample],reason='same BioSample',value=sample))
            sample_seen[sample] = accession
        # Conflicting strain declarations are not used to bridge distinct samples.
        keys = {normalized(x) for x in (aliases or [row['resolved_strain']])} if not row['strain_labels_conflict'] else set()
        for key in keys:
            if key and key not in {'s.cerevisiae','saccharomycescerevisiae','yeast','brewersyeast'}:
                if key in strain_seen:
                    union(accession, strain_seen[key]); edges.append(dict(a=accession,b=strain_seen[key],reason='exact normalized strain/isolate label or explicit semicolon alias',value=key))
                strain_seen[key] = accession
        enriched.append(row)
    groups = defaultdict(list)
    for row in enriched:
        groups[find(row['assembly_accession'])].append(row)
    selected, decisions = [], []
    for group in groups.values():
        winner = max(group, key=rank)
        selected.append(winner)
        aliases = sorted({r['resolved_strain'] for r in group})
        for row in group:
            decisions.append(dict(assembly_accession=row['assembly_accession'], resolved_strain=row['resolved_strain'],
                                  biosample=row.get('biosample') or 'N/A', selected=row is winner,
                                  representative=winner['assembly_accession'], duplicate_group_size=len(group),
                                  strain_aliases='; '.join(aliases), strain_resolution_source=row['strain_resolution_source'],
                                  strain_labels_conflict=row['strain_labels_conflict'],
                                  representative_policy='reference; Full; not extremely sparse; assembly level; fewer Ns; sequence N50; release date; accession'))
    selected.sort(key=lambda r:r['assembly_accession'])
    write_json(OUT / 'selected_assemblies.json', selected)
    write_json(OUT / 'strain_deduplication_edges.json', edges)
    write_json(OUT / 'strain_deduplication.json', decisions)
    for name, data in [('strain_deduplication.csv', decisions), ('selected_assemblies.csv', selected)]:
        with (OUT / name).open('w',newline='',encoding='utf-8-sig') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    summary = dict(input_assemblies=len(rows), retained_assemblies=len(selected), removed_duplicate_assemblies=len(rows)-len(selected),
                   retained_protein_records=sum(r['number_of_proteins'] for r in selected),
                   retained_without_resolved_strain=sum(r['resolved_strain']=='N/A' for r in selected),
                   retained_sparse_proteomes=[r['assembly_accession'] for r in selected if r['number_of_proteins']<1000],
                   conflicting_strain_metadata=[r['assembly_accession'] for r in enriched if r['strain_labels_conflict']],
                   scope='All distinct resolvable strain/sample units from original cohort; no QC exclusion of retained strains',
                   unresolved_identity='Unknowns are retained by accession/BioSample, never merged as N/A',
                   caveat='Exact strain labels and BioSamples are operational identities, not sequence-verified clone identities; alias-only duplicates can remain')
    write_json(OUT / 'cohort_summary.json', summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
