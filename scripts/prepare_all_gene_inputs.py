"""Fetch every retained assembly and export every available coding locus.

All source protein records remain in source FASTA. For orthogroup copy counts,
select the longest translated isoform per mapped gene; retain unmapped proteins
as explicitly marked separate units. Never collapse distinct loci by sequence.
"""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import read_fasta, read_gff, write_fasta, write_json, sha256

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/orthology'
OUT=ROOT/'outputs/orthology_20261030'
SUFFIXES={'genome.fna.gz':'genomic.fna.gz','annotation.gff.gz':'genomic.gff.gz',
          'cds.fna.gz':'cds_from_genomic.fna.gz','protein.faa.gz':'protein.faa.gz'}


def fetch_assembly(row):
    accession=row['assembly_accession']
    inventory=ROOT/'data/inventory/s_cerevisiae_2026-10-30'
    original=next(r for r in json.loads((inventory/'file_audit.json').read_text()) if r['assembly_accession']==accession)
    urlbase=original['ftp_path'].rstrip('/')
    prefix=urlbase.split('/')[-1]
    checksums={Path(line.split()[-1]).name:line.split()[0] for line in
               (inventory/'file_manifests'/(accession+'.md5.txt')).read_text().splitlines() if line.strip()}
    directory=BASE/'source'/accession
    directory.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for local,suffix in SUFFIXES.items():
        target=directory/local
        filename=prefix+'_'+suffix
        expected=checksums[filename]
        if not target.exists() or hashlib.md5(target.read_bytes()).hexdigest()!=expected:
            if shutil.disk_usage(ROOT).free<8*1024**3:
                raise RuntimeError('Less than 8 GiB free before source download')
            error=None
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(urlbase+'/'+filename,timeout=120) as response:
                        body=response.read()
                    if hashlib.md5(body).hexdigest()!=expected:
                        raise ValueError('NCBI source snapshot changed: '+filename)
                    partial=target.with_suffix(target.suffix+'.partial')
                    partial.write_bytes(body)
                    os.replace(partial,target)
                    error=None
                    break
                except Exception as exc:
                    error=exc
                    time.sleep(attempt+1)
            if error:
                raise error
        manifest.append(dict(file=local,url=urlbase+'/'+filename,md5=expected,sha256=sha256(target)))
    write_json(directory/'source_manifest.json',manifest)
    return accession


def prepare(row):
    accession=row['assembly_accession']
    directory=BASE/'source'/accession
    proteins=read_fasta(directory/'protein.faa.gz')
    if len(proteins)!=row['number_of_proteins']:
        raise ValueError('Protein record count changed: '+accession)
    features=read_gff(directory/'annotation.gff.gz')
    by_id={f['attrs']['ID']:f for f in features if 'ID' in f['attrs']}
    protein_genes={}
    gene_meta={}
    def roots(identifier,visited=None):
        visited=set() if visited is None else visited
        if identifier in visited:
            raise ValueError('Cyclic GFF parent chain')
        visited=visited|{identifier}
        feature=by_id.get(identifier)
        if feature is None:
            return set()
        if feature['type']=='gene':
            gene_meta[identifier]=feature
            return {identifier}
        result=set()
        for parent in feature['attrs'].get('Parent','').split(','):
            if parent:
                result.update(roots(parent,visited))
        return result
    for feature in features:
        if feature['type']!='CDS':
            continue
        attrs=feature['attrs']
        protein=attrs.get('protein_id')
        if not protein and attrs.get('Name') in proteins:
            protein=attrs['Name']
        if not protein:
            continue
        genes=set()
        for parent in attrs.get('Parent','').split(','):
            if parent:
                genes.update(roots(parent))
        if not genes:
            # Gene-less submissions still contain coding loci. Group by stable CDS
            # ID, not peptide identity, and keep this mapping distinction explicit.
            cds_id=attrs.get('ID') or protein
            gene_meta.setdefault(cds_id,feature)
            genes={cds_id}
        protein_genes.setdefault(protein,set()).update(genes)
    candidates={}
    for protein in proteins:
        for gene in protein_genes.get(protein,{'unmapped:'+protein}):
            candidates.setdefault(gene,[]).append(protein)
    fasta={}
    mapping=[]
    for index,(gene,isoforms) in enumerate(sorted(candidates.items()),1):
        chosen=max(isoforms,key=lambda p:(len(proteins[p].rstrip('*')),p))
        sequence=proteins[chosen].rstrip('*')
        if not sequence:
            raise ValueError('Empty representative protein after terminal stop removal')
        token=accession.replace('.','_')+'__g%07d'%index
        fasta[token]=sequence
        feature=gene_meta.get(gene,{})
        attrs=feature.get('attrs',{})
        mapping.append(dict(sequence_id=token,assembly_accession=accession,gene_id=gene,
                            ncbi_protein_id=chosen,all_isoform_protein_ids=';'.join(sorted(isoforms)),
                            mapping_type='unmapped_protein' if gene.startswith('unmapped:') else
                                         'gene' if feature.get('type')=='gene' else 'CDS_without_gene_feature',
                            gene_name=attrs.get('Name',''),locus_tag=attrs.get('locus_tag',''),
                            seqid=feature.get('seqid',''),start=feature.get('start',''),end=feature.get('end',''),
                            strand=feature.get('strand',''),protein_length=len(sequence),
                            internal_stop='*' in sequence,terminal_stop_removed=proteins[chosen]!=sequence,
                            isoform_policy='longest protein per coding locus; original records retained'))
    protein_dir=BASE/'proteomes'
    protein_dir.mkdir(parents=True,exist_ok=True)
    write_fasta(protein_dir/(accession+'.faa'),fasta)
    map_dir=BASE/'gene_maps'
    map_dir.mkdir(exist_ok=True)
    with (map_dir/(accession+'.tsv')).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(mapping[0]),delimiter='\t');writer.writeheader();writer.writerows(mapping)
    annotated_genes={f['attrs'].get('ID') for f in features if f['type']=='gene' and
                     f['attrs'].get('gene_biotype')=='protein_coding'}
    absent=sorted(annotated_genes-set(candidates))
    write_json(map_dir/(accession+'.genes_without_protein.json'),absent)
    summary=dict(assembly_accession=accession,resolved_strain=row['resolved_strain'],
                 source_protein_records=len(proteins),orthology_coding_units=len(fasta),
                 unmapped_protein_units=sum(r['mapping_type']=='unmapped_protein' for r in mapping),
                 coding_genes_without_protein=len(absent),internal_stop_units=sum(r['internal_stop'] for r in mapping),
                 proteome_sha256=sha256(protein_dir/(accession+'.faa')),
                 sparse_annotation_warning=len(proteins)<1000,
                 input_scope='All coding loci with available protein records; no locus subsampling')
    write_json(map_dir/(accession+'.summary.json'),summary)
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--download',action='store_true')
    args=parser.parse_args()
    rows=json.loads((OUT/'selected_assemblies.json').read_text())
    if args.download:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            for index,_ in enumerate(executor.map(fetch_assembly,rows),1):
                if index%10==0:
                    print('Assembly source packages ready: %d/%d'%(index,len(rows)),flush=True)
    summary=[]
    for index,row in enumerate(rows,1):
        summary.append(prepare(row))
        if index%10==0:
            print('All-gene proteomes prepared: %d/%d'%(index,len(rows)),flush=True)
    write_json(OUT/'all_gene_input_summary.json',summary)
    with (OUT/'all_gene_input_summary.csv').open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    print(json.dumps(dict(assemblies=len(summary),source_protein_records=sum(r['source_protein_records'] for r in summary),
                          orthology_coding_units=sum(r['orthology_coding_units'] for r in summary),
                          unmapped_units=sum(r['unmapped_protein_units'] for r in summary),
                          genes_without_protein=sum(r['coding_genes_without_protein'] for r in summary)),indent=2),flush=True)


if __name__=='__main__':
    main()
