#!/usr/bin/env python3
import csv,json,pathlib,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff
rows=json.loads((ROOT/'data/bacteria/cohort.json').read_text());out=ROOT/'data/bacteria/gene_maps';out.mkdir(exist_ok=True)
for row in rows:
 acc=row['assembly_accession'];source=ROOT/row['source_directory'];proteins=read_fasta(source/'protein.faa.gz');features=read_gff(source/'annotation.gff.gz');parents={f['attrs'].get('ID'):f for f in features if f['type']=='gene'};mapped={}
 for f in features:
  if f['type']!='CDS':continue
  pid=f['attrs'].get('protein_id')
  if pid not in proteins:continue
  parent=f['attrs'].get('Parent','').split(',')[0];attrs=dict(parents.get(parent,{}).get('attrs',{}));attrs.update(f['attrs']);locus=attrs.get('locus_tag',parent)
  mapped.setdefault(pid,dict(sequence_id=pid,assembly_accession=acc,gene_id=parent or f['attrs'].get('ID',pid),ncbi_protein_id=pid,gene_name=attrs.get('gene',locus),locus_tag=locus,seqid=f['seqid'],start=f['start'],end=f['end'],strand=f['strand'],protein_length=len(proteins[pid])))
 if set(proteins)-set(mapped):raise ValueError('Unmapped bacterial protein IDs: '+str(list(set(proteins)-set(mapped))[:5]))
 with (out/(acc+'.tsv')).open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(next(iter(mapped.values()))),delimiter='\t');writer.writeheader();writer.writerows(mapped.values())
 row.update(resolved_strain=row['strain'],reference_accession='GCF_000005845.2')
 for key in ['sequencing_group','raw_read_status','biosample']:row.setdefault(key,'N/A: not surveyed')
 print(acc,len(mapped),'coding proteins mapped')
(ROOT/'data/bacteria/cohort.json').write_text(json.dumps(rows,indent=2)+'\n')
