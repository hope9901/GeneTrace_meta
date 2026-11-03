#!/usr/bin/env python3
"""Sequential low-resource proteome-only domain scan for the active bacterial panel."""
import collections,csv,datetime,json,pathlib,shutil,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.domains import scan
from genetrace.io import read_fasta,write_fasta,sha256
OUT=ROOT/'outputs/ecoli100_domains_20261102';RUN=ROOT/'runs/ecoli100_domains_20261102'
def atomic(p,obj):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(p)
def table(p,rows,fields):
 t=p.with_suffix('.tmp')
 with t.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True)
 cohort=json.loads((ROOT/'data/bacteria/cohort.json').read_text());models=json.loads((ROOT/'data/domain_models/manifest.json').read_text());counts=[];members=[];done=0
 def status(stage,**kw):atomic(OUT/'status.json',dict(stage=stage,completed=done,total=len(cohort),cpu_threads=1,scope='all annotated proteins; six selected Pfam models; no genomic ORF scan or raw reads',updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw))
 try:
  for row in cohort:
   if shutil.disk_usage(ROOT).free<6*1024**3:status('storage_required');return
   acc=row['assembly_accession'];status('scanning',current=acc);d=RUN/acc;d.mkdir(exist_ok=True);source=ROOT/row['source_directory']/'protein.faa.gz';proteins=read_fasta(source);fasta=d/'proteins.faa';source_hash=sha256(source);signature=d/'protein_source.json'
   if not fasta.exists() or not signature.exists() or json.loads(signature.read_text()).get('sha256')!=source_hash:write_fasta(fasta,proteins);atomic(signature,dict(sha256=source_hash))
   hits=scan(fasta,d/'proteome')
   for h in hits:members.append(dict(assembly_accession=acc,strain=row['strain'],species='Escherichia coli',**h))
   for model in models:
    profile=model['accession'].split('.')[0];hh=[h for h in hits if h['profile']==profile]
    counts.append(dict(assembly_accession=acc,strain=row['strain'],species='Escherichia coli',profile=profile,proteins_with_domain=len({h['protein_id'] for h in hh}),domain_instances=len(hh),near_complete_instances=sum(h['coverage_class']=='near_complete' for h in hh),partial_instances=sum(h['coverage_class']=='partial' for h in hh),proteome_records=len(proteins),annotation_scope_warning=False,orthology_inferred=False))
   table(OUT/'domain_counts.csv',counts,list(counts[0]));
   if members:table(OUT/'domain_members.csv',members,list(members[0]))
   done+=1;status('scanning',current=acc);print(done,'/',len(cohort),acc,len(hits),flush=True)
  atomic(OUT/'summary.json',dict(assemblies=done,proteome_records=sum(r['number_of_proteins'] for r in cohort),profiles=models,protein_domain_counts=counts,genomic_scan_performed=False,raw_reads_examined=False,orthology_inferred=False,biological_loss_inferred=False))
  lines=['# E. coli 100 proteomes: selected Pfam domain panel','',f'All {done} proteomes were scanned with six fixed Pfam profiles using HMMER --cut_ga. Domain instances and proteins carrying domains are reported separately. Coverage >=70% marks near-complete instances; shorter GA-passing hits remain.','','| Profile | Proteins carrying domain per assembly (min–max) | Domain instances (min–max) |','|---|---|---|']
  for model in models:
   pf=model['accession'].split('.')[0];rr=[r for r in counts if r['profile']==pf];a=[r['proteins_with_domain'] for r in rr];b=[r['domain_instances'] for r in rr];lines.append(f'| {pf} / {model["name"]} | {min(a)}–{max(a)} | {min(b)}–{max(b)} |')
  lines+=['','Counts are descriptive; shared domains do not establish orthology, copy-number variation, function, biological loss or sequencing error. Only six profiles were tested, not the whole Pfam database. No 100-genome six-frame ORF search or bacterial raw read validation was performed. Commands, HMMER version and input hashes are in per-assembly scan manifests.']
  (OUT/'report.md').write_text('\n'.join(lines)+'\n');status('complete')
 except Exception as e:status('failed_resumable',error=str(e));raise
if __name__=='__main__':main()
