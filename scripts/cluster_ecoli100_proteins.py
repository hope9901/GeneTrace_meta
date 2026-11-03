#!/usr/bin/env python3
"""All-protein homolog screening; not a phylogenetically inferred orthogroup call."""
import collections,csv,datetime,json,pathlib,shutil,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,write_fasta,sha256
OUT=ROOT/'outputs/ecoli100_gene_groups_20261102';RUN=ROOT/'runs/ecoli100_gene_groups_20261102'
def dump(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def table(p,rows,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True);cohort=json.loads((ROOT/'data/bacteria/cohort.json').read_text());records={};maps={}
 def status(stage,**kw):dump(OUT/'status.json',dict(stage=stage,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw))
 try:
  status('preparing',assemblies=len(cohort));inputs=[]
  for r in cohort:
   acc=r['assembly_accession'];p=ROOT/r['source_directory']/'protein.faa.gz';inputs.append(dict(assembly=acc,sha256=sha256(p)))
   with (ROOT/'data/bacteria/gene_maps'/f'{acc}.tsv').open() as f:gm={g['ncbi_protein_id']:g for g in csv.DictReader(f,delimiter='\t')}
   for pid,seq in read_fasta(p).items():
    key=acc+'|'+pid;records[key]=seq;maps[key]=gm[pid]
  fasta=RUN/'all_proteins.faa';write_fasta(fasta,records);binary=ROOT/'tools/orthofinder-env/bin/mmseqs';cmd=[str(binary),'easy-linclust',str(fasta),str(RUN/'groups'),str(RUN/'tmp'),'--min-seq-id','0.9','-c','0.8','--cov-mode','0','--threads','2','--split-memory-limit','2G','--remove-tmp-files','1']
  signature=dict(inputs=inputs,fasta_sha256=sha256(fasta),command=cmd,version=subprocess.check_output([str(binary),'version'],text=True).strip(),note='90% AA identity, 80% bidirectional coverage homolog clusters; paralogy/phylogeny unresolved; not final orthogroups')
  dump(OUT/'manifest.json',signature);status('clustering',assemblies=len(cohort),proteins=len(records),threads=2)
  if shutil.disk_usage(ROOT).free<7*1024**3:raise RuntimeError('7 GiB reserve required before clustering')
  subprocess.run(cmd,stdout=(RUN/'mmseqs.stdout').open('w'),stderr=(RUN/'mmseqs.stderr').open('w'),check=True)
  groups=collections.defaultdict(list)
  with (RUN/'groups_cluster.tsv').open() as f:
   for line in f:
    rep,member=line.rstrip().split('\t');groups[rep].append(member)
  assigned=[m for members in groups.values() for m in members]
  if len(assigned)!=len(records) or len(set(assigned))!=len(records) or set(assigned)!=set(records):raise ValueError('Protein accounting failed')
  accs=[r['assembly_accession'] for r in cohort];summary=[];members_out=[];candidates=[];core=[]
  for index,(rep,members) in enumerate(sorted(groups.items()),1):
   gid=f'HG{index:06d}';counts=collections.Counter(m.split('|')[0] for m in members);refmembers=[m for m in members if m.startswith('GCF_000005845.2|')]
   summary.append(dict(group=gid,representative=rep,protein_count=len(members),assembly_count=len(counts),single_copy=all(n==1 for n in counts.values()),reference_members=';'.join(refmembers),**{acc:counts[acc] for acc in accs}))
   for member in members:members_out.append(dict(group=gid,member=member,**maps[member]))
   if len(counts)==100 and all(n==1 for n in counts.values()):core.append(dict(group=gid,members=members,reference_gene=maps[refmembers[0]] if len(refmembers)==1 else None))
   if len(refmembers)==1 and 95<=len(counts)<100 and all(n==1 for n in counts.values()):
    for acc in accs:
     if acc not in counts:candidates.append(dict(group=gid,assembly_accession=acc,reference_member=refmembers[0],reference_gene=maps[refmembers[0]]['gene_name'],reference_locus=maps[refmembers[0]]['locus_tag'],group_prevalence=len(counts),status='protein_cluster_missing_not_gene_loss',biological_loss_inferred=False,error_inferred=False))
  table(OUT/'group_counts.csv',summary,list(summary[0]));table(OUT/'group_members.csv',members_out,list(members_out[0]));table(OUT/'missing_candidates.csv',candidates,['group','assembly_accession','reference_member','reference_gene','reference_locus','group_prevalence','status','biological_loss_inferred','error_inferred']);dump(OUT/'single_copy_core.json',core)
  dump(OUT/'summary.json',dict(assemblies=len(cohort),proteins=len(records),accounted_proteins=len(assigned),homolog_groups=len(groups),single_copy_core_groups=len(core),missing_protein_cluster_candidates=len(candidates),orthogroup_inference_complete=False,phylogeny_complete=False,biological_loss_inferred=False,limitations=['Identity/coverage clusters are homolog groups, not proven orthogroups.','Fragmented genes, divergent alleles, paralogs and clustering thresholds can create apparent absences.','A missing protein cluster requires genome search, synteny and raw read checks.']))
  status('complete',assemblies=len(cohort),proteins=len(records),homolog_groups=len(groups),single_copy_core_groups=len(core),missing_candidates=len(candidates),orthogroup_inference_complete=False)
 except Exception as e:status('failed_resumable',error=str(e));raise
if __name__=='__main__':main()
