#!/usr/bin/env python3
"""Provisional core-protein tree and near-neighbor context; recombination unmodeled."""
import collections,csv,hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,write_fasta,sha256
OUT=ROOT/'outputs/ecoli100_context_20261102';RUN=ROOT/'runs/ecoli100_context_20261102';GROUP=ROOT/'outputs/ecoli100_gene_groups_20261102'
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True);dump(OUT/'status.json',dict(stage='aligning_core_markers'))
 try:
  records=read_fasta(ROOT/'runs/ecoli100_gene_groups_20261102/all_proteins.faa');core=json.loads((GROUP/'single_copy_core.json').read_text());eligible=[g for g in core if all(100<=len(records[m])<=500 for m in g['members']) and len({records[m] for m in g['members']})>1];selected=sorted(eligible,key=lambda g:hashlib.sha256(g['group'].encode()).hexdigest())[:30]
  if len(selected)<10:raise ValueError('Insufficient variable core markers')
  accs=sorted(m.split('|')[0] for m in selected[0]['members']);combined={a:'' for a in accs};manifest=[]
  for i,g in enumerate(selected,1):
   fasta=RUN/(g['group']+'.faa');aln=RUN/(g['group']+'.aln.faa');write_fasta(fasta,{m.split('|')[0]:records[m] for m in g['members']});cmd=[str(ROOT/'tools/orthofinder-env/bin/mafft'),'--auto','--thread','1',str(fasta)]
   with aln.open('w') as output,(RUN/(g['group']+'.stderr')).open('w') as err:subprocess.run(cmd,stdout=output,stderr=err,check=True,timeout=180)
   aligned=read_fasta(aln);length=len(next(iter(aligned.values())));keep=[p for p in range(length) if sum(aligned[a][p] in 'ACDEFGHIKLMNPQRSTVWY' for a in accs)>=95]
   for a in accs:combined[a]+=''.join(aligned[a][p] for p in keep)
   manifest.append(dict(group=g['group'],reference_gene=g.get('reference_gene'),retained_columns=len(keep),command=cmd,input_sha256=sha256(fasta)));dump(OUT/'status.json',dict(stage='aligning_core_markers',completed=i,total=len(selected)));print(i,'/',len(selected),flush=True)
  concat=RUN/'concatenated.faa';write_fasta(concat,combined);treefile=OUT/'core30.tree';cmd=[str(ROOT/'tools/orthofinder-env/bin/FastTree'),'-wag','-gamma',str(concat)]
  with treefile.open('w') as out,(RUN/'fasttree.stderr').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=300)
  from Bio import Phylo
  tree=Phylo.read(treefile,'newick');metadata={r['assembly_accession']:r for r in json.loads((ROOT/'data/bacteria/cohort.json').read_text())};neighbors={}
  for a in accs:
   others=sorted([(b,tree.distance(a,b)) for b in accs if b!=a],key=lambda x:(x[1],x[0]));neighbors[a]=[dict(assembly=b,distance=d,bioproject=metadata[b].get('bioproject','N/A'),sequencing_group=metadata[b].get('sequencing_group','N/A')) for b,d in others[:5]]
  dump(OUT/'nearest_neighbors.json',neighbors);dump(OUT/'manifest.json',dict(markers=manifest,tree_command=cmd,alignment_sha256=sha256(concat),tree_sha256=sha256(treefile),selection='30 variable single-copy homolog clusters, lengths 100–500 AA; deterministic SHA256 group ordering',limitations=['Provisional 30-protein tree, not a recombination-corrected genome phylogeny.','Homolog clusters have not undergone gene-tree orthology reconciliation.','Near ties/identical marker sequences cannot resolve close strains.','Tree context is not independent validation of groups derived from the same proteins.']))
  dump(OUT/'status.json',dict(stage='complete',assemblies=len(accs),markers=len(selected),alignment_columns=len(next(iter(combined.values()))),phylogeny='provisional_core30_protein_tree'))
 except Exception as e:dump(OUT/'status.json',dict(stage='failed',error=str(e)));raise
if __name__=='__main__':main()
