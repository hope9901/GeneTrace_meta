#!/usr/bin/env python3
"""Proper core+assign OrthoFinder for all bacterial coding loci, disk-guarded."""
import collections,csv,json,os,pathlib,shutil,subprocess,sys,datetime
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from genetrace.io import read_fasta,write_fasta,sha256
import run_scalable_orthology as workflow
OUT=ROOT/'outputs/ecoli100_orthology_20261102';DATA=ROOT/'data/ecoli100_orthology_20261102';CORE=ROOT/'runs/ecoli100_orthofinder_core_20261102'
def dump(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True);(DATA/'core').mkdir(exist_ok=True);(DATA/'additional').mkdir(exist_ok=True);workflow.OUT=OUT;workflow.STATUS=OUT/'status.json'
 metadata=json.loads((ROOT/'data/bacteria/cohort.json').read_text());accs={r['assembly_accession'] for r in metadata};statefile=OUT/'state.json'
 if statefile.exists():state=json.loads(statefile.read_text())
 else:
  from Bio import Phylo
  tree=Phylo.read(ROOT/'outputs/ecoli100_context_20261102/core30.tree','newick');selected=['GCF_000005845.2']
  while len(selected)<16:selected.append(max(sorted(accs-set(selected)),key=lambda a:min(tree.distance(a,b) for b in selected)))
  manifest=[];expected={}
  for row in metadata:
   acc=row['assembly_accession'];proteins=read_fasta(ROOT/row['source_directory']/'protein.faa.gz');groups=collections.defaultdict(list)
   for g in csv.DictReader((ROOT/'data/bacteria/gene_maps'/f'{acc}.tsv').open(),delimiter='\t'):groups[g['gene_id']].append(g['ncbi_protein_id'])
   records={acc+'__'+max(pids,key=lambda p:(len(proteins[p]),p)):proteins[max(pids,key=lambda p:(len(proteins[p]),p))] for pids in groups.values()};p=DATA/('core' if acc in selected else 'additional')/(acc+'.faa');write_fasta(p,records);expected[acc]=sorted(records);manifest.append(dict(assembly=acc,source_records=len(proteins),coding_loci=len(records),longest_isoform_per_gene=True,path=str(p),sha256=sha256(p)))
  dump(OUT/'inputs.json',manifest);dump(OUT/'expected_genes.json',expected);state=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),retained_assemblies=100,coding_units=sum(len(x) for x in expected.values()),core_size=16,additional_size=84,core=selected,additional=sorted(accs-set(selected)),commands=[],core_selection='Farthest-first distances in provisional 30-protein tree; gene input not subsampled')
  dump(statefile,state)
 env=dict(os.environ);env.update(PATH=str(ROOT/'tools/orthofinder32-env/bin')+os.pathsep+str(ROOT/'tools/orthofinder-env/bin')+os.pathsep+env.get('PATH',''),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');base=[str(ROOT/'tools/orthofinder32-env/bin/python'),str(ROOT/'tools/orthofinder32-env/bin/orthofinder')];opts=['-X','-S','diamond','-t','2','-a','1']
 try:
  if 'core_results' not in state:
   if CORE.exists():
    dirs=sorted(CORE.glob('Results_*'))
    if not dirs:raise RuntimeError('Core results directory missing')
    cmd=base+['-b',str(dirs[-1])]+opts
   else:cmd=base+['-f',str(DATA/'core'),'-o',str(CORE)]+opts
   state['commands'].append(cmd);dump(statefile,state);workflow.run_stage('core_inference',cmd,env,state,CORE,2);state['core_results']=str(workflow.final_result(CORE,set(state['core']),require_tree=True));dump(statefile,state)
  core=pathlib.Path(state['core_results'])
  if 'full_results' not in state:
   dirs=[p for p in core.parent.glob('Results_*') if p!=core]
   cmd=base+(['-b',str(sorted(dirs)[-1])] if dirs else ['--assign',str(DATA/'additional'),'--core',str(core)])+opts;state['commands'].append(cmd);dump(statefile,state);workflow.run_stage('all_genome_assignment',cmd,env,state,core.parent,2);state['full_results']=str(workflow.final_result(core.parent,accs));dump(statefile,state)
  result=pathlib.Path(state['full_results']);expected=json.loads((OUT/'expected_genes.json').read_text());observed=collections.defaultdict(list)
  for filename in ['Orthogroups.tsv','Orthogroups_UnassignedGenes.tsv']:
   with (result/'Orthogroups'/filename).open() as f:
    reader=csv.DictReader(f,delimiter='\t')
    if set(reader.fieldnames[1:])!=accs:raise ValueError('Final cohort mismatch')
    for row in reader:
     for acc in accs:observed[acc].extend(x.strip() for x in row[acc].split(',') if x.strip())
  for acc in accs:
   if collections.Counter(observed[acc])!=collections.Counter(expected[acc]):raise ValueError('Gene accounting mismatch: '+acc)
  dump(OUT/'status.json',dict(status='complete',retained_assemblies=100,coding_units=state['coding_units'],results=str(result),all_gene_accounting_verified=True))
 except Exception as e:dump(OUT/'status.json',dict(status='failed_resumable',error=str(e),state=str(statefile)));raise
if __name__=='__main__':main()
