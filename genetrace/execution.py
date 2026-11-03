"""Bounded, real local reanalysis of a fixed public-data case. No model-supplied commands."""
import collections,datetime,gzip,json,pathlib,shutil,subprocess,time
from .agent import ROOT,model_config,request_model
from .io import read_fasta,read_gff,sha256
from .regional_qc import analyze_sam
from .read_provenance import assess_read_provenance
PACK=ROOT/'demo/execution/nimt'
LABELS={'inspect_annotation':'Inspect current structural annotation','align_candidate_reads':'Align selected reads and check nimT','check_neighbor_regions':'Test neighboring genes as local controls','inspect_read_provenance':'Audit read structure and source linkage'}

def binary(name):
 p=ROOT/'tools/readqc-env/bin'/name
 found=str(p) if p.is_file() else shutil.which(name)
 if not found:raise ValueError(name+' is required for execution mode')
 return found

def availability():
 try:
  binary('minimap2');binary('samtools')
  return (PACK/'manifest.json').is_file()
 except ValueError:return False

def verify_pack():
 m=json.loads((PACK/'manifest.json').read_text())
 for name,entry in m['files'].items():
  if name not in ['genome.fna.gz','annotation.gff.gz','reads_1.fastq.gz','reads_2.fastq.gz']:raise ValueError('Unexpected execution input')
  if sha256(PACK/name)!=entry['sha256']:raise ValueError('Execution input checksum mismatch')
 return m

def compact(r):
 return {k:r[k] for k in ['name','region','signals','metrics','partial_input','coherent_substitutions','indels','causal_diagnosis','limitations']}

def verdict(evidence):
 by={e['tool']:e for e in evidence};cand=by.get('align_candidate_reads');neighbors=by.get('check_neighbor_regions');prov=by.get('inspect_read_provenance')
 ids=[e['id'] for e in evidence]
 if not cand:return dict(status='incomplete',observation='No new read alignment result yet.',possible_causes=[],limitations=['Required checks are incomplete.'],next_checks=['Run candidate read alignment.'],evidence_ids=ids)
 c=cand['result']['candidate'];signal=any(x.startswith('coherent_') for x in c['signals'])
 broader=bool(neighbors and any(any(s.startswith('coherent_') for s in r['signals']) for r in neighbors['result']['neighbor_regions']))
 complete=bool(neighbors and prov) if signal else bool(prov)
 return dict(status='unresolved' if complete else 'provisional',observation=('Candidate and neighboring genes disagree with selected reads; a local-only explanation is insufficient.' if broader else 'Candidate read disagreement requires context.' if signal else 'No coherent disagreement flag under the current selection and thresholds.'),possible_causes=['Biological sequence differences','Assembly or systematic read/alignment effects','Read dataset / assembly-version mismatch'],limitations=['Reference-selected read subset; depth is not an unbiased whole-run estimate.','Exact run-to-assembly-version linkage unverified.','No independent physical read validation or confirmed biological loss.']+([] if complete else ['Required follow-up checks are incomplete.']),next_checks=['Verify original run and assembly version records.','Confirm with independently generated, correctly linked reads.','Reassess the locus with independent assembly/alignment evidence.'],evidence_ids=ids,local_only_explanation_supported=False if broader else None,causal_diagnosis='unresolved',biological_loss_inferred=False)

class Executor:
 def __init__(self,folder,notify=lambda event:None):
  self.folder=pathlib.Path(folder);self.folder.mkdir(parents=True,exist_ok=True);self.notify=notify;self.manifest=verify_pack();self.genome=read_fasta(PACK/'genome.fna.gz');self.regions=self.manifest['regions'];self.commands=[]
 def run(self,args,output=None):
  start=time.monotonic();rec=dict(command=[str(x).replace(str(ROOT)+'/','') for x in args],started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
  if output:
   with pathlib.Path(output).open('w') as f:r=subprocess.run([str(a) for a in args],stdout=f,stderr=subprocess.PIPE,text=True,timeout=120)
   result=''
  else:r=subprocess.run([str(a) for a in args],capture_output=True,text=True,timeout=120);result=r.stdout
  rec.update(returncode=r.returncode,elapsed_seconds=round(time.monotonic()-start,3))
  if '--version' in args:rec['version_output']=result[:1500]
  self.commands.append(rec)
  if r.returncode:raise RuntimeError('Allowlisted analysis command failed')
  return result
 def analyze(self,region):
  sam=self.run([binary('samtools'),'view',self.folder/'aligned.bam',f"{region['seqid']}:{max(1,region['start']-1000)}-{region['end']+1000}"])
  r=analyze_sam(sam,self.genome[region['seqid']],region['seqid'],region['start'],region['end'],partial_input=True);r['name']=region['name'];return compact(r)
 def execute(self,tool):
  if tool not in LABELS:raise ValueError('Tool not allowed')
  self.notify(dict(kind='tool_started',tool=tool,label=LABELS[tool]));start=time.monotonic();self.commands=[]
  if tool=='inspect_annotation':
   target=self.regions[0];features=read_gff(PACK/'annotation.gff.gz');overlap=[{k:f[k] for k in ['seqid','type','start','end','strand','attrs']} for f in features if f['seqid']==target['seqid'] and f['type'] in ['gene','pseudogene','CDS'] and f['start']<=target['end'] and f['end']>=target['start']]
   data=dict(target=target,overlapping_features=overlap,recorded_disruption=any('frameshift' in json.dumps(f).lower() or 'pseudo' in json.dumps(f).lower() for f in overlap),operation='Fresh GFF parsing; no annotation prediction performed')
  elif tool=='align_candidate_reads':
   genome=self.folder/'reference.fna'
   from .io import write_fasta
   write_fasta(genome,self.genome)
   self.run([binary('minimap2'),'--version']);self.run([binary('samtools'),'--version'])
   self.run([binary('minimap2'),'-ax','sr','-t','2',genome,PACK/'reads_1.fastq.gz',PACK/'reads_2.fastq.gz'],self.folder/'alignment.sam')
   self.run([binary('samtools'),'sort','-@','1','-o',self.folder/'aligned.bam',self.folder/'alignment.sam']);self.run([binary('samtools'),'index',self.folder/'aligned.bam'])
   data=dict(candidate=self.analyze(self.regions[0]),read_pairs=self.manifest['source_fastq'][0]['selected_records'],selection=self.manifest['selection'],fresh_alignment=True,reference_scope='Whole assembly; reference-selected read subset',bam_sha256=sha256(self.folder/'aligned.bam'))
  elif tool=='check_neighbor_regions':
   if not (self.folder/'aligned.bam').exists():raise ValueError('Candidate alignment required')
   data=dict(neighbor_regions=[self.analyze(r) for r in self.regions[1:]],fresh_region_analysis=True,shared_alignment=True)
  else:
   if not (self.folder/'aligned.bam').exists():raise ValueError('Candidate alignment required')
   qualities=collections.Counter()
   for mate in [1,2]:
    with gzip.open(PACK/f'reads_{mate}.fastq.gz','rt') as f:
     for i,line in enumerate(f):
      if i%4==3:qualities.update(line.strip())
   sam=self.run([binary('samtools'),'view',self.folder/'aligned.bam'])
   p=assess_read_provenance(sam,qualities);p.pop('note',None)
   data=dict(run=self.manifest['run'],assembly=self.manifest['assembly'],read_platform=self.manifest['read_platform'],assembly_platform=self.manifest['assembly_platform'],exact_assembly_run_use_verified=False,read_provenance=p,partial_input=True,interpretation='Fresh subset audit; a negative subset flag cannot clear whole-run provenance concerns.')
  result=dict(result=data,acquisition='fresh_local_execution',execution=dict(tool=tool,label=LABELS[tool],elapsed_seconds=round(time.monotonic()-start,3),commands=self.commands,input_manifest_sha256=sha256(PACK/'manifest.json'),python_operation=tool,source_code_sha256=sha256(pathlib.Path(__file__)),input_sha256s={n:v['sha256'] for n,v in self.manifest['files'].items()}))
  dest=self.folder/(tool+'.json');dest.write_text(json.dumps(result,indent=2)+'\n');result['sources']=[dict(path=str(dest.relative_to(ROOT)),sha256=sha256(dest))];self.notify(dict(kind='tool_completed',tool=tool,label=LABELS[tool],elapsed_seconds=result['execution']['elapsed_seconds']));return result

def available_tools(done,evidence):
 if 'inspect_annotation' not in done:return ['inspect_annotation']
 if 'align_candidate_reads' not in done:return ['align_candidate_reads']
 tools=[]
 if 'inspect_read_provenance' not in done:tools.append('inspect_read_provenance')
 c=next(e['result']['candidate'] for e in evidence if e['tool']=='align_candidate_reads')
 if any(s.startswith('coherent_') for s in c['signals']) and 'check_neighbor_regions' not in done:tools.append('check_neighbor_regions')
 return tools

def investigate(report,case,persist):
 start=time.monotonic();cfg=model_config();cfg['response_language']='en';report.update(status='running',pipeline='execute',model_id=cfg['model'] if report['mode']=='meta' else None,analysis_execution='Fresh local execution on a reference-selected public read subset; population context is a saved analysis.')
 def notify(event):report['trace'].append(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**event));persist(report)
 executor=Executor(ROOT/'runs/execution'/report['id'],notify)
 population=next(t for t in case['tools'] if t['tool']=='lookup_cohort_comparison')
 report['evidence'].append(dict(id='E001',tool=population['tool'],result=population['data'],sources=population['sources'],acquisition='lookup_completed_analysis_snapshot'))
 done=set();persist(report)
 while True:
  available=available_tools(done,report['evidence'])
  if not available:break
  context=dict(question='Investigate nimT reliability in GCF_001558995.2. Select and execute the next check, then reassess.',execution_mode='Real allowlisted tools, fixed inputs. A reference-selected partial read subset cannot establish unbiased depth, absence, or whole-run quality.',available_tool_descriptions={t:LABELS[t] for t in available},current_judgment=verdict(report['evidence']))
  if report['mode']=='meta':
   step,usage=request_model(context,available,report['evidence'],cfg);report['model_calls']+=1;report['hypotheses']=step['hypotheses'];notify(dict(kind='live_meta_planner',step=step,usage=usage));choice=step['next_tool']
   if choice is None:choice=available[0];notify(dict(kind='guard',reason='Required evidence is incomplete; enforce the fixed minimum investigation policy.',tool=choice))
  else:choice=available[0];notify(dict(kind='baseline_fixed_order',tool=choice,reason='Fixed minimum-check policy; no model'))
  result=executor.execute(choice);report['evidence'].append(dict(id=f"E{len(report['evidence'])+1:03d}",tool=choice,**result));done.add(choice)
  report['judgment']=verdict(report['evidence']);notify(dict(kind='judgment_updated',status=report['judgment']['status'],observation=report['judgment']['observation']))
 if report['mode']=='meta':
  step,usage=request_model(dict(question=case['question'],instruction='All required fresh checks completed. Explain observation, possible causes, limitations, and next verification. next_tool must be null. Neighbor results must inform the conclusion; no causal diagnosis or accuracy claim.',guarded_judgment=report['judgment']),[],report['evidence'],cfg);report['model_calls']+=1;report['hypotheses']=step['hypotheses'];notify(dict(kind='live_meta_conclusion',step=step,usage=usage))
 report.update(status='complete',next_check='Verify source run and assembly version; obtain independent linked reads.',causal_diagnosis='unresolved',biological_loss_inferred=False,elapsed_seconds=round(time.monotonic()-start,3));persist(report)
