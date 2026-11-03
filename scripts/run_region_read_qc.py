#!/usr/bin/env python3
"""Evidence-only HN1 YJU3 regional mapping; partial reads never support absence."""
import argparse,collections,csv,datetime,gzip,hashlib,json,pathlib,re,shutil,subprocess,sys,time,os
os.environ.setdefault("MPLCONFIGDIR","/tmp/genetrace-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME","/tmp/genetrace-cache")
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff,write_fasta,sha256
from genetrace.read_evidence import quality_encoding_assessment,regional_evidence_assessment
ACC='GCA_903819125.2';SEQ='LR813493.2';START=287238;END=288176

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--reads',type=pathlib.Path,required=True);ap.add_argument('--output',type=pathlib.Path,required=True);ap.add_argument('--partial',action='store_true');a=ap.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);status=out/'status.json';commands=[]
 def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
 def stage(s):save('status.json',dict(stage=s,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),partial_input=a.partial,finished=False))
 def run(args,name):
  commands.append(args)
  with (out/(name+'.stdout')).open('w') as so,(out/(name+'.stderr')).open('w') as se:
   subprocess.run(args,stdout=so,stderr=se,check=True)
  save('commands.json',commands);return (out/(name+'.stdout')).read_text()
 try:
  stage('checking_inputs')
  if shutil.disk_usage(ROOT).free<6*1024**3:raise RuntimeError('Insufficient disk reserve')
  if not a.reads.is_file():raise FileNotFoundError(a.reads)
  st=str(ROOT/'tools/readqc-env/bin/samtools');mm=str(ROOT/'tools/readqc-env/bin/minimap2')
  save('tool_versions.json',dict(minimap2=subprocess.check_output([mm,'--version'],text=True).strip(),samtools=subprocess.check_output([st,'--version'],text=True).splitlines()[0]))
  verification=a.reads.parent/(a.reads.name+'.verified.json')
  if not a.partial:
   if not verification.is_file():raise RuntimeError('Full-run QC requires a verified whole-file download manifest')
   manifest=json.loads(verification.read_text());observed=hashlib.md5()
   with a.reads.open('rb') as f:
    for block in iter(lambda:f.read(1024*1024),b''):observed.update(block)
   if not manifest.get('verified') or a.reads.stat().st_size!=manifest['expected_bytes'] or observed.hexdigest()!=manifest['expected_md5']:raise RuntimeError('Input no longer matches verified full FASTQ')
   (ROOT/'outputs/read_dataset_20261030/local_read_assets.json').write_text(json.dumps(dict(run='ERR4194695',path=str(a.reads.resolve()),whole_file_md5_verified=True,download_complete=True,assembly_use_verified=False),indent=2)+'\n')
  source=ROOT/'data/orthology/source'/ACC/'genome.fna.gz';genome=out/'genome.fna'
  sequences=read_fasta(source)
  if not genome.exists():write_fasta(genome,sequences)
  run([st,'faidx',str(genome)],'faidx')
  # FASTQ record integrity and quality representation are checked before filtering.
  stage('fastq_validation');n=bp=0;lengths=[];qual=collections.Counter();bad_names=0
  opener=gzip.open if str(a.reads).endswith('.gz') else open
  with opener(a.reads,'rt') as f:
   while True:
    h=f.readline()
    if not h:break
    s=f.readline().rstrip();plus=f.readline();q=f.readline().rstrip()
    if not h.startswith('@') or not plus.startswith('+') or not s or len(s)!=len(q):raise ValueError('Malformed FASTQ')
    n+=1;bp+=len(s);lengths.append(len(s));qual.update(q)
    if not re.match(r'[^/]+/\d+/\d+_\d+',h[1:].split()[0]):bad_names+=1
  quality=dict(reads=n,bases=bp,min_length=min(lengths),max_length=max(lengths),median_length=sorted(lengths)[len(lengths)//2],quality_character_counts=dict(qual),all_quality_values_zero=set(qual)=={'!'},original_pacbio_zmw_names_preserved=(bad_names==0))
  save('fastq_qc.json',quality)
  save('inputs.json',dict(assembly=ACC,run='ERR4194695',reads=str(a.reads.resolve()),reads_sha256=sha256(a.reads),genome_sha256=sha256(source),gff_sha256=sha256(ROOT/'data/orthology/source'/ACC/'annotation.gff.gz'),partial_input=a.partial,provenance='publication strain/project supported; exact run-to-version use unverified',biological_loss_inferred=False))
  bam=out/'aligned.bam'
  if not bam.exists():
   stage('mapping')
   args=[mm,'-ax','map-pb','--MD','-t','2','-K','50M',str(genome),str(a.reads.resolve())]
   sorter=[st,'sort','-@','1','-m','256M','-T',str(out/'sorttmp'),'-o',str(out/'aligned.bam.part'),'-']
   commands.extend([args,sorter]);save('commands.json',commands)
   with (out/'minimap2.stderr').open('w') as e1,(out/'sort.stderr').open('w') as e2:
    mapper=subprocess.Popen(args,stdout=subprocess.PIPE,stderr=e1)
    sorting=subprocess.Popen(sorter,stdin=mapper.stdout,stderr=e2);mapper.stdout.close()
    sortcode=sorting.wait();mapcode=mapper.wait()
   if sortcode or mapcode:raise RuntimeError(f'Mapping failed ({mapcode}, {sortcode}); no biological inference')
   (out/'aligned.bam.part').rename(bam)
  run([st,'quickcheck','-v',str(bam)],'quickcheck');run([st,'index',str(bam)],'index')
  run([st,'flagstat','-O','json',str(bam)],'flagstat')
  stage('regional_qc')
  run([st,'coverage','-q','20','-Q','0','--ff','3844',str(bam)],'whole_genome_coverage_mq20_bq0')
  lo=max(1,START-25000);hi=min(len(sequences[SEQ]),END+25000);region=f'{SEQ}:{lo}-{hi}'
  profiles={}
  for name,mq,bq in [('all_mappings',0,0),('mq20_bq0',20,0),('mq20_bq10',20,10)]:
   txt=run([st,'depth','-aa','-r',region,'-Q',str(mq),'-q',str(bq),'-G','3844',str(bam)],'depth_'+name)
   values={int(x.split('\t')[1]):int(x.split('\t')[2]) for x in txt.splitlines()}
   if len(values)!=hi-lo+1:raise ValueError('Depth did not return entire requested interval')
   profiles[name]=values
  features=read_gff(ROOT/'data/orthology/source'/ACC/'annotation.gff.gz')
  genes=[g for g in features if g['seqid']==SEQ and g['type']=='gene' and g['end']>=lo and g['start']<=hi]
  left=sorted([g for g in genes if g['end']<START],key=lambda g:g['end'],reverse=True)[:2]
  right=sorted([g for g in genes if g['start']>END],key=lambda g:g['start'])[:2]
  intervals=[dict(name='YJU3_homolog_candidate',start=START,end=END,role='candidate')]+[dict(name=g['attrs'].get('Name',g['attrs'].get('ID','gene')),start=max(lo,g['start']),end=min(hi,g['end']),role='neighbor_control') for g in left+right]
  stats=[]
  for interval in intervals:
   for name,depth in profiles.items():
    v=[depth[p] for p in range(interval['start'],interval['end']+1)]
    stats.append(dict(**interval,filter=name,bases=len(v),mean_depth=sum(v)/len(v),median_depth=sorted(v)[len(v)//2],min_depth=min(v),covered_fraction=sum(x>0 for x in v)/len(v),fraction_ge_5=sum(x>=5 for x in v)/len(v)))
  with (out/'region_depth_summary.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(stats[0]));w.writeheader();w.writerows(stats)
  reads=run([st,'view','-F','3844',str(bam),f'{SEQ}:{START}-{END}'],'candidate_primary_alignments')
  spans=[]
  for line in reads.splitlines():
   c=line.split('\t');flag=int(c[1]);pos=int(c[3]);mq=int(c[4]);cigar=re.findall(r'(\d+)([MIDNSHP=X])',c[5]);ref_len=sum(int(n) for n,op in cigar if op in 'MDN=X');end=pos+ref_len-1
   nm=next((int(x[5:]) for x in c[11:] if x.startswith('NM:i:')),None)
   spans.append(dict(read=c[0],start=pos,end=end,mapq=mq,reverse=bool(flag&16),nm=nm,cigar=c[5],spans_candidate_plus_100bp=(pos<=START-100 and end>=END+100),has_internal_reference_skip_or_deletion_ge100=any(int(n)>=100 and op in 'DN' for n,op in cigar)))
  save('candidate_read_support.json',spans)
  strong=[r for r in spans if r['mapq']>=20 and r['spans_candidate_plus_100bp'] and not r['has_internal_reference_skip_or_deletion_ge100']]
  candidate=next(x for x in stats if x['role']=='candidate' and x['filter']=='mq20_bq0')
  result=dict(assembly=ACC,run='ERR4194695',region=f'{SEQ}:{START}-{END}',partial_input=a.partial,fastq=quality,candidate_depth=candidate,primary_alignments_overlapping=len(spans),mq20_reads_spanning_plus100bp=len(strong),spanning_forward=sum(not x['reverse'] for x in strong),spanning_reverse=sum(x['reverse'] for x in strong),independent_zmw_count='N/A: archive FASTQ read names do not preserve original ZMW IDs' if bad_names else 'not_counted',assembly_use_verified=False,functional_status='not_assessed',biological_loss_inferred=False,interpretation='Read-supported candidate sequence presence; annotation-gap hypothesis remains; functional integrity and per-base consensus unverified' if strong else 'Insufficient positive spanning-read evidence; absence cannot be inferred',limitations=['Prefix preview is not random or full-depth data; no low-coverage or absence diagnosis is allowed.'] if a.partial else ['Mapping support is not a functional assay or proof of error-free assembly.','PacBio subreads may not be independent molecules.'])
  result['quality_encoding_assessment']=quality_encoding_assessment(dict(qual))
  result['regional_evidence_assessment']=regional_evidence_assessment(candidate['covered_fraction'],len(strong),a.partial,(not a.partial and (a.reads.parent/(a.reads.name+'.verified.json')).exists()))
  save('report.json',result)
  # Preserve a small auditable profile and compact scientific plot.
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
  fig,axes=plt.subplots(2,1,figsize=(11,5),sharex=True,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
  for name,color in [('mq20_bq0','#0072B2'),('mq20_bq10','#D55E00')]:
   x=list(range(lo,hi+1,100));y=[sum(profiles[name][p] for p in range(b,min(b+100,hi+1)))/len(range(b,min(b+100,hi+1))) for b in x];axes[0].plot(x,y,label=name,color=color)
  axes[0].axvspan(START,END,color='#E69F00',alpha=.3);axes[0].set_ylabel('Aligned base depth');axes[0].legend();axes[0].set_title('HN1 YJU3 homolog: '+('PREFIX PREVIEW (not random / incomplete)' if a.partial else 'full PacBio run'))
  for g in genes:
   axes[1].plot([g['start'],g['end']],[.5,.5],lw=6,color='#777777')
  axes[1].plot([START,END],[.8,.8],lw=6,color='#E69F00');axes[1].set_yticks([]);axes[1].set_xlabel(SEQ+' coordinate (1-based)');axes[1].set_ylim(0,1)
  fig.savefig(out/'regional_read_support.png',dpi=160);plt.close(fig)
  stage('complete');save('status.json',dict(stage='complete',finished=True,partial_input=a.partial,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  print(json.dumps(result,indent=2),flush=True)
 except Exception as e:
  save('status.json',dict(stage='failed',finished=False,error=str(e),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));raise
if __name__=='__main__':main()
