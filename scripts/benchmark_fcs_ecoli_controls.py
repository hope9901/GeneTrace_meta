#!/usr/bin/env python3
"""Compare label-blind population alignability flags with external NCBI FCS spans.
Not an independent contamination classifier; does not tune thresholds to FCS labels.
"""
import collections,csv,gzip,hashlib,json,pathlib,shutil,subprocess,sys,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,write_fasta,sha256
D=ROOT/'data/contamination_validation_20261102';OUT=ROOT/'outputs/contamination_validation_20261102';RUN=ROOT/'runs/contamination_validation_20261102'
def dump(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def fetch(url,p):
 if p.exists():return p.read_bytes()
 t=p.with_suffix(p.suffix+'.part');subprocess.run(['curl','-fLsS','--retry','2','--connect-timeout','20','--max-time','180',url,'-o',str(t)],check=True);t.replace(p);return p.read_bytes()
def coverage(intervals):
 n=0;end=-1
 for a,b in sorted(intervals):
  if b>end:n+=b-max(a,end);end=b
 return n

def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True);controls=json.loads((D/'confirmed_ecoli_fcs_positive.json').read_text())['ecoli_controls'];controls=[c for c in controls if 'contaminated' in c['report']['assembly_info'].get('atypical',{}).get('warnings',[])][:3];dump(OUT/'selection.json',dict(controls=controls,selection_bias='First three confirmed E. coli controls from FCS high-contamination candidate screen; not random',labels_for_evaluation_only=True));cohort=json.loads((ROOT/'data/bacteria/cohort.json').read_text());unique={};mapping=[]
 try:
  dump(OUT/'status.json',dict(stage='preparing_population_reference'))
  for r in cohort:
   for seqid,seq in read_fasta(ROOT/r['source_directory']/'genome.fna.gz').items():
    digest=hashlib.sha256(seq.encode()).hexdigest();unique.setdefault(digest,seq);mapping.append(dict(assembly=r['assembly_accession'],seqid=seqid,reference_id=digest))
  ref=RUN/'population_unique_contigs.fna';write_fasta(ref,unique);dump(RUN/'population_reference_manifest.json',dict(assemblies=100,unique_sequences=len(unique),sha256=sha256(ref),mapping=mapping));evaluations=[]
  for index,c in enumerate(controls,1):
   r=c['report'];acc=r['accession'];name=acc+'_'+r['assembly_info']['assembly_name'];digits=acc.split('_')[1].split('.')[0];base='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/'+'/'.join(digits[i:i+3] for i in (0,3,6))+'/'+urllib.parse.quote(name,safe='')+'/';d=D/acc;d.mkdir(exist_ok=True);dump(OUT/'status.json',dict(stage='downloading_control',assembly=acc,completed=index-1,total=len(controls)))
   checks={x.split()[1].removeprefix('./'):x.split()[0] for x in fetch(base+'md5checksums.txt',d/'md5checksums.txt').decode().splitlines() if len(x.split())==2};files=[]
   for suffix,local in [('genomic.fna.gz','genome.fna.gz'),('fcs_report.txt','fcs_report.txt')]:
    remote=name+'_'+suffix;p=d/local;raw=fetch(base+urllib.parse.quote(remote,safe=''),p)
    if remote in checks and hashlib.md5(raw).hexdigest()!=checks[remote]:raise ValueError('Control file checksum mismatched '+remote)
    if local=='genome.fna.gz' and remote not in checks:raise ValueError('Genome lacks publisher checksum')
    files.append(dict(path=str(p.relative_to(ROOT)),md5=checks.get(remote),publisher_md5_available=remote in checks,sha256=sha256(p),note='Publisher MD5 absent for some auxiliary FCS reports; SHA256 records the retrieved snapshot.'))
   sequences=read_fasta(d/'genome.fna.gz');query=RUN/(acc+'.fna');write_fasta(query,sequences);paf=RUN/(acc+'.paf');cmd=[str(ROOT/'tools/readqc-env/bin/minimap2'),'-x','asm20','--secondary=no','-t','2',str(ref),str(query)];dump(OUT/'status.json',dict(stage='population_alignment',assembly=acc,completed=index-1,total=len(controls)))
   if shutil.disk_usage(ROOT).free<6*1024**3:raise RuntimeError('6 GiB disk reserve reached')
   with paf.open('w') as out,(RUN/(acc+'.stderr')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
   intervals=collections.defaultdict(list)
   for line in paf.read_text().splitlines():
    f=line.split('\t')
    if int(f[10])>=1000 and int(f[9])/int(f[10])>=.9:intervals[f[0]].append((int(f[2]),int(f[3])))
   # Label-blind output is finalized before the external FCS report is read.
   flags=[dict(seqid=seqid,length=len(seq),aligned_fraction=coverage(intervals[seqid])/len(seq),flag_low_population_alignability=len(seq)>=2000 and coverage(intervals[seqid])/len(seq)<.5) for seqid,seq in sequences.items()];dump(OUT/(acc+'_blind_flags.json'),flags)
   fcs=[]
   for line in (d/'fcs_report.txt').read_text().splitlines():
    if line.startswith('#') or not line:continue
    f=line.split('\t')
    if len(f)>=8:fcs.append(dict(seqid=f[0],start=int(f[1]),end=int(f[2]),action=f[4],contam_type=f[5],details=f[7]))
   positive={x['seqid'] for x in fcs if x['action'] in ['EXCLUDE','TRIM','FIX']};evaluated=[dict(**f,external_FCS_action_on_contig=f['seqid'] in positive) for f in flags];positives=[f for f in evaluated if f['external_FCS_action_on_contig']];nonlabeled=[f for f in evaluated if not f['external_FCS_action_on_contig']]
   result=dict(assembly=acc,ncbi_contaminated_warning=True,fcs_summary=c['fcs'],files=files,command=cmd,contigs=evaluated,fcs_spans=fcs,FCS_action_contigs=len(positives),FCS_action_contigs_flagged=sum(f['flag_low_population_alignability'] for f in positives),non_FCS_action_contigs=len(nonlabeled),non_FCS_action_contigs_flagged=sum(f['flag_low_population_alignability'] for f in nonlabeled),interpretation='Contig-level concordance, not sensitivity/specificity: FCS labels are incomplete and overlap is not span-resolution validation; unmatched genuine accessory/mobile DNA can be flagged.',independent_contamination_diagnosis=False)
   dump(OUT/(acc+'_evaluation.json'),result);evaluations.append({k:v for k,v in result.items() if k not in ['contigs','fcs_spans']});dump(OUT/'summary.json',dict(evaluations=evaluations,thresholds=dict(min_contig_bp=2000,min_alignment_bp=1000,min_alignment_identity=.9,flag_aligned_fraction_below=.5),thresholds_tuned_to_labels=False,FCS_used_for_inference=False,causal_diagnosis=False));print(acc,'FCS contigs',len(positives),'flagged',result['FCS_action_contigs_flagged'],flush=True)
  dump(OUT/'status.json',dict(stage='complete',completed=len(evaluations),total=len(controls),independent_contamination_diagnosis=False))
 except Exception as e:dump(OUT/'status.json',dict(stage='failed_resumable',error=str(e)));raise
if __name__=='__main__':main()
