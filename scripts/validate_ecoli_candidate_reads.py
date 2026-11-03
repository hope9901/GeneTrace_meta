#!/usr/bin/env python3
"""Whole-run checksum-verified paired-read validation for two candidate genomes."""
import collections,datetime,gzip,hashlib,json,pathlib,shutil,subprocess,sys,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff,write_fasta,sha256
from genetrace.regional_qc import analyze_sam
from genetrace.read_provenance import assess_read_provenance
OUT=ROOT/'outputs/ecoli100_raw_validation_20261102';RUN=ROOT/'runs/ecoli100_raw_validation_20261102';DATA=ROOT/'data/ecoli100_reads_20261102/fastq'
def dump(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def main():
 for d in [OUT,RUN,DATA]:d.mkdir(parents=True,exist_ok=True)
 audits=json.loads((ROOT/'outputs/ecoli100_candidate_audit_20261102/results.json').read_text());metadata={r['assembly_accession']:r for r in json.loads((ROOT/'data/bacteria/cohort.json').read_text())};targets=['GCF_013112335.1','GCF_001558995.2'];summaries=[]
 try:
  for acc in targets:
   candidates=[r for r in audits if r['assembly_accession']==acc and r.get('both_flanks') and r['best_genome_hit']];run=candidates[0]['run_candidate'];rid=run['run_accession'];d=RUN/acc;d.mkdir(exist_ok=True);dump(OUT/'status.json',dict(stage='downloading',assembly=acc,run=rid));reads=[];verification=[];commands=[]
   for ftp,md5,size in zip(run['fastq_ftp'].split(';'),run['fastq_md5'].split(';'),run['fastq_bytes'].split(';')):
    if shutil.disk_usage(ROOT).free<6*1024**3:raise RuntimeError('6 GiB disk reserve reached')
    url=ftp if ftp.startswith('https://') else 'https://'+ftp;p=DATA/pathlib.Path(urllib.parse.urlsplit(url).path).name
    if not p.exists():
     temp=p.with_suffix(p.suffix+'.part');cmd=['curl','-fLsS','--retry','0','--connect-timeout','20','--max-time','900',url,'-o',str(temp)];cmd=cmd[:1]+(['--continue-at','-'] if temp.exists() and temp.stat().st_size else [])+cmd[1:];commands.append(cmd);result=subprocess.run(cmd)
     if result.returncode==6:
      lookup=json.loads((ROOT/'data/ecoli100_reads_20261102/ftp_dns_lookup.json').read_text());ip=next(x['data'] for x in lookup.get('Answer',[]) if x['type']==1);cmd=cmd[:1]+['--resolve','ftp.sra.ebi.ac.uk:443:'+ip]+cmd[1:];commands.append(cmd);subprocess.run(cmd,check=True)
     elif result.returncode:raise RuntimeError('FASTQ download failed')
     temp.replace(p)
    digest=hashlib.md5(p.read_bytes()).hexdigest()
    if digest!=md5 or p.stat().st_size!=int(size):raise ValueError('FASTQ checksum/size mismatch '+str(p))
    n=bases=0;qual=collections.Counter()
    with gzip.open(p,'rt') as f:
     while True:
      h=f.readline()
      if not h:break
      seq=f.readline().strip();plus=f.readline();q=f.readline().strip()
      if not h.startswith('@') or not plus.startswith('+') or len(seq)!=len(q) or not seq:raise ValueError('Malformed FASTQ')
      n+=1;bases+=len(seq);qual.update(q)
    reads.append(str(p));verification.append(dict(path=str(p.relative_to(ROOT)),md5=digest,expected_md5=md5,bytes=int(size),records=n,bases=bases,quality_character_counts=dict(qual),verified=True))
   if len(verification)==2 and verification[0]['records']!=verification[1]['records']:raise ValueError('Paired FASTQ record count mismatch')
   dump(d/'fastq_verification.json',verification);genome=read_fasta(ROOT/metadata[acc]['source_directory']/'genome.fna.gz');fasta=d/'genome.fna';write_fasta(fasta,genome);mm=ROOT/'tools/readqc-env/bin/minimap2';st=ROOT/'tools/readqc-env/bin/samtools';bam=d/'aligned.bam';signature=dict(reference_sha256=sha256(fasta),reads=[v['md5'] for v in verification],preset='sr');sidecar=d/'alignment_inputs.json'
   dump(OUT/'status.json',dict(stage='mapping',assembly=acc,run=rid))
   if not bam.exists() or not sidecar.exists() or json.loads(sidecar.read_text())!=signature:
    sam=d/'aligned.sam';cmd=[str(mm),'-ax','sr','-t','2',str(fasta)]+reads;commands.append(cmd)
    with sam.open('w') as out,(d/'minimap2.stderr').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
    cmd=[str(st),'sort','-@','1','-o',str(bam),str(sam)];commands.append(cmd);subprocess.run(cmd,check=True);dump(sidecar,signature);sam.unlink()
   subprocess.run([str(st),'index',str(bam)],check=True);features=read_gff(ROOT/metadata[acc]['source_directory']/'annotation.gff.gz');regions=[]
   for c in candidates:
    h=c['best_genome_hit'];overlap=[f for f in features if f['type'] in ['gene','pseudogene'] and f['seqid']==h['sseqid'] and f['start']<=h['end'] and f['end']>=h['start']];lo=min([h['start']]+[f['start'] for f in overlap]);hi=max([h['end']]+[f['end'] for f in overlap]);regions.append(dict(name=c['reference_gene'],seqid=h['sseqid'],start=lo,end=hi,role='candidate_annotation_bounded_region' if overlap else 'candidate_HSP_span_not_full_gene',protein_query_coverage=h['query_coverage']))
   # Add two annotated neighbors per genome as local context, not a universal normal baseline.
   first=regions[0];neighbors=sorted([f for f in features if f['type']=='CDS' and f['seqid']==first['seqid'] and (f['end']<first['start'] or f['start']>first['end'])],key=lambda f:min(abs(f['start']-first['end']),abs(f['end']-first['start'])))[:2]
   regions += [dict(name=f['attrs'].get('gene',f['attrs'].get('locus_tag','neighbor')),seqid=f['seqid'],start=f['start'],end=f['end'],role='local_neighbor_control') for f in neighbors]
   first=regions[0];sample_cmd=[str(st),'view',str(bam),f"{first['seqid']}:{max(1,first['start']-10000)}-{first['end']+10000}"];commands.append(sample_cmd);sample_sam=subprocess.check_output(sample_cmd,text=True);qualities=collections.Counter()
   for v in verification:qualities.update(v['quality_character_counts'])
   provenance=assess_read_provenance(sample_sam,qualities);dump(d/'provenance_review.json',provenance)
   reports=[]
   for i,r in enumerate(regions):
    cmd=[str(st),'view',str(bam),f"{r['seqid']}:{max(1,r['start']-1000)}-{r['end']+1000}"];commands.append(cmd);samtext=subprocess.check_output(cmd,text=True);report=analyze_sam(samtext,genome[r['seqid']],r['seqid'],r['start'],r['end']);report.update(read_provenance=provenance,**r,assembly=acc,run=rid,exact_biosample_match=True,assembly_run_version_use_verified=False,causal_diagnosis='unresolved',pair_independence='Paired records are not independent fragments; metrics count alignment records, not molecules.');dump(OUT/(acc+'_'+str(i)+'.json'),report);reports.append(dict(**r,signals=report['signals'],metrics=report['metrics'],coherent_substitution_count=len(report['coherent_substitutions']),coherent_indel_count=sum(x['coherent'] for x in report['indels']),quality_information_unresolved=report['quality_information_unresolved']))
   dump(d/'manifest.json',dict(assembly=acc,run_metadata=run,fastq=verification,alignment=signature,commands=commands,version=subprocess.check_output([str(mm),'--version'],text=True).strip(),causal_diagnosis='unresolved'));summaries.append(dict(assembly=acc,run=rid,regions=reports,whole_run_files_verified=True,read_provenance=provenance,assembly_run_version_use_verified=False));dump(OUT/'summary.json',dict(assemblies=summaries,causal_error_diagnosis=False,biological_loss_inferred=False));print(acc,rid,'validated',flush=True)
  dump(OUT/'status.json',dict(stage='complete',assemblies=len(summaries),causal_diagnosis='unresolved'))
 except Exception as e:dump(OUT/'status.json',dict(stage='failed_resumable',error=str(e)));raise
if __name__=='__main__':main()
