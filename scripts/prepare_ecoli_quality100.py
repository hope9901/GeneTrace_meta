#!/usr/bin/env python3
"""Collect complete RefSeq inventory, select auditable quality panel and verify sources."""
import concurrent.futures,csv,datetime,gzip,hashlib,json,pathlib,shutil,subprocess,sys,time,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,sha256
OUT=ROOT/'outputs/ecoli_quality100_20261102';INV=ROOT/'data/bacteria_full/inventory';ANCHOR='GCF_000005845.2'
PARAMS={'page_size':1000,'filters.assembly_version':'current','filters.assembly_level':'complete_genome','filters.assembly_source':'refseq','filters.exclude_atypical':'true'}
def atomic(p,obj):
 tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def status(stage,**kw):atomic(OUT/'status.json',dict(stage=stage,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw))
def fetch(url,p):
 if p.exists():return p.read_bytes()
 if shutil.disk_usage(ROOT).free<6*1024**3:raise RuntimeError('Storage reserve reached (6 GiB)')
 p.parent.mkdir(parents=True,exist_ok=True);temp=p.with_suffix(p.suffix+'.part')
 r=subprocess.run(['curl','-fLsS','--retry','3','--retry-delay','2','--connect-timeout','20','--max-time','180',url,'-o',str(temp)],capture_output=True,text=True)
 if r.returncode:raise RuntimeError(url+': '+r.stderr[-300:])
 temp.replace(p);return p.read_bytes()
def quality(r):
 c=r.get('checkm_info',{});a=r.get('assembly_stats',{});total=int(a.get('total_sequence_length',0));ungap=int(a.get('total_ungapped_length',0))
 return (float(c.get('contamination',1e6)),-float(c.get('completeness',0)),(total-ungap)/total if total else 1,-int(a.get('contig_n50',0)),r['accession'])
def eligible(r):
 c=r.get('checkm_info',{});a=r.get('assembly_stats',{});ann=r.get('annotation_info',{});ani=r.get('average_nucleotide_identity',{})
 return r['accession'].startswith('GCF_') and r.get('assembly_info',{}).get('assembly_level')=='Complete Genome' and c.get('completeness',0)>=99 and c.get('contamination',100)<=1 and ann.get('stats',{}).get('gene_counts',{}).get('protein_coding',0)>0 and ani.get('taxonomy_check_status') in ('OK',None) and int(a.get('total_sequence_length',0))>0

def download(r):
 acc=r['accession'];info=r['assembly_info'];name=acc+'_'+info['assembly_name'];digits=acc.split('_')[1].split('.')[0];base='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/'+'/'.join(digits[i:i+3] for i in (0,3,6))+'/'+urllib.parse.quote(name,safe='')+'/'
 d=ROOT/'data/bacteria/source'/acc;d.mkdir(parents=True,exist_ok=True)
 checks={x.split()[1].removeprefix('./'):x.split()[0] for x in fetch(base+'md5checksums.txt',d/'md5checksums.txt').decode().splitlines() if len(x.split())==2};files=[]
 for suffix,local in [('genomic.fna.gz','genome.fna.gz'),('genomic.gff.gz','annotation.gff.gz'),('cds_from_genomic.fna.gz','cds.fna.gz'),('protein.faa.gz','protein.faa.gz'),('assembly_report.txt','assembly_report.txt')]:
  remote=name+'_'+suffix;p=d/local;data=fetch(base+urllib.parse.quote(remote,safe=''),p)
  if remote not in checks or hashlib.md5(data).hexdigest()!=checks[remote]:raise ValueError('Checksum mismatch: '+str(p))
  if local.endswith('.gz'):gzip.decompress(data)
  files.append(dict(file=local,url=base+urllib.parse.quote(remote,safe=''),md5=checks[remote],sha256=sha256(p)))
 genome=read_fasta(d/'genome.fna.gz');proteins=read_fasta(d/'protein.faa.gz');bases=''.join(genome.values());n=sum(bases.count(x) for x in 'ACGT');org=r['organism'];sample=info.get('biosample',{});strain=org.get('infraspecific_names',{}).get('strain') or sample.get('strain') or acc
 if acc==ANCHOR:strain='MG1655'
 row=dict(assembly_accession=acc,strain=strain,resolved_strain=strain,organism=org['organism_name'],assembly_level=info['assembly_level'],genome_size_bp=len(bases),gc_percent=100*(bases.count('G')+bases.count('C'))/n if n else None,number_of_proteins=len(proteins),sequence_count=len(genome),reference=acc==ANCHOR,reference_accession=ANCHOR,source_directory=str(d.relative_to(ROOT)),source_files=files,raw_reads_downloaded=False,raw_read_status='N/A: not surveyed',sequencing_group=info.get('sequencing_tech') or 'N/A',biosample=sample.get('accession','N/A'),bioproject=info.get('bioproject_accession','N/A'),submitter=info.get('submitter','N/A'),checkm_info=r.get('checkm_info',{}),assembly_method=info.get('assembly_method') or 'N/A',annotation_provider=r.get('annotation_info',{}).get('provider','N/A'),genome_fasta_sha256=files[0]['sha256'])
 atomic(d/'quality100_source_manifest.json',row);return row

def main():
 OUT.mkdir(parents=True,exist_ok=True);INV.mkdir(parents=True,exist_ok=True);reports=[];token=None;page=1
 try:
  while True:
   params=dict(PARAMS)
   if token:params['page_token']=token
   url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/562/dataset_report?'+urllib.parse.urlencode(params)
   raw=fetch(url,INV/f'complete_refseq_page_{page:05d}.json');s=json.loads(raw);reports.extend(s['reports']);token=s.get('next_page_token');status('collecting_candidates',records=len(reports),ncbi_total=s['total_count'],page=page);print('candidates',len(reports),'/',s['total_count'],flush=True)
   if not token:break
   page+=1;time.sleep(.4)
  atomic(OUT/'candidate_reports.json',reports)
  candidates=sorted([r for r in reports if eligible(r)],key=quality);selected=[];seen_samples=set();seen_strains=set()
  anchor=next(r for r in candidates if r['accession']==ANCHOR)
  for r in [anchor]+[r for r in candidates if r['accession']!=ANCHOR]:
   info=r['assembly_info'];sample=info.get('biosample',{}).get('accession');strain=r['organism'].get('infraspecific_names',{}).get('strain','').strip().casefold()
   if sample and sample in seen_samples:continue
   if strain and strain not in ['na','n/a','unknown','not provided','missing'] and strain in seen_strains:continue
   selected.append(r)
   if sample:seen_samples.add(sample)
   if strain:seen_strains.add(strain)
   if len(selected)==100:break
  if len(selected)!=100:raise ValueError('Insufficient eligible distinct strains')
  atomic(OUT/'selection.json',dict(candidate_count=len(reports),eligible_count=len(candidates),selected_count=len(selected),criteria=dict(scope=PARAMS,minimum_checkm_completeness=99,maximum_checkm_contamination=1,rank_order=['contamination ascending','completeness descending','gap fraction ascending','contig N50 descending','accession'],anchor_exception='MG1655 explicitly retained as reference; other 99 ranked',deduplication='BioSample or nonmissing casefolded strain name; unknown strains kept separately; not an ANI clone deduplication',limitations='Quality proxies, not ground truth; not population representative; no universal global best ranking'),assemblies=selected))
  with (OUT/'selection.csv').open('w',newline='') as f:
   w=csv.writer(f);w.writerow(['rank','accession','strain','biosample','completeness','contamination','contig_n50','genome_size_bp','anchor'])
   for rank,r in enumerate(selected,1):w.writerow([rank,r['accession'],r['organism'].get('infraspecific_names',{}).get('strain','N/A'),r['assembly_info'].get('biosample',{}).get('accession','N/A'),r['checkm_info']['completeness'],r['checkm_info']['contamination'],r['assembly_stats'].get('contig_n50','N/A'),r['assembly_stats']['total_sequence_length'],r['accession']==ANCHOR])
  complete=[];errors=[]
  with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
   futures={pool.submit(download,r):r['accession'] for r in selected}
   for future in concurrent.futures.as_completed(futures):
    try:complete.append(future.result());print('verified',len(complete),'/100',futures[future],flush=True)
    except Exception as e:errors.append(dict(accession=futures[future],error=str(e)));print('ERROR',futures[future],str(e),flush=True)
    atomic(OUT/'downloaded.json',complete);atomic(OUT/'errors.json',errors);status('downloading',completed=len(complete),total=100,failed=len(errors))
  if errors:raise RuntimeError('Some downloads failed; rerun to resume verified files')
  old=ROOT/'data/bacteria/cohort.json';archive=OUT/'previous_cohort.json'
  if not archive.exists():shutil.copy2(old,archive)
  pilot=ROOT/'outputs/cohort_comparison_20261102/b3987';backup=OUT/'pilot4_comparison'
  if pilot.exists() and not backup.exists():shutil.copytree(pilot,backup)
  atomic(old,sorted(complete,key=lambda r:r['assembly_accession']));status('preparing_gene_maps',completed=100,total=100)
  subprocess.run([sys.executable,str(ROOT/'scripts/prepare_bacterial_gene_maps.py')],check=True,cwd=ROOT)
  status('comparing_rpoB',completed=100,total=100)
  subprocess.run([sys.executable,str(ROOT/'scripts/compare_locus_cohort.py'),'--assembly',ANCHOR,'--gene','rpoB','--workers','2'],check=True,cwd=ROOT)
  subprocess.run([sys.executable,str(ROOT/'scripts/summarize_ecoli_quality100.py')],check=True,cwd=ROOT)
  status('complete',completed=100,total=100,notes='100 sequence/annotation sets verified; rpoB population comparison complete; raw reads and causal error diagnosis not performed')
 except Exception as e:status('failed_resumable',error=str(e));raise
if __name__=='__main__':main()
