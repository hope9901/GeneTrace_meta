#!/usr/bin/env python3
"""Build a sample-linked read inventory; sample identity is not assembly-use proof."""
import argparse, concurrent.futures, csv, datetime, hashlib, io, json, pathlib, subprocess, urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1]
FIELDS='run_accession,sample_accession,secondary_sample_accession,study_accession,experiment_accession,instrument_platform,instrument_model,library_strategy,library_source,library_selection,library_layout,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes,sample_title'
def dump(p,x):
 p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def table(p,rows):
 if not rows:return
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);a=ap.parse_args()
 cohort=json.loads((ROOT/'outputs/orthology_20261030/selected_assemblies.json').read_text())
 out=ROOT/'outputs/read_dataset_20261030';out.mkdir(parents=True,exist_ok=True)
 cache=ROOT/'data/read_dataset_20261030/ena_reports';cache.mkdir(parents=True,exist_ok=True)
 def fetch(sample):
  p=cache/(sample+'.tsv');url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode(dict(accession=sample,result='read_run',fields=FIELDS,format='tsv'))
  if not p.exists():
   q=subprocess.run(['curl','-fLsS','--retry','2','--max-time','75',url],capture_output=True)
   if q.returncode:return sample,[],dict(status='query_failed',error=q.stderr.decode()[-500:],url=url)
   p.write_bytes(q.stdout)
  try:
   text=p.read_text();rows=list(csv.DictReader(io.StringIO(text),delimiter='\t'))
   if not text.startswith('run_accession\t'):raise ValueError('Unexpected API response')
   return sample,rows,dict(status='queried',url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
  except Exception as e:return sample,[],dict(status='query_failed',error=str(e),url=url)
 samples=sorted({r['biosample'] for r in cohort if r.get('biosample')})
 results={}
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  for i,(s,r,m) in enumerate(pool.map(fetch,samples),1):
   results[s]=(r,m)
   if i%30==0:print(f'Queried {i}/{len(samples)} samples',flush=True)
 runs=[]; assemblies=[];files=[]
 for row in cohort:
  acc=row['assembly_accession'];sample=row.get('biosample');rr,meta=results.get(sample,([],{'status':'no_biosample','url':'N/A'}))
  matched=[r for r in rr if sample in (r.get('sample_accession'),r.get('secondary_sample_accession'))]
  genomic=[r for r in matched if r['library_source']=='GENOMIC' and r['library_strategy']=='WGS']
  base=dict(row)
  for label,filename in [('genome','genome.fna.gz'),('gff','annotation.gff.gz'),('cds','cds.fna.gz'),('protein','protein.faa.gz')]:
   local=ROOT/'data/orthology/source'/acc/filename
   base[label+'_local_path']=str(local.relative_to(ROOT))
   base[label+'_local_exists']=local.exists()
  base.update(raw_query_status=meta['status'],raw_read_status=('sample_linked_WGS_available' if genomic else 'sample_linked_other_libraries_only' if matched else 'no_runs_returned' if meta['status']=='queried' else 'N/A'),raw_run_count=len(matched),wgs_run_count=len(genomic),raw_platforms=';'.join(sorted({r['instrument_platform'] for r in genomic})) or 'N/A',raw_run_accessions=';'.join(r['run_accession'] for r in genomic) or 'N/A',assembly_read_use_verified=False,raw_reads_downloaded=False,read_link_evidence='exact BioSample accession; assembly use unverified' if matched else 'N/A',raw_query_url=meta['url'])
  assemblies.append(base)
  for r in matched:
   runs.append(dict(assembly_accession=acc,assembly_biosample=sample,link_evidence='exact_BioSample_only',assembly_use_verified=False,ncbi_sra_url='https://www.ncbi.nlm.nih.gov/sra/'+r['run_accession'],**r))
   urls=r['fastq_ftp'].split(';') if r['fastq_ftp'] else []
   md5=r['fastq_md5'].split(';');sizes=r['fastq_bytes'].split(';')
   if urls and not(len(urls)==len(md5)==len(sizes)):raise ValueError('FASTQ metadata length mismatch')
   for u,h,n in zip(urls,md5,sizes):
    files.append(dict(assembly_accession=acc,biosample=sample,run_accession=r['run_accession'],library_strategy=r['library_strategy'],library_source=r['library_source'],platform=r['instrument_platform'],url='https://'+u,md5=h,bytes=int(n),downloaded=False,assembly_use_verified=False))
 dump(out/'assemblies.json',assemblies);table(out/'assemblies.csv',assemblies);table(out/'read_runs.csv',runs);table(out/'read_files.csv',files);dump(out/'query_provenance.json',{s:m for s,(_,m) in results.items()})
 summary=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),assemblies=len(assemblies),samples_queried=len(samples),query_failures=sum(m['status']=='query_failed' for _,m in results.values()),assemblies_with_sample_linked_wgs=sum(r['wgs_run_count']>0 for r in assemblies),assemblies_with_any_sample_linked_reads=sum(r['raw_run_count']>0 for r in assemblies),unique_runs=len({r['run_accession'] for r in runs}),unique_wgs_runs=len({r['run_accession'] for r in runs if r['library_source']=='GENOMIC' and r['library_strategy']=='WGS'}),fastq_bytes_unique=sum({r['url']:r['bytes'] for r in files}.values()),assembly_use_verified=0,fastq_downloaded=0,notes=['NCBI assembly cohort joined by exact BioSample to ENA INSDC read reports.','No run returned does not prove raw data absent. Failed queries are N/A.','Same BioSample does not prove these reads generated this assembly.','No biological loss or sequencing error can be diagnosed from this inventory alone.'])
 dump(out/'summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
