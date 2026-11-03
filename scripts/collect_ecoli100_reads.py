#!/usr/bin/env python3
"""Exact-BioSample WGS run inventory; assembly-version usage remains unresolved."""
import concurrent.futures,csv,datetime,hashlib,io,json,pathlib,subprocess,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/ecoli100_reads_20261102';CACHE=ROOT/'data/ecoli100_reads_20261102/ena_reports'
FIELDS='run_accession,sample_accession,secondary_sample_accession,study_accession,instrument_platform,instrument_model,library_strategy,library_source,library_layout,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes'
def dump(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True);cohort=json.loads((ROOT/'data/bacteria/cohort.json').read_text());runs=[];results=[]
 def fetch(row):
  sample=row.get('biosample','N/A');acc=row['assembly_accession']
  if not sample.startswith(('SAMN','SAMEA','SAMD')):return dict(assembly_accession=acc,biosample=sample,status='no_biosample',runs=[])
  p=CACHE/(sample+'.tsv');url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode(dict(accession=sample,result='read_run',fields=FIELDS,format='tsv'))
  if not p.exists():
   q=subprocess.run(['curl','-fLsS','--retry','2','--connect-timeout','15','--max-time','75',url],capture_output=True)
   if q.returncode:return dict(assembly_accession=acc,biosample=sample,status='query_failed',error=q.stderr.decode()[-300:],runs=[])
   if not q.stdout.startswith(b'run_accession\t'):return dict(assembly_accession=acc,biosample=sample,status='invalid_response',runs=[])
   p.write_bytes(q.stdout)
  rr=list(csv.DictReader(io.StringIO(p.read_text()),delimiter='\t'));matched=[r for r in rr if sample in [r.get('sample_accession'),r.get('secondary_sample_accession')] and r['library_strategy']=='WGS' and r['library_source']=='GENOMIC']
  return dict(assembly_accession=acc,biosample=sample,status='sample_linked_WGS_available' if matched else 'no_exact_sample_WGS_returned',runs=matched,assembly_read_use_verified=False,query_url=url,response_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  for result in pool.map(fetch,cohort):
   results.append(result)
   for run in result['runs']:runs.append(dict(assembly_accession=result['assembly_accession'],assembly_biosample=result['biosample'],assembly_read_use_verified=False,**run))
   dump(OUT/'assemblies.json',results);dump(OUT/'status.json',dict(stage='querying',completed=len(results),total=100,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
   print(len(results),'/100',result['status'],flush=True)
 if runs:
  with (OUT/'runs.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(runs[0]));w.writeheader();w.writerows(runs)
 dump(OUT/'summary.json',dict(assemblies=len(results),assemblies_with_WGS=sum(bool(r['runs']) for r in results),runs=len(runs),query_failures=sum(r['status'] in ['query_failed','invalid_response'] for r in results),raw_reads_downloaded=False,assembly_read_use_verified=False))
 dump(OUT/'status.json',dict(stage='complete',completed=len(results),total=100,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
if __name__=='__main__':main()
