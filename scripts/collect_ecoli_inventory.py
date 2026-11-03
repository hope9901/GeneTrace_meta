#!/usr/bin/env python3
"""Resumable full NCBI taxon 562 inventory; retain records, collapse paired assemblies.

No sequence download until a storage destination has been selected. Accession pairs
are not strain deduplication; BioSample and strain information remain available.
"""
import argparse,csv,datetime,gzip,json,pathlib,shutil,sqlite3,subprocess,time,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/562/dataset_report'
def atomic(path,obj):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--directory',type=pathlib.Path,default=ROOT/'data/bacteria_full/inventory');a=ap.parse_args();d=a.directory;d.mkdir(parents=True,exist_ok=True)
 db=sqlite3.connect(d/'inventory.sqlite');db.execute('CREATE TABLE IF NOT EXISTS records(accession TEXT PRIMARY KEY, canonical TEXT, report TEXT)');db.execute('CREATE INDEX IF NOT EXISTS canonical_idx ON records(canonical)');db.commit()
 checkpoint=json.loads((d/'checkpoint.json').read_text()) if (d/'checkpoint.json').exists() else {'next_page':1,'next_page_token':None,'total_count':None}
 page=checkpoint['next_page'];token=checkpoint['next_page_token'];started=time.monotonic()
 def status(stage,**kw):
  n=db.execute('SELECT COUNT(*) FROM records').fetchone()[0];pairs=db.execute('SELECT COUNT(DISTINCT canonical) FROM records').fetchone()[0]
  atomic(d/'status.json',dict(stage=stage,taxon=562,records_collected=n,paired_assembly_groups=pairs,ncbi_total_records=checkpoint.get('total_count'),next_page=page,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),free_bytes=shutil.disk_usage(d).free,**kw))
 if checkpoint.get('complete'):status('complete');return
 try:
  while True:
   if shutil.disk_usage(d).free<6*1024**3:status('storage_required',reason='6 GiB reserve reached; inventory checkpoint retained');return
   archive=d/f'page_{page:05d}.json.gz';plain=d/f'page_{page:05d}.json';params={'page_size':1000,'filters.assembly_version':'current','filters.exclude_atypical':'false'}
   if token:params['page_token']=token
   url=BASE+'?'+urllib.parse.urlencode(params)
   if archive.exists():raw=gzip.decompress(archive.read_bytes())
   elif plain.exists():raw=plain.read_bytes()
   else:
    temp=d/'response.part';p=subprocess.run(['curl','-fLsS','--retry','4','--retry-delay','3','--connect-timeout','20','--max-time','180',url,'-o',str(temp)],capture_output=True,text=True)
    if p.returncode:raise RuntimeError(p.stderr[-600:])
    raw=temp.read_bytes();temp.unlink()
   data=json.loads(raw)
   if 'reports' not in data:raise ValueError('NCBI response has no reports')
   if not archive.exists():
    part=archive.with_suffix('.part');part.write_bytes(gzip.compress(raw));part.replace(archive)
   for r in data['reports']:
    acc=r['accession'];paired=r.get('paired_accession') or r.get('assembly_info',{}).get('paired_assembly',{}).get('accession');canonical=min([acc,paired]) if paired else acc
    db.execute('INSERT OR REPLACE INTO records VALUES(?,?,?)',(acc,canonical,json.dumps(r,separators=(',',':'))))
   db.commit();token=data.get('next_page_token');page+=1
   checkpoint=dict(next_page=page,next_page_token=token,total_count=data.get('total_count'),complete=not bool(token),source_endpoint=BASE,scope='current taxon 562 descendants; includes atypical; all assembly levels; paired records retained',started_utc=checkpoint.get('started_utc',datetime.datetime.now(datetime.timezone.utc).isoformat()))
   atomic(d/'checkpoint.json',checkpoint);status('collecting' if token else 'exporting');print('page',page-1,'records',db.execute('SELECT COUNT(*) FROM records').fetchone()[0],'/',checkpoint['total_count'],flush=True)
   if not token:break
   time.sleep(.4)
  columns=['assembly_accession','paired_accession','strain','biosample','assembly_level','genome_size_bp','gc_percent','protein_coding_genes','annotation_provider','sequencing_technology','assembly_method','bioproject','organism_name']
  with (d/'assemblies.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
   # Prefer RefSeq for each explicitly paired group, retaining original records in SQLite.
   query='SELECT r.report FROM records r WHERE r.accession=(SELECT MAX(s.accession) FROM records s WHERE s.canonical=r.canonical) ORDER BY r.canonical'
   for (raw,) in db.execute(query):
    r=json.loads(raw);info=r.get('assembly_info',{});stats=r.get('assembly_stats',{});ann=r.get('annotation_info',{});org=r.get('organism',{});sample=info.get('biosample',{})
    w.writerow(dict(assembly_accession=r['accession'],paired_accession=r.get('paired_accession','N/A'),strain=org.get('infraspecific_names',{}).get('strain','N/A'),biosample=sample.get('accession','N/A'),assembly_level=info.get('assembly_level','N/A'),genome_size_bp=stats.get('total_sequence_length','N/A'),gc_percent=stats.get('gc_percent','N/A'),protein_coding_genes=ann.get('stats',{}).get('gene_counts',{}).get('protein_coding','N/A'),annotation_provider=ann.get('provider','N/A'),sequencing_technology=info.get('sequencing_tech','N/A'),assembly_method=info.get('assembly_method','N/A'),bioproject=info.get('bioproject_accession','N/A'),organism_name=org.get('organism_name','N/A')))
  status('complete',elapsed_seconds=time.monotonic()-started)
 except Exception as e:status('failed_resumable',error=str(e));raise
 finally:db.close()
if __name__=='__main__':main()
