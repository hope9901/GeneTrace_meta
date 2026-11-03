#!/usr/bin/env python3
"""Continue full-run regional QC once the background download has been verified."""
import datetime,json,pathlib,subprocess,sys,time,shutil
root=pathlib.Path(__file__).resolve().parents[1];out=root/'outputs/read_qc_HN1_20261030';out.mkdir(parents=True,exist_ok=True)
reads=root/'data/read_dataset_20261030/fastq/ERR4194695_subreads.fastq.gz'
def status(stage,**extra):
 (out/'pipeline_status.json').write_text(json.dumps(dict(stage=stage,pid=__import__('os').getpid(),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**extra),indent=2)+'\n')
try:
 while not (reads.parent/(reads.name+'.verified.json')).exists():
  parts=reads.parent/(reads.name+'.chunks');complete=sum(p.stat().st_size for p in parts.glob('*') if p.name.isdigit()) if parts.exists() else 0
  status('waiting_for_verified_full_fastq',complete_chunk_bytes=complete,expected_bytes=927324001)
  if shutil.disk_usage(root).free<6*1024**3:raise RuntimeError('Disk reserve below 6 GiB')
  time.sleep(15)
 status('full_read_mapping')
 with (out/'pipeline.log').open('w') as log:
  subprocess.run([sys.executable,str(root/'scripts/run_region_read_qc.py'),'--reads',str(reads),'--output',str(root/'runs/read_qc_HN1_full_20261030')],stdout=log,stderr=subprocess.STDOUT,check=True)
 for name in ['report.json','region_depth_summary.csv','regional_read_support.png','fastq_qc.json','inputs.json','commands.json','tool_versions.json']:
  shutil.copy2(root/'runs/read_qc_HN1_full_20261030'/name,out/name)
 status('complete',finished=True)
except Exception as e:status('failed',error=str(e));raise
