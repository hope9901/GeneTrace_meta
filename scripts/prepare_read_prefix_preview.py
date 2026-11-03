#!/usr/bin/env python3
"""Extract complete records from a contiguous download prefix, explicitly incomplete."""
import gzip,hashlib,json,pathlib,shutil
root=pathlib.Path(__file__).resolve().parents[1];parts=root/'data/read_dataset_20261030/fastq/ERR4194695_subreads.fastq.gz.chunks';out=root/'runs/read_qc_HN1_prefix_20261030';out.mkdir(parents=True,exist_ok=True)
chunks=[]
while (parts/f'{len(chunks):05d}').is_file():chunks.append(parts/f'{len(chunks):05d}')
if not chunks:raise RuntimeError('No contiguous completed chunks')
prefix=out/'compressed_prefix.gz.partial'
with prefix.open('wb') as target:
 for p in chunks:
  with p.open('rb') as src:shutil.copyfileobj(src,target)
count=0
with gzip.open(prefix,'rt') as src,(out/'reads.prefix.fastq').open('w') as dst:
 try:
  for _ in range(20000):
   lines=[src.readline() for _ in range(4)]
   if not all(lines) or not lines[3].endswith('\n'):break
   if len(lines[1].strip())!=len(lines[3].strip()):raise ValueError('Invalid record')
   dst.writelines(lines);count+=1
 except EOFError:pass
(out/'prefix_provenance.json').write_text(json.dumps(dict(run='ERR4194695',records=count,selection='first 20000 complete records or available prefix, whichever shorter; NOT RANDOM',full_run_md5_verified=False,compressed_prefix_bytes=prefix.stat().st_size,compressed_prefix_sha256=hashlib.sha256(prefix.read_bytes()).hexdigest(),chunks=[dict(name=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in chunks],absence_inference_allowed=False),indent=2)+'\n')
print(count,'complete preview records')
