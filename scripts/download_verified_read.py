#!/usr/bin/env python3
"""Resumable, bounded HTTP range downloads with whole-file MD5 verification."""
import argparse,concurrent.futures,hashlib,json,pathlib,re,shutil,subprocess,time

def main():
 p=argparse.ArgumentParser();p.add_argument('--url',required=True);p.add_argument('--bytes',type=int,required=True);p.add_argument('--md5',required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
 a.output.parent.mkdir(parents=True,exist_ok=True);parts=a.output.parent/(a.output.name+'.chunks');parts.mkdir(exist_ok=True)
 chunk=4*1024*1024;total=(a.bytes+chunk-1)//chunk
 def md5(path):
  h=hashlib.md5()
  with path.open('rb') as f:
   for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
  return h.hexdigest()
 if a.output.exists():
  assert a.output.stat().st_size==a.bytes and md5(a.output)==a.md5
  print('Already verified',a.output);return
 def fetch(i):
  target=parts/f'{i:05d}';begin=i*chunk;end=min(a.bytes,begin+chunk)-1;expected=end-begin+1
  if target.exists() and target.stat().st_size==expected:return expected
  for attempt in range(8):
   # Fresh files prevent a failed new response from reusing an older response body.
   import uuid
   token=uuid.uuid4().hex;piece=parts/f'{i:05d}.{token}.piece';headers=parts/f'{i:05d}.{token}.headers'
   result=subprocess.run(['curl','-fLsS','--connect-timeout','15','--max-time','60','--range',f'{begin}-{end}','-D',str(headers),a.url,'-o',str(piece)],capture_output=True)
   header_text=headers.read_text(errors='replace') if headers.exists() else ''
   ranges=re.findall(r'content-range:\s*bytes\s+(\d+)-(\d+)/(\d+)',header_text,re.I)
   n=piece.stat().st_size if piece.exists() else 0
   if result.returncode==0 and re.search(r'HTTP/[^ ]+ 206 ',header_text) and ranges and tuple(map(int,ranges[-1]))==(begin,end,a.bytes) and n==expected:
    piece.replace(target);headers.replace(parts/f'{i:05d}.headers');return expected
   # Failed attempt files are never appended or reused.
   if piece.exists():piece.unlink()
   if headers.exists():headers.unlink()
   time.sleep(min(1+attempt,8))
  raise RuntimeError(f'Chunk {i} failed after retries: {result.stderr.decode()[-200:]}')
 done=0;last=time.monotonic()
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  futures=[pool.submit(fetch,i) for i in range(total)]
  for future in concurrent.futures.as_completed(futures):
   done+=future.result()
   if time.monotonic()-last>10:
    print(f'{done/a.bytes:.1%} complete ({done/1e6:.1f}/{a.bytes/1e6:.1f} MB)',flush=True);last=time.monotonic()
 combined=a.output.with_suffix(a.output.suffix+'.assembling')
 with combined.open('wb') as dst:
  for i in range(total):
   with (parts/f'{i:05d}').open('rb') as src:shutil.copyfileobj(src,dst)
 observed=md5(combined)
 if combined.stat().st_size!=a.bytes or observed!=a.md5:raise RuntimeError('Full-file size/MD5 mismatch')
 combined.rename(a.output)
 (a.output.parent/(a.output.name+'.verified.json')).write_text(json.dumps(dict(url=a.url,expected_bytes=a.bytes,expected_md5=a.md5,observed_md5=observed,verified=True),indent=2)+'\n')
 print('MD5 verified:',a.output,flush=True)
if __name__=='__main__':main()
