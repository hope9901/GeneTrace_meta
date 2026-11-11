#!/usr/bin/env python3
"""Build an explicit allowlisted source/evidence archive; never include local .env/data/tools."""
import hashlib,json,pathlib,zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
 out=ROOT/'dist';out.mkdir(exist_ok=True);paths=[]
 for name in ['genetrace','web','demo','docs','scripts','tests','assets']:
  paths.extend(p for p in (ROOT/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc','.zip','.xlsx','.pdf'] and p.stat().st_size<3000000)
 paths.extend(ROOT/name for name in ['README.md','.env.example','.gitignore','requirements-analysis.txt','NOTICE.md'] if (ROOT/name).exists());paths=sorted(set(paths));files=[]
 for p in paths:
  rel=p.relative_to(ROOT)
  if p.name=='.env' or any(x in rel.parts for x in ['.venv','tools','runs','data','.git']):raise ValueError('Excluded path '+str(rel))
  files.append(dict(path=str(rel),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 archive=out/'GeneTrace-submission-source.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
  for p in paths:z.write(p,'GeneTrace/'+str(p.relative_to(ROOT)))
 (out/'manifest.json').write_text(json.dumps(dict(files=files,archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),source_only=True,includes_live_Meta_verification=(ROOT/'demo/live_meta_verification.json').exists(),license_owner_decision_pending=True),indent=2)+'\n');print(archive,'files',len(files),'bytes',archive.stat().st_size)
if __name__=='__main__':main()
