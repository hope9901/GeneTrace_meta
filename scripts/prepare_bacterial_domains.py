#!/usr/bin/env python3
"""Retrieve checksum-verified model E. coli strains and selected Pfam profiles."""
import concurrent.futures,csv,gzip,hashlib,json,pathlib,re,subprocess,sys,datetime
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,sha256
MODELS=['PF00005','PF00072','PF00512','PF02518','PF00069','PF00076']
STRAINS={'MG1655':'GCF_000005845.2','W3110':'GCF_000010245.2','REL606':'GCF_000017985.1','BL21_DE3':'GCF_000022665.1'}
def fetch(url,path):
 path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():return path.read_bytes()
 temp=path.with_suffix(path.suffix+'.part');p=subprocess.run(['curl','-fLsS','--retry','2','--connect-timeout','15','--max-time','120',url,'-o',str(temp)],capture_output=True)
 if p.returncode:raise RuntimeError(url+': '+p.stderr.decode()[-200:])
 temp.replace(path);return path.read_bytes()
def assembly(item):
 strain,acc=item;directory=ROOT/'data/bacteria/source'/acc
 digits=acc.split('_')[1].split('.')[0];parent='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/'+ '/'.join(digits[i:i+3] for i in (0,3,6))+'/'
 index=fetch(parent,directory/'ftp_directory.html').decode();folders=re.findall(r'href="('+re.escape(acc)+r'_[^"/]+/)',index)
 if not folders:raise ValueError('Requested version not found at NCBI: '+acc)
 folder=sorted(set(folders))[0].rstrip('/');base=parent+folder+'/'
 md5=fetch(base+'md5checksums.txt',directory/'md5checksums.txt').decode();checks={line.split()[1].removeprefix('./'):line.split()[0] for line in md5.splitlines() if len(line.split())==2}
 files=[]
 for suffix,local in [('genomic.fna.gz','genome.fna.gz'),('genomic.gff.gz','annotation.gff.gz'),('cds_from_genomic.fna.gz','cds.fna.gz'),('protein.faa.gz','protein.faa.gz'),('assembly_report.txt','assembly_report.txt')]:
  name=folder+'_'+suffix;url=base+name;p=directory/local;data=fetch(url,p)
  if hashlib.md5(data).hexdigest()!=checks[name]:raise ValueError('MD5 mismatch '+name)
  if local.endswith('.gz'):gzip.decompress(data)
  files.append(dict(file=local,url=url,md5=checks[name],sha256=sha256(p)))
 text=(directory/'assembly_report.txt').read_text();metadata={}
 for line in text.splitlines():
  if line.startswith('# ') and ':' in line:
   k,v=line[2:].split(':',1);metadata[k.strip()]=v.strip()
 if not metadata.get('Organism name','').startswith('Escherichia coli'):raise ValueError('Wrong species')
 genome=read_fasta(directory/'genome.fna.gz');proteins=read_fasta(directory/'protein.faa.gz');bases=''.join(genome.values());acgt=sum(bases.count(x) for x in 'ACGT')
 row=dict(assembly_accession=acc,strain=strain,organism=metadata['Organism name'],assembly_level=metadata.get('Assembly level','N/A'),genome_size_bp=len(bases),gc_percent=100*(bases.count('G')+bases.count('C'))/acgt,number_of_proteins=len(proteins),sequence_count=len(genome),reference=(strain=='MG1655'),source_directory=str(directory.relative_to(ROOT)),source_files=files,raw_reads_downloaded=False)
 (directory/'source_manifest.json').write_text(json.dumps(row,indent=2)+'\n');print(strain,acc,'verified',flush=True);return row
def model(acc):
 directory=ROOT/'data/domain_models';url=f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=hmm';raw=fetch(url,directory/(acc+'.hmm.download'));data=gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
 text=data.decode()
 if not text.startswith('HMMER3/') or not re.search(r'^ACC\s+'+acc+r'\.',text,re.M) or not re.search(r'^GA\s+',text,re.M):raise ValueError('Invalid model or missing Pfam GA cutoff: '+acc)
 (directory/(acc+'.hmm')).write_bytes(data)
 row=dict(accession=re.search(r'^ACC\s+(\S+)',text,re.M)[1],name=re.search(r'^NAME\s+(\S+)',text,re.M)[1],length=int(re.search(r'^LENG\s+(\d+)',text,re.M)[1]),ga=re.search(r'^GA\s+(.+)',text,re.M)[1],url=url,sha256=sha256(directory/(acc+'.hmm')))
 print(row['accession'],row['name'],'verified',flush=True);return row
def main():
 out=ROOT/'outputs/domain_bacteria_20261101';out.mkdir(parents=True,exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(assembly,STRAINS.items()))
 (ROOT/'data/bacteria/cohort.json').write_text(json.dumps(rows,indent=2)+'\n')
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:models=list(pool.map(model,MODELS))
 (ROOT/'data/domain_models/manifest.json').write_text(json.dumps(models,indent=2)+'\n')
 (ROOT/'data/domain_models/panel.hmm').write_text(''.join((ROOT/'data/domain_models'/(m+'.hmm')).read_text() for m in MODELS))
 (out/'source_summary.json').write_text(json.dumps(dict(retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),bacterial_assemblies=rows,profiles=models),indent=2)+'\n')
 print('All source files and models verified',flush=True)
if __name__=='__main__':main()
