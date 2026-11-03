"""Build a small public-data execution case from MD5-verified original FASTQ records."""
import gzip,hashlib,json,pathlib,re,shutil,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import sha256
OUT=ROOT/'demo/execution/nimt';OUT.mkdir(parents=True,exist_ok=True)
acc='GCF_001558995.2';source=ROOT/'data/bacteria/source'/acc
summary=json.loads((ROOT/'outputs/ecoli100_raw_validation_20261102/summary.json').read_text());a=next(x for x in summary['assemblies'] if x['assembly']==acc)
regions=[{k:r[k] for k in ['name','seqid','start','end','role']} for r in a['regions']]
bam=ROOT/'runs/ecoli100_raw_validation_20261102'/acc/'aligned.bam';st=ROOT/'tools/readqc-env/bin/samtools'
selection=f"{regions[0]['seqid']}:{min(r['start'] for r in regions)-5000}-{max(r['end'] for r in regions)+5000}"
sam=subprocess.check_output([str(st),'view',str(bam),selection],text=True)
def name(s):return re.sub(r'/[12]$','',s.split()[0].lstrip('@'))
names={name(l.split('\t')[0]) for l in sam.splitlines()};counts=[];source_files=[]
verification=json.loads((bam.parent/'fastq_verification.json').read_text())
for mate in [1,2]:
 p=ROOT/f'data/ecoli100_reads_20261102/fastq/SRR2762258_{mate}.fastq.gz';md5=hashlib.md5(p.read_bytes()).hexdigest();v=verification[mate-1];assert md5==v['expected_md5']
 out=OUT/f'reads_{mate}.fastq.gz';n=0
 with gzip.open(p,'rt') as inp,gzip.open(out,'wt') as dest:
  while True:
   row=[inp.readline() for _ in range(4)]
   if not row[0]:break
   assert all(row) and row[0].startswith('@') and len(row[1].strip())==len(row[3].strip())
   if name(row[0]) in names:dest.writelines(row);n+=1
 counts.append(n);source_files.append(dict(run='SRR2762258',mate=mate,whole_file_md5=md5,whole_file_bytes=p.stat().st_size,selected_records=n))
assert counts[0]==counts[1] and counts[0]>0
for name_ in ['genome.fna.gz','annotation.gff.gz']:shutil.copyfile(source/name_,OUT/name_)
manifest=dict(assembly=acc,strain='JEONG-1266',run='SRR2762258',sample='SAMN04191546',read_platform='Illumina MiSeq',assembly_platform='PacBio',regions=regions,selection=dict(method='Original paired FASTQ records selected by names overlapping a precomputed BAM window; sequence and quality retained unchanged',window=selection,not_an_independent_sample=True,partial_input=True,limitations='Reference-dependent selection can miss unmapped/divergent reads; not whole-run QC, unbiased depth, absence evidence, or independent validation.'),source_fastq=source_files,exact_assembly_run_use_verified=False,files={p.name:dict(sha256=sha256(p),bytes=p.stat().st_size) for p in OUT.iterdir() if p.name!='manifest.json'})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(read_pairs=counts[0],total_bytes=sum(p.stat().st_size for p in OUT.iterdir())),indent=2))
