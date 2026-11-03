#!/usr/bin/env python3
"""Audit HN1 candidate and neighboring controls using the existing full-run BAM."""
import csv,datetime,json,pathlib,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,sha256,write_json
from genetrace.regional_qc import analyze_sam
OUT=ROOT/'outputs/regional_validation_20261102';RUN=ROOT/'runs/regional_validation_20261102'

def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True)
 source=ROOT/'runs/read_qc_HN1_full_20261030';bam=source/'aligned.bam';reference=read_fasta(source/'genome.fna');samtools=ROOT/'tools/readqc-env/bin/samtools';regions=[]
 for row in csv.DictReader((source/'region_depth_summary.csv').open()):
  if row['filter']=='mq20_bq0':regions.append(dict(name=row['name'],start=int(row['start']),end=int(row['end']),role=row['role']))
 summaries=[];commands=[];begin=time.monotonic()
 for index,region in enumerate(regions):
  seqid='LR813493.2';start,end=region['start'],region['end'];cmd=[str(samtools),'view',str(bam),f'{seqid}:{max(1,start-1000)}-{end+1000}'];commands.append(cmd)
  p=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=60);(RUN/f'{index}.sam').write_text(p.stdout)
  report=analyze_sam(p.stdout,reference[seqid],seqid,start,end);report.update(name=region['name'],role=region['role'],assembly='GCA_903819125.2',run='ERR4194695',source_bam=str(bam.relative_to(ROOT)),source_scope='full run, whole-file MD5 verified; exact read use in assembly version unverified',analysis_elapsed_seconds=time.monotonic()-begin)
  write_json(OUT/f'region_{index}.json',report)
  summaries.append(dict(region_file=f'region_{index}.json',name=region['name'],role=region['role'],signals=report['signals'],quality_information_unresolved=report['quality_information_unresolved'],metrics=report['metrics'],coherent_substitution_count=len(report['coherent_substitutions']),coherent_indel_count=sum(r['coherent'] for r in report['indels']),mixed_allele_count=len(report['mixed_alleles'])))
  print(region['name'],report['signals'],report['metrics']['mismatch_fraction'],flush=True)
 write_json(OUT/'real_region_summary.json',dict(regions=summaries,elapsed_seconds=time.monotonic()-begin,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),causal_diagnosis='unresolved'))
 write_json(OUT/'real_input_manifest.json',dict(bam=str(bam.relative_to(ROOT)),bam_sha256=sha256(bam),reference_sha256=sha256(source/'genome.fna'),commands=commands,samtools_version=subprocess.check_output([str(samtools),'--version'],text=True).splitlines()[0]))
if __name__=='__main__':main()
