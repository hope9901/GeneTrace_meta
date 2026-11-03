#!/usr/bin/env python3
"""Attach provenance guard to existing reports without repeating alignment."""
import collections,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.read_provenance import assess_read_provenance
OUT=ROOT/'outputs/ecoli100_raw_validation_20261102'
def main():
 p=OUT/'summary.json';summary=json.loads(p.read_text())
 for assembly in summary['assemblies']:
  acc=assembly['assembly'];d=ROOT/'runs/ecoli100_raw_validation_20261102'/acc;r=assembly['regions'][0];cmd=[str(ROOT/'tools/readqc-env/bin/samtools'),'view',str(d/'aligned.bam'),f"{r['seqid']}:{max(1,r['start']-10000)}-{r['end']+10000}"];sam=subprocess.check_output(cmd,text=True);verification=json.loads((d/'fastq_verification.json').read_text());qualities=collections.Counter()
  for v in verification:qualities.update(v['quality_character_counts'])
  review=assess_read_provenance(sam,qualities);(d/'provenance_review.json').write_text(json.dumps(dict(**review,command=cmd),indent=2)+'\n');assembly['read_provenance']=review;assembly['independent_physical_read_validation_established']=False
  for path in OUT.glob(acc+'_*.json'):
   report=json.loads(path.read_text());report['read_provenance']=review;report['causal_diagnosis']='unresolved';path.write_text(json.dumps(report,indent=2)+'\n')
  print(acc,'provenance_review',review['provenance_review_required'],review['modal_start_step_fraction'],review['modal_template_length_fraction'])
 p.write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
