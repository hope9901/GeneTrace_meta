#!/usr/bin/env python3
"""Audit full overlapping annotated gene bounds, not just a translated HSP."""
import json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff,sha256
from genetrace.regional_qc import analyze_sam
OUT=ROOT/'outputs/ecoli100_raw_validation_20261102'
def main():
 p=OUT/'summary.json';summary=json.loads(p.read_text());metadata={r['assembly_accession']:r for r in json.loads((ROOT/'data/bacteria/cohort.json').read_text())}
 for assembly in summary['assemblies']:
  acc=assembly['assembly'];d=ROOT/'runs/ecoli100_raw_validation_20261102'/acc;source=ROOT/metadata[acc]['source_directory'];features=read_gff(source/'annotation.gff.gz');genome=read_fasta(source/'genome.fna.gz')
  for i,r in enumerate(assembly['regions']):
   if not r['role'].startswith('candidate'):continue
   overlap=[f for f in features if f['type'] in ['gene','pseudogene'] and f['seqid']==r['seqid'] and f['start']<=r['end'] and f['end']>=r['start']]
   if not overlap:continue
   start=min([r['start']]+[f['start'] for f in overlap]);end=max([r['end']]+[f['end'] for f in overlap]);cmd=[str(ROOT/'tools/readqc-env/bin/samtools'),'view',str(d/'aligned.bam'),f"{r['seqid']}:{max(1,start-1000)}-{end+1000}"];sam=subprocess.check_output(cmd,text=True);report=analyze_sam(sam,genome[r['seqid']],r['seqid'],start,end);report.update(name=r['name'],role='candidate_annotation_bounded_region',assembly=acc,run=assembly['run'],source_gff_sha256=sha256(source/'annotation.gff.gz'),command=cmd,original_HSP_span=[r['start'],r['end']],causal_diagnosis='unresolved',exact_assembly_run_use_verified=False)
   (OUT/(acc+'_'+str(i)+'.json')).write_text(json.dumps(report,indent=2)+'\n');r.update(start=start,end=end,role='candidate_annotation_bounded_region',metrics=report['metrics'],signals=report['signals'],coherent_substitution_count=len(report['coherent_substitutions']),coherent_indel_count=sum(x['coherent'] for x in report['indels']),quality_information_unresolved=report['quality_information_unresolved'])
 p.write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
