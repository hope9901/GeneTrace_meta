#!/usr/bin/env python3
"""Deterministic, labeled simulation benchmark; not a natural-error validation set."""
import collections,csv,datetime,json,pathlib,random,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff,write_fasta,write_json,sha256
from genetrace.regional_qc import analyze_sam,POLICY
OUT=ROOT/'outputs/regional_validation_20261102';RUN=ROOT/'runs/regional_validation_20261102/benchmark'
SCENARIOS=['normal','annotation_removed','reference_substitution','reference_indel','assembly_gap','misjoin','noisy_reads','dropout','mixed_alleles']
EXPECTED={'normal':'read_assembly_supported','annotation_removed':'annotation_gap_candidate','reference_substitution':'coherent_substitution_disagreement','reference_indel':'coherent_indel_disagreement','assembly_gap':'assembly_gap_or_ambiguity','misjoin':'junction_discontinuity_signal','noisy_reads':'diffuse_read_disagreement','dropout':'insufficient_read_support','mixed_alleles':'mixed_allele_signal'}
def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def other(c):return next(b for b in 'ACGT' if b!=c)
def windows():
 result=[];rng=random.Random(20261102)
 for source,acc in [('data/orthology/source/GCA_903819125.2','GCA_903819125.2'),('data/bacteria/source/GCF_000005845.2','GCF_000005845.2')]:
  sequences=read_fasta(ROOT/source/'genome.fna.gz');features=read_gff(ROOT/source/'annotation.gff.gz');eligible=[f for f in features if f['type']=='CDS' and f['end']-f['start']>=600 and (f['start']+f['end'])//2>3100 and (f['start']+f['end'])//2<len(sequences[f['seqid']])-3100];rng.shuffle(eligible);chosen=[]
  for f in eligible:
   center=(f['start']+f['end'])//2
   if any(seqid==f['seqid'] and abs(center-c)<6500 for seqid,c in chosen):continue
   seq=sequences[f['seqid']][center-3000:center+3000]
   if set(seq)-set('ACGT'):continue
   chosen.append((f['seqid'],center));result.append(dict(assembly=acc,source_seqid=f['seqid'],source_start=center-2999,source_end=center+3000,source_gene=f['attrs'].get('gene',f['attrs'].get('locus_tag','N/A')),sequence=seq,split='calibration' if len(chosen)<=3 else 'holdout'))
   if len(chosen)==12:break
  if len(chosen)!=12:raise ValueError('Insufficient nonoverlapping source loci')
 return result

def reads(dna,scenario,rng):
 out=[]
 for index,start in enumerate(range(0,5001,50)):
  if scenario=='dropout' and start<3200 and start+1000>2800:continue
  for strand in ['+','-']:
   seq=list(dna[start:start+1000]);rate=.12 if scenario=='noisy_reads' else .003
   for i,c in enumerate(seq):
    if rng.random()<rate:seq[i]=rng.choice([b for b in 'ACGT' if b!=c])
   if scenario=='mixed_alleles' and index%2==0 and start<=3000<start+1000:seq[3000-start]=other(dna[3000])
   seq=''.join(seq);seq=seq if strand=='+' else reverse(seq)
   out.append((f'r{index:03d}_{strand}',seq,('#' if scenario=='noisy_reads' else 'D')*len(seq)))
 return out

def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True);cases=[];start_time=time.monotonic();mm=ROOT/'tools/readqc-env/bin/minimap2';sources=windows();commands=[]
 write_json(OUT/'benchmark_design.json',dict(seed=20261102,scenarios=SCENARIOS,expected_signals=EXPECTED,policy_frozen_before_run=POLICY,loci=[{k:v for k,v in w.items() if k!='sequence'} for w in sources],read_model='1000bp single segments, both orientations, 50bp start spacing, background independent substitution probability 0.003; noisy scenario 0.12; not a calibrated instrument simulator',split='3 calibration + 9 holdout nonoverlapping loci per organism; no threshold fitting on either split'))
 for wi,window in enumerate(sources):
  dna=window['sequence']
  for si,scenario in enumerate(SCENARIOS):
   case=f'locus{wi:02d}_{scenario}';directory=RUN/case;directory.mkdir(parents=True,exist_ok=True);reference=dna
   if scenario=='reference_substitution':reference=dna[:3000]+other(dna[3000])+dna[3001:]
   elif scenario=='reference_indel':reference=dna[:3000]+'TACGTAC'+dna[3000:]
   elif scenario=='assembly_gap':reference=dna[:2950]+'N'*100+dna[3050:]
   elif scenario=='misjoin':reference=dna[:3000]+dna[3800:]
   records=reads(dna,scenario,random.Random(20261102+wi*100+si));write_fasta(directory/'reference.fna',{'region':reference})
   with (directory/'reads.fastq').open('w') as f:
    for name,seq,qual in records:f.write('@'+name+'\n'+seq+'\n+\n'+qual+'\n')
   gff='##gff-version 3\n' if scenario=='annotation_removed' else '##gff-version 3\nregion\tsimulated\tCDS\t2850\t3150\t.\t+\t0\tID=control-CDS\n'
   (directory/'annotation.gff3').write_text(gff)
   cmd=[str(mm),'-ax','map-ont','--MD','-t','1','--secondary=no',str(directory/'reference.fna'),str(directory/'reads.fastq')]
   signature=dict(reference_sha256=sha256(directory/'reference.fna'),reads_sha256=sha256(directory/'reads.fastq'),mapping_preset='map-ont --MD --secondary=no')
   cached=json.loads((directory/'alignment_inputs.json').read_text()) if (directory/'alignment_inputs.json').exists() else None
   if not (directory/'alignments.sam').exists() or cached!=signature:
    with (directory/'alignments.sam').open('w') as so,(directory/'minimap2.stderr').open('w') as se:subprocess.run(cmd,stdout=so,stderr=se,check=True,timeout=60)
    write_json(directory/'alignment_inputs.json',signature)
   report=analyze_sam((directory/'alignments.sam').read_text(),reference,'region',2901,3100)
   if not read_gff(directory/'annotation.gff3') and report['signals']==['read_assembly_supported']:report['signals']=['annotation_gap_candidate']
   priority=['assembly_gap_or_ambiguity','insufficient_read_support','junction_discontinuity_signal','coherent_indel_disagreement','coherent_substitution_disagreement','diffuse_read_disagreement','mixed_allele_signal','annotation_gap_candidate','read_assembly_supported']
   predicted=next(s for s in priority if s in report['signals']);expected=EXPECTED[scenario]
   report.update(scenario=scenario,expected_signal=expected,synthetic_control=True);write_json(directory/'report.json',report)
   cases.append(dict(case=case,locus=wi,split=window['split'],source_assembly=window['assembly'],scenario=scenario,expected=expected,predicted=predicted,expected_signal_found=expected in report['signals'],mismatch_fraction=report['metrics']['mismatch_fraction'],mean_depth=report['metrics']['mean_depth'],spanning_records=report['metrics']['spanning_plus100_records'],boundary_clipped_fraction=report['metrics']['boundary_clipped_fraction'],coherent_snp_count=len(report['coherent_substitutions']),coherent_indel_count=sum(r['coherent'] for r in report['indels']),mixed_count=len(report['mixed_alleles']),biological_loss_inferred=False))
   commands.append(dict(case=case,command=cmd,reference_sha256=sha256(directory/'reference.fna'),reads_sha256=sha256(directory/'reads.fastq'),annotation_sha256=sha256(directory/'annotation.gff3')))
  write_json(OUT/'benchmark_status.json',dict(stage='running',completed_loci=wi+1,total_loci=len(sources),completed_cases=len(cases),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));print(f'{wi+1}/{len(sources)} loci / {len(cases)} cases',flush=True)
 with (OUT/'benchmark_cases.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(cases[0]));w.writeheader();w.writerows(cases)
 groups=[]
 for split in ['calibration','holdout']:
  for scenario in SCENARIOS:
   group=[r for r in cases if r['split']==split and r['scenario']==scenario];groups.append(dict(split=split,scenario=scenario,cases=len(group),expected_signal_detected=sum(r['expected_signal_found'] for r in group),primary_label_correct=sum(r['predicted']==r['expected'] for r in group)))
 write_json(OUT/'benchmark_summary.json',dict(cases=len(cases),loci=len(sources),elapsed_seconds=time.monotonic()-start_time,groups=groups,holdout_primary_accuracy=sum(r['predicted']==r['expected'] for r in cases if r['split']=='holdout')/sum(r['split']=='holdout' for r in cases),holdout_expected_signal_recall=sum(r['expected_signal_found'] for r in cases if r['split']=='holdout')/sum(r['split']=='holdout' for r in cases),natural_error_detection_performance='not_established',limitations=['Simulation-only pattern recovery; not causal diagnosis or real-world clinical/research performance.','Single-segment read model, no paired-end library, no realistic CLR indel error model.','Small non-overlapping locus sample; no repetitive genome-wide competitor sequences in the mapping reference.','Control labels are hidden from the signal analyzer; annotation inspection is a separate layer.']))
 write_json(OUT/'benchmark_commands.json',dict(minimap2_version=subprocess.check_output([str(mm),'--version'],text=True).strip(),commands=commands));write_json(OUT/'benchmark_status.json',dict(stage='complete',cases=len(cases)))
if __name__=='__main__':main()
