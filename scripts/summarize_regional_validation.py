#!/usr/bin/env python3
import collections,csv,json,pathlib,sys,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/genetrace-matplotlib');os.environ.setdefault('XDG_CACHE_HOME','/tmp/genetrace-cache')
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.regional_qc import analyze_sam
from genetrace.io import read_fasta,write_json,sha256
OUT=ROOT/'outputs/regional_validation_20261102';RUN=ROOT/'runs/regional_validation_20261102'
def main():
 real=json.loads((OUT/'real_region_summary.json').read_text());bench=json.loads((OUT/'benchmark_summary.json').read_text());cases=list(csv.DictReader((OUT/'benchmark_cases.csv').open()))
 # Identical observations can arise from different causes; demonstrate the limit explicitly.
 d=RUN/'benchmark/locus00_reference_substitution';text=(d/'alignments.sam').read_text();reference=read_fasta(d/'reference.fna')['region'];report=analyze_sam(text,reference,'region',2901,3100)
 write_json(OUT/'causal_ambiguity_demo.json',dict(observed_sam_sha256=sha256(d/'alignments.sam'),same_observations_hypotheses=['The reference/assembly base is wrong for the sampled genome.','Reads come from a genuine variant sample compared with another reference.','A systematic sequencing/conversion error produces the same alternative base.'],observed_signals=report['signals'],causal_diagnosis=report['causal_diagnosis'],distinguished_by_this_evidence=False,additional_evidence_needed=['Exact run-to-assembly and sample provenance','Independent sequencing library/technology','Original quality metadata, local repeat and haplotype context']))
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;import numpy as np
 scenarios=list(dict.fromkeys(r['scenario'] for r in cases));groups=[g for g in bench['groups'] if g['split']=='holdout'];fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
 axes[0].barh(range(len(groups)),[g['expected_signal_detected'] for g in groups],color='#317362');axes[0].set_yticks(range(len(groups)),[g['scenario'] for g in groups],fontsize=9);axes[0].invert_yaxis();axes[0].set_xlim(0,20);axes[0].set_xlabel('Expected pattern recovered / 18 held-out loci');axes[0].set_title('Simulation controls (not real-error accuracy)')
 for i,g in enumerate(groups):axes[0].text(g['expected_signal_detected']+.2,i,f"{g['expected_signal_detected']}/{g['cases']}",va='center',fontsize=9)
 rr=real['regions'];axes[1].bar(range(len(rr)),[r['metrics']['mismatch_fraction']*100 for r in rr],color=['#c18a48']+['#6e9f8b']*4);axes[1].set_xticks(range(len(rr)),['YJU3']+[r['name'].replace('PACBIOSEQ_LOCUS','locus ') for r in rr[1:]],rotation=25,ha='right');axes[1].set_ylabel('Local substitution mismatch (%)');axes[1].set_title('HN1 full PacBio run / MQ>=20, BQ>=0');axes[1].set_ylim(0,max(r['metrics']['mismatch_fraction']*100 for r in rr)*1.3)
 fig.savefig(OUT/'validation_overview.png',dpi=170);plt.close(fig)
 candidate=rr[0];m=candidate['metrics'];normal=[r for r in cases if r['split']=='holdout' and r['scenario']=='normal'];fp=sum(r['predicted']!='read_assembly_supported' for r in normal);matched=sum(g['expected_signal_detected'] for g in groups);all_test=sum(g['cases'] for g in groups)
 text=f'''# Regional read disagreement and control validation — 2026-11-02

## Implemented

The CIGAR parser compares SAM SEQ directly to reference bases within the requested 1-based region. Whole-read NM is not reused as a local error rate. Primary MQ>=20 alignments are examined with BQ>=0; MAPQ 255 is excluded as unavailable. Supplementary alignments are counted separately, not counted as independent primary support. M, =, X, I, D, N, S, H and P operations are handled explicitly. Reference skips are not deletions, and reverse-strand SAM SEQ is not reverse-complemented again.

Per-base support, bidirectional substitutions, aligned insertion/deletion representations, soft/hard clipping near boundaries and actual aligned support at both outer anchors are recorded. Indel denominators require aligned bases at both event flanks; incomplete read-end events are excluded from coherent support. Indels are NOT left-normalized or reported as variant calls.

Flags use provisional thresholds: depth>=10, alternate records>=5, support fraction>=0.8, alternate records on both strands>=2. Mixed support>=0.2 but<0.8 is a separate flag. Clipping and diffuse mismatch flags are investigation cues, not causal diagnoses. These thresholds are not universally calibrated across sequencing technologies.

## HN1 YJU3 and neighboring controls

Five intervals were checked using the existing full-run BAM. YJU3 has mean aligned depth {m['mean_depth']:.2f}, local substitution mismatch {m['mismatch_fraction']*100:.3f}%, and {m['spanning_plus100_records']} primary records covering both extra 100bp anchors without an internal D/N gap>=50bp. This is a stricter span definition than the earlier 141-record bounding-coordinate check; the counts must not be interpreted as a loss of read data.

The four neighboring control intervals have local substitution mismatch {min(r['metrics']['mismatch_fraction'] for r in rr[1:])*100:.3f}%–{max(r['metrics']['mismatch_fraction'] for r in rr[1:])*100:.3f}%. All five intervals have `read_assembly_supported` as the current flag, without substitutions or indels meeting the coherent-disagreement rule. Random per-read indel burden is still present. This does NOT establish an error-free assembly or intact gene function.

All original FASTQ qualities encode `!` and independent PacBio ZMW identity is unavailable. Exact run use in this assembly version is unverified. Accordingly no variant calls, sequencing-failure diagnosis, biological-loss or definitive assembly-error calls are made. The result strengthens the local sequence-presence evidence underlying the YJU3 annotation-gap candidate.

## Labeled controls

Twenty-four non-overlapping 6kb windows (12 from HN1 and 12 from E. coli MG1655) were sampled with fixed seed20261102. Per organism three windows were assigned to a development/calibration split and nine to holdout. Thresholds were fixed before the benchmark; neither split was used to fit them. Nine scenarios yield216 cases: normal, annotation removal, reference substitution, reference insertion, assembly N gap, misjoin, diffuse read substitution noise, read dropout, and mixed alleles.

The held-out set recovered the expected signal in **{matched}/{all_test}** cases. Normal held-out controls had **{fp}/{len(normal)}** flagged deviations. These are simulation pattern-recovery results, NOT an estimate of accuracy or sensitivity on natural sequencing/assembly errors. The reads are1000bp single segments with simple independent substitutions, no paired-end library or realistic CLR indel error model. Alignment references contain a single cropped window, without genome-wide repeat competitors. The case-level table records extra signals and primary labels for audit.

The causal ambiguity demonstration shows that the same coherent alternate-base pattern can be compatible with a wrong assembly/reference base, a real sample variant, or systematic read error. This evidence alone does not distinguish them. Exact sample provenance, another library/technology, quality and repeat/haplotype context are required.

## Artifacts and reproduction

- `real_region_summary.json`, `region_0.json`–`region_4.json`: actual HN1 metrics and per-base/read observations.
- `benchmark_design.json`, `benchmark_cases.csv`, `benchmark_summary.json`: fixed design, all216 cases and split-specific results.
- `causal_ambiguity_demo.json`: what local alignments cannot decide.
- `validation_overview.png`: held-out recovery and real mismatch comparison.
- `real_input_manifest.json`, `benchmark_commands.json`: source checksums, commands and versions.

```sh
.venv/bin/python scripts/analyze_regional_disagreement.py
.venv/bin/python scripts/benchmark_regional_qc.py
.venv/bin/python scripts/summarize_regional_validation.py
```

Parsing follows the [SAM specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.tex). Original genomes, annotations and the full BAM are not modified. No additional whole-genome alignment is needed for the real HN1 audit.
'''
 (OUT/'report.md').write_text(text);(ROOT/'docs/REGIONAL_VALIDATION.md').write_text(text);print('Final report:',matched,'/',all_test,'holdout expected patterns; normal flags',fp,'/',len(normal))
if __name__=='__main__':main()
