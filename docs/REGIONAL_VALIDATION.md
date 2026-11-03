# Regional read disagreement and control validation — 2026-11-02

## Implemented

The CIGAR parser compares SAM SEQ directly to reference bases within the requested 1-based region. Whole-read NM is not reused as a local error rate. Primary MQ>=20 alignments are examined with BQ>=0; MAPQ 255 is excluded as unavailable. Supplementary alignments are counted separately, not counted as independent primary support. M, =, X, I, D, N, S, H and P operations are handled explicitly. Reference skips are not deletions, and reverse-strand SAM SEQ is not reverse-complemented again.

Per-base support, bidirectional substitutions, aligned insertion/deletion representations, soft/hard clipping near boundaries and actual aligned support at both outer anchors are recorded. Indel denominators require aligned bases at both event flanks; incomplete read-end events are excluded from coherent support. Indels are NOT left-normalized or reported as variant calls.

Flags use provisional thresholds: depth>=10, alternate records>=5, support fraction>=0.8, alternate records on both strands>=2. Mixed support>=0.2 but<0.8 is a separate flag. Clipping and diffuse mismatch flags are investigation cues, not causal diagnoses. These thresholds are not universally calibrated across sequencing technologies.

## HN1 YJU3 and neighboring controls

Five intervals were checked using the existing full-run BAM. YJU3 has mean aligned depth 184.89, local substitution mismatch 2.587%, and 124 primary records covering both extra 100bp anchors without an internal D/N gap>=50bp. This is a stricter span definition than the earlier 141-record bounding-coordinate check; the counts must not be interpreted as a loss of read data.

The four neighboring control intervals have local substitution mismatch 2.404%–2.660%. All five intervals have `read_assembly_supported` as the current flag, without substitutions or indels meeting the coherent-disagreement rule. Random per-read indel burden is still present. This does NOT establish an error-free assembly or intact gene function.

All original FASTQ qualities encode `!` and independent PacBio ZMW identity is unavailable. Exact run use in this assembly version is unverified. Accordingly no variant calls, sequencing-failure diagnosis, biological-loss or definitive assembly-error calls are made. The result strengthens the local sequence-presence evidence underlying the YJU3 annotation-gap candidate.

## Labeled controls

Twenty-four non-overlapping 6kb windows (12 from HN1 and 12 from E. coli MG1655) were sampled with fixed seed20261102. Per organism three windows were assigned to a development/calibration split and nine to holdout. Thresholds were fixed before the benchmark; neither split was used to fit them. Nine scenarios yield216 cases: normal, annotation removal, reference substitution, reference insertion, assembly N gap, misjoin, diffuse read substitution noise, read dropout, and mixed alleles.

The held-out set recovered the expected signal in **162/162** cases. Normal held-out controls had **0/18** flagged deviations. These are simulation pattern-recovery results, NOT an estimate of accuracy or sensitivity on natural sequencing/assembly errors. The reads are1000bp single segments with simple independent substitutions, no paired-end library or realistic CLR indel error model. Alignment references contain a single cropped window, without genome-wide repeat competitors. The case-level table records extra signals and primary labels for audit.

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
