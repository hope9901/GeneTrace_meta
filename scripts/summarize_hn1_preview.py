#!/usr/bin/env python3
import json,pathlib,shutil
root=pathlib.Path(__file__).resolve().parents[1];src=root/'runs/read_qc_HN1_prefix_20261030';out=root/'outputs/read_qc_HN1_prefix_20261030';out.mkdir(parents=True,exist_ok=True)
r=json.loads((src/'report.json').read_text())
for name in ['report.json','region_depth_summary.csv','regional_read_support.png','fastq_qc.json','inputs.json','tool_versions.json','prefix_provenance.json']:
 shutil.copy2(src/name,out/name)
text=f'''# HN1 YJU3 raw-read evidence — prefix preview, 2026-10-30

## Finding

The previously unannotated YJU3 homolog in GCA_903819125.2 at {r['region']} has positive support from candidate PacBio run ERR4194695. This is provisional evidence of sequence presence and strengthens the annotation-gap hypothesis. It does not establish function, biological gene loss, an error-free assembly, or exact run use in this assembly version.

Only the first {r['fastq']['reads']:,} complete FASTQ records ({r['fastq']['bases']:,} bases) were examined. This prefix is **not a random subsample**, is incomplete and has **not passed the full-run MD5 check**. Its observed coverage is not an estimate of final sequencing coverage. Do not interpret missing or uneven prefix coverage as a sequencing/assembly defect.

| Candidate interval evidence | Result |
| --- | --- |
| Interval length | 939 bp |
| Covered bases, MAPQ >=20 / BQ >=0 | 100% |
| Mean aligned base depth in this preview | {r['candidate_depth']['mean_depth']:.2f} |
| Minimum aligned base depth | {r['candidate_depth']['min_depth']} |
| MAPQ >=20 alignments spanning both extra 100bp flanks | {r['mq20_reads_spanning_plus100bp']} read records |
| Forward / reverse spanning read records | {r['spanning_forward']} / {r['spanning_reverse']} |
| Coverage with BQ >=10 instead | 0% |

All {r['fastq']['bases']:,} observed FASTQ quality characters were `!` (encoded Phred 0). Their original meaning remains unresolved. Treating them as informative quality scores removes every observed base at BQ>=10 and creates an apparent zero-coverage region. This observation alone does not mean poor sequencing: original submission or conversion metadata must be checked. Mapping quality is a separate alignment confidence measure; the regional comparison uses BQ>=0 and records that decision explicitly.

Archive-generated FASTQ names do not preserve original PacBio ZMW identifiers. Therefore the {r['mq20_reads_spanning_plus100bp']} spanning records must not be described as that many independent molecules. Per-base consensus, small indels and functional integrity remain untested.

## Provenance

The paper's Table S5 associates HN1 with PacBio project PRJEB38713 and reused Illumina project PRJNA396809. The corresponding run candidates are ERR4194695 and SRR7851909. JXXY16.1 candidates are ERR4194697 and SRR7851726. Assembly and raw reads use different BioSample accessions. Across eight study assemblies, locally downloaded nuclear and mitochondrial lengths match Table S1 exactly, corroborating the study/assembly association; matching lengths are not sequence identity proof or an exact run/version provenance record.

See `../read_dataset_20261030/publication_supported_read_links.csv` (14 candidate run links across eight assemblies) and `publication_assembly_metric_checks.csv` (8/8 exact nuclear/mitochondrial length matches). Short-read lookup for the other six matching strains uses exact sample aliases, avoiding HN1/HN18 confusion. Short-read projects for Y55 and EM14S01-3B remain to be checked.

Sources: [paper](https://doi.org/10.1093/gbe/evab001), [supplementary Tables S1/S5](https://www.pure.ed.ac.uk/ws/portalfiles/portal/480257147/evab001_supplementary_data.pdf), saved metadata reports under `data/read_dataset_20261030/`.

## Reproducibility and pending work

`runs/read_qc_HN1_prefix_20261030/` retains the FASTQ prefix, input hashes, alignments, per-position depths, commands and tool versions. Preset `map-pb` is used for the PacBio subreads file; no HiFi status is inferred. Secondary, supplementary, unmapped, QC-failed and duplicate-flagged alignments are excluded from the depth comparison. Four neighboring annotated gene intervals provide local controls in `region_depth_summary.csv`.

The full 927,324,001-byte FASTQ download continues in resumable chunks. After whole-file MD5 verification, `scripts/finish_hn1_read_qc.py` automatically runs full-read regional mapping and exports to `outputs/read_qc_HN1_20261030/`. Latest state: `pipeline_status.json` in that directory. Full-run results must be reviewed before making a full-depth assessment. Genome-wide orthogroup inference continues separately.

The interpretation guards have five new tests: all-zero quality encoding does not prove sequencing failure; incomplete or unverified inputs cannot support low-coverage diagnosis; positive spanning reads support candidate sequence presence only.
'''
(out/'report.md').write_text(text);(root/'docs/HN1_RAW_READ_PREVIEW.md').write_text(text)
print('Preview report and evidence exported')
