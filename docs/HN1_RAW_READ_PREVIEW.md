# HN1 YJU3 raw-read evidence — prefix preview, 2026-10-30

## Finding

The previously unannotated YJU3 homolog in GCA_903819125.2 at LR813493.2:287238-288176 has positive support from candidate PacBio run ERR4194695. This is provisional evidence of sequence presence and strengthens the annotation-gap hypothesis. It does not establish function, biological gene loss, an error-free assembly, or exact run use in this assembly version.

Only the first 20,000 complete FASTQ records (83,047,354 bases) were examined. This prefix is **not a random subsample**, is incomplete and has **not passed the full-run MD5 check**. Its observed coverage is not an estimate of final sequencing coverage. Do not interpret missing or uneven prefix coverage as a sequencing/assembly defect.

| Candidate interval evidence | Result |
| --- | --- |
| Interval length | 939 bp |
| Covered bases, MAPQ >=20 / BQ >=0 | 100% |
| Mean aligned base depth in this preview | 19.19 |
| Minimum aligned base depth | 12 |
| MAPQ >=20 alignments spanning both extra 100bp flanks | 17 read records |
| Forward / reverse spanning read records | 8 / 9 |
| Coverage with BQ >=10 instead | 0% |

All 83,047,354 observed FASTQ quality characters were `!` (encoded Phred 0). Their original meaning remains unresolved. Treating them as informative quality scores removes every observed base at BQ>=10 and creates an apparent zero-coverage region. This observation alone does not mean poor sequencing: original submission or conversion metadata must be checked. Mapping quality is a separate alignment confidence measure; the regional comparison uses BQ>=0 and records that decision explicitly.

Archive-generated FASTQ names do not preserve original PacBio ZMW identifiers. Therefore the 17 spanning records must not be described as that many independent molecules. Per-base consensus, small indels and functional integrity remain untested.

## Provenance

The paper's Table S5 associates HN1 with PacBio project PRJEB38713 and reused Illumina project PRJNA396809. The corresponding run candidates are ERR4194695 and SRR7851909. JXXY16.1 candidates are ERR4194697 and SRR7851726. Assembly and raw reads use different BioSample accessions. Across eight study assemblies, locally downloaded nuclear and mitochondrial lengths match Table S1 exactly, corroborating the study/assembly association; matching lengths are not sequence identity proof or an exact run/version provenance record.

See `../read_dataset_20261030/publication_supported_read_links.csv` (14 candidate run links across eight assemblies) and `publication_assembly_metric_checks.csv` (8/8 exact nuclear/mitochondrial length matches). Short-read lookup for the other six matching strains uses exact sample aliases, avoiding HN1/HN18 confusion. Short-read projects for Y55 and EM14S01-3B remain to be checked.

Sources: [paper](https://doi.org/10.1093/gbe/evab001), [supplementary Tables S1/S5](https://www.pure.ed.ac.uk/ws/portalfiles/portal/480257147/evab001_supplementary_data.pdf), saved metadata reports under `data/read_dataset_20261030/`.

## Reproducibility and pending work

`runs/read_qc_HN1_prefix_20261030/` retains the FASTQ prefix, input hashes, alignments, per-position depths, commands and tool versions. Preset `map-pb` is used for the PacBio subreads file; no HiFi status is inferred. Secondary, supplementary, unmapped, QC-failed and duplicate-flagged alignments are excluded from the depth comparison. Four neighboring annotated gene intervals provide local controls in `region_depth_summary.csv`.

The full 927,324,001-byte FASTQ download continues in resumable chunks. After whole-file MD5 verification, `scripts/finish_hn1_read_qc.py` automatically runs full-read regional mapping and exports to `outputs/read_qc_HN1_20261030/`. Latest state: `pipeline_status.json` in that directory. Full-run results must be reviewed before making a full-depth assessment. Genome-wide orthogroup inference continues separately.

The interpretation guards have five new tests: all-zero quality encoding does not prove sequencing failure; incomplete or unverified inputs cannot support low-coverage diagnosis; positive spanning reads support candidate sequence presence only.
