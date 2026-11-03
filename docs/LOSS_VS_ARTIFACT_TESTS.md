# First tests: biological loss versus data/annotation artifacts

Date: 2026-10-30. The full 279-genome OrthoFinder computation remains running separately. None of its source or prepared protein inputs were changed by these tests.

## What can be claimed now

The data contain real examples where an apparent protein/gene-content difference cannot be interpreted as biological sequence loss. Within-sample assembly contrasts, annotation-scope audits, direct sequence searches and neighborhood evidence support this conclusion. **The assembly-level boxplots themselves do not identify a biological-loss mechanism or an independent causal effect of assembly level.** No functional gene losses were confirmed.

The supported research claim is: public assembly representation, annotation scope and search settings can alter observed gene-content estimates; sequence and neighborhood checks can identify specific misleading absence signals and separate them from unresolved disruption candidates.

## 1. Same registered biological sample, alternative assemblies

The pre-deduplication inventory was reused only as a paired control, not reintroduced as independent genomes into the 279-genome analysis. One record per BioSample/assembly-level was selected deterministically, avoiding double-counting duplicated deposits.

- 23 BioSamples have a Contig–Scaffold pair.
- All 23 pairs share the reported sequencing group and annotation provider.
- The median absolute difference in protein FASTA records is **772**.
- Scaffold has more records in 10 pairs and fewer in 13; median Scaffold-minus-Contig is −260. There is no consistent direction that would justify treating assembly level as a simple biological ranking.
- Largest contrast: SAMEA12946385, GCA_949123435.1 (`ANM_HP.nuclear_genome.ScRAP`, Contig) has 57.099456 Mb and 27,828 protein records; GCA_949124205.1 (`ANM_collapsed.nuclear_genome.ScRAP`, Scaffold) has 11.949778 Mb and 5,893 records. These are alternative representations of the same registered sample, not independent strains showing evolutionary gene expansion/loss. HP terminology is not used as an independently verified ploidy measurement.

These observations are not a randomized assembly experiment, nor proof that either representation is wrong. Completeness, phasing/collapsing and annotation can contribute. Same BioSample is metadata linkage, not sequence-verified clone identity.

Outputs: `same_biosample_cross_level_pairs.csv`, `same_biosample_pair_summary.json`, `same_sample_assembly_pairs.png/.pdf`.

## 2. Whole-cohort scope and confounding audit

All 279 retained assembly FASTA/GFF/gene maps were audited. GC measurements were checked against the previously computed whole-FASTA values. Explicitly labeled organellar/plasmid sequences were separated; unlabeled contigs were retained as nuclear-or-unclassified, not guessed to be organellar.

- 116 assemblies include explicitly labeled organellar sequence.
- JXXY16.1, GCA_903819185.2, includes **11,692,380 bp of nuclear sequence**, but has **zero submitted nuclear coding loci and 16 mitochondrial protein records**. The original GFF contains nuclear `region` entries, not nuclear gene models. This is an annotation-scope mismatch, not evidence that the nuclear genes disappeared.
- Of 34 annotation-provider × reported-sequencing strata, only 2 span multiple assembly levels. Most level contrasts are confounded with data source/technology. Within-stratum matching does not fully isolate assembly level, and provider is not asserted to be the actual software pipeline.
- Complete Genome versus Chromosome median GC gap is about **0.303 percentage points** using all FASTA records and **0.319 points** after excluding explicitly labeled non-nuclear sequence. Organelle inclusion changes measurements but does **not** remove this group difference.
- After strain deduplication, protein medians remain 5,556.5 / 5,412 / 5,725 / 5,712 for Complete Genome / Chromosome / Scaffold / Contig. Remaining differences are not automatically biological or technical; lineage, representation and annotation remain entangled.

Outputs: `annotation_scope_audit.csv`, `assembly_level_scope_summary.csv`, `provider_sequencing_strata.csv`, `gc_scope_comparison.png/.pdf`.

## 3. Direct tests on five unmodified assemblies

A fixed panel was chosen **before examining target search results**: four coordinate-quantile loci per S288C nuclear chromosome, yielding 64 distinct reference genes. They have a single reference CDS, 150–800 aa and no X/internal-stop residues. This controls assay complexity and does not represent all intron-rich, short or highly divergent genes. The full all-gene orthology run is unchanged.

Each reference protein was searched against both the submitted protein FASTA (BLASTP) and the complete genome (TBLASTN). Strong-hit criteria: identity ≥90%, single-HSP query coverage ≥90%, E-value ≤1e-20. Results are query-homology counts, **not independently confirmed ortholog counts or functional-gene counts**.

| Assembly / strain | Submitted protein records | Queries positive in protein FASTA | Queries with strong genomic hit | Additional partial/divergent evidence |
|---|---:|---:|---:|---:|
| GCA_000146045.2 / S288C reference | 6,002 | 64 | 64 | 0 |
| GCA_903819185.2 / JXXY16.1 | 16 | 0 | 63 | 1 |
| GCA_903819125.2 / HN1 | 5,722 | 63 | 64 | 0 |
| GCA_000190255.1 / FostersB | 3,689 | 30 | 38 | 20 |
| GCA_000326005.1 / FostersO | 3,575 | 32 | 37 | 23 |

The cases were chosen as diagnostic extremes plus reference/source-matched controls. They do not estimate artifact prevalence across all 279 genomes. Multiple queries may recognize homologous loci; do not count them as distinct functional genes without orthogroup and locus checks.

### Calibration and search-setting sensitivity

With default SEG/composition correction, IES6 and ILM1 did not pass 90% coverage even against their own reference sequences. Their 100%-identity alignments covered about 84% and 83% of the queries. Turning SEG and composition correction off recovered both. The sensitivity mode recovered all 64 reference queries in both assays.

Final positive evidence is the union of settings that pass the reference control for that query. Partial-hit evidence uses ≥70% identity, ≥50% single-HSP coverage, E≤1e-20 and remains unresolved. No absent-call is made from failure to pass a threshold. Removing masking can increase false-positive risk; strict thresholds, independent coordinates and neighborhood checks are retained, and homology is not declared to be orthology/function.

JXXY16.1's remaining strict-negative reference query is MRPL10. It has a **98.846% identity, 80.745% coverage** genomic hit. It is a partial/coverage-limited homolog candidate, not an established deletion or intact functional gene.

Outputs: `adjudicated_case_summary.csv`, `adjudicated_locus_results.csv`, `search_calibration_summary.json`, `protein_vs_genome_evidence.png/.pdf`. Default and sensitivity outputs are retained separately. The early `natural_case_reports.json` and per-case raw diagnostic files are preliminary search-only observations; the adjudicated files add gene-model context and supersede their candidate classification labels.

## 4. Neighborhood and disrupted-model checks

Two reference flanking coding loci were searched independently for nine selected candidates. A supported neighborhood requires both anchors on the target contig, each on its expected side and orientation, within 25 kb, identity ≥70%, coverage ≥50% and E≤1e-20. Four of nine passed both-anchor checks. Failure of this check is unresolved, not gene-loss evidence.

| Candidate | Evidence | Current interpretation |
|---|---|---|
| JXXY16.1 ACS1 | 99.439% identity, 100% query coverage; PEX22/FLC2 flanks support location; no gene model at locus | Unannotated syntenic homolog; intact function untested |
| JXXY16.1 PSY4 | 95.556% identity, 100% coverage; EDE1/COR1 flanks support location; no gene model | Unannotated syntenic homolog; intact function untested |
| HN1 YJU3 | 99.042% identity, 100% coverage at LR813493.2:287238–288176; YJU2/MBR1 flanks; no overlapping gene/CDS | Strong annotation-gap candidate in an otherwise extensively annotated genome |
| FostersB SWC3 | 91.401% identity, 100% coverage; CYS3/MDM10 flanks; original gene Note says `CDS interrupted by frameshift`, with six `heterozygous` sequence-alteration records | **Disrupted-model candidate, not a simple annotation omission.** Biological pseudogenization, assembly/base error and haplotype handling remain unresolved |

Across the 64-query pilot, 3 FostersB and 2 FostersO strong-DNA/protein-negative cases have explicit disruption/pseudogene-related model evidence. Those are classified separately rather than mislabeled as unannotated genes. All overlapping gene models and structural notes are retained in the adjudicated CSV.

The strongest annotation-gap candidates retain homologous DNA and neighborhood context, but no expression/protein-function assay or raw-read confirmation was performed. True functional loss can occur without physical deletion, so DNA homology alone cannot rule out functional loss.

Outputs: `synteny_candidate_tests.csv/.json`. Exact commands, source hashes, full alignments and assessment records live under `runs/loss_artifact_tests_20261030/`.

## 5. Safeguards added for the future agent

- Check whether the relevant compartment is supplied and annotated before interpreting group zeros.
- Run a self/reference control and search-setting sensitivity check before trusting a threshold-negative result.
- Inspect all overlapping gene models and structural notes, not only the presence of a CDS.
- Separate unannotated homologs, disrupted models and no/partial-hit unresolved cases.
- Add two-flank neighborhood evidence where possible; wait for full orthogroup evidence for difficult paralog cases.
- Reserve biological-loss claims for adequate locus localization, sequence/ORF inspection, raw-read or independent-assembly support and an explicit distinction between physical deletion and functional inactivation.

The background grouping input is unchanged. Its eventual postprocessor now adds `presence_absence_scope_aware.tsv` **alongside** the original observed matrices. Missing nuclear annotation or missing supplied organelle scope can produce N/A rather than a misleading zero. A retained zero still means unobserved annotation membership, not validated biological absence.

## Reproduce

```sh
.venv/bin/python scripts/audit_annotation_scope.py
python3 scripts/test_natural_annotation_gaps.py
python3 scripts/test_natural_annotation_gaps.py --sensitivity
.venv/bin/python scripts/summarize_loss_artifact_tests.py
python3 scripts/test_candidate_synteny.py
.venv/bin/python scripts/test_same_sample_assembly_pairs.py
.venv/bin/python -m unittest discover -s tests -p 'test_loss_artifact_evidence.py' -v
```

Existing successful natural-case searches are cached only for the same panel SHA-256. Use a new run path to change thresholds or the query panel. Do not overwrite an active OrthoFinder proteome.
