# Domain profiles and bacterial transfer test — 2026-11-01

## Scope

All annotated proteins from 279 representative S. cerevisiae assemblies and four E. coli laboratory model strains were scanned with six selected Pfam HMM profiles: PF00005 (ABC_tran), PF00072 (Response_reg), PF00512 (HisKA), PF02518 (HATPase_c), PF00069 (Pkinase), PF00076 (RRM_1). This is a selected-panel test, not a whole-Pfam annotation. Current profile versions and SHA256 are in `data/domain_models/manifest.json`.

Bacterial genome, GFF, CDS and protein files were downloaded from NCBI and checked against the NCBI MD5 manifests. Strains are MG1655, W3110, REL606 and BL21(DE3), all complete assemblies. Their source reports and exact accessions are in `data/bacteria/cohort.json`.

## Method

HMMER 3.4 `hmmsearch --cpu 1 --cut_ga` searches each complete proteome. Both sequence and domain scores must meet profile-specific Pfam gathering thresholds. HMM coverage >=70% labels a near-complete domain; shorter GA-passing hits remain as partial evidence. Per-profile count tables distinguish proteins carrying the domain from domain instances; repeated domains in one protein are not extra genes. Architecture strings describe only these six profiles and are not full protein architectures.

For bacteria, all six reading frames are translated, split at stop codons and ORFs of at least 30 amino acids are searched with the same profiles. Internal codon translations follow table 11; these ORFs are not gene predictions and alternative start codons, complete CDS structure and function are not inferred. Hits are projected to 1-based genomic coordinates and compared with same-strand CDS. Pseudogene, frameshift and interrupted model records are inspected before suggesting annotation omission. Whole-genome profile searches have not been run on intron-containing yeast genomes.

## Findings and controls

MG1655 has 77 proteins carrying PF00005, with 94 domain instances. This illustrates why domain instance count and gene count must remain separate. A genomic PF00005 hit outside protein-coding annotation overlaps the annotated pseudogene **ycjV**: classify it as `pseudogene_or_disrupted_model_unresolved`, not an omitted functional gene.

The synthetic positive control removes MG1655 protein NP_416685.1 and its structural annotation from copied input files. Its PF00005 hit disappears from the protein search while the unchanged genome still contains the domain, recovering `unannotated_genomic_domain_candidate`. Original source files remain unchanged. This establishes recovery of one known annotation omission, not sensitivity to natural errors generally.

The existing whole-protein algorithm was also run on rpoB in all four strains: a full-length DNA hit at 100% amino-acid identity was detected in each, along with protein/structural annotation evidence. This is a conserved-gene positive control, not proof that every bacterial gene or artifact type is handled correctly.

JXXY16.1 has zero selected-profile hits in the supplied 16-protein organelle-only annotation. This does not establish nuclear domain or gene loss. Pkinase/RRM_1 profile absence in E. coli does not establish absence of every kinase or RNA-binding protein family; these profiles are specific families.

## Outputs

- `domain_counts.csv`: per-assembly, per-profile distinct protein / domain / partial counts.
- `domain_members.csv`: domain coordinates, scores, HMM coverage and protein IDs.
- `panel_architectures.csv`: order and repeats of the selected profiles on proteins.
- `genomic_domain_candidates.csv`: bacterial DNA coordinates and structural-annotation adjudication.
- `controls.json`: synthetic annotation-omission control outcome.
- `bacterial_domain_counts.png`: bacterial profile count comparison.

Raw HMMER tables, command/version/hash manifests, six-frame maps and control inputs are under `runs/domain_bacteria_20261101/`. No bacterial raw reads were downloaded; no sequencing error, biological loss or functional gene status is confirmed by these tests. Natural differences between the four assemblies require locus-specific read/synteny evidence and provenance checks.

References: [InterPro API](https://interpro-documentation.readthedocs.io/en/latest/api.html), [HMMER documentation](https://hmmer-web-docs.readthedocs.io/en/latest/searches.html), [NCBI MG1655 assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/).
