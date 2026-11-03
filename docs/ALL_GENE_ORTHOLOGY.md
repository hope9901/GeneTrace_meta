# All-gene, strain-deduplicated S. cerevisiae orthology cohort

This replaces the earlier small-gene pilot scope for cohort-wide analysis. All retained genomes and all available protein-coding loci participate. The ACT1 controls remain an independent validation fixture.

## Cohort and duplicate decisions

- Starting cohort: 435 current GenBank assembly accessions with GFF3, CDS and protein FASTA available. The original GenBank/RefSeq pair was already collapsed in the inventory.
- Queried all 283 non-empty NCBI BioSample accessions for additional strain/isolate metadata.
- Selected **279 representative assemblies**, recording 156 duplicate assembly removals.
- Merge evidence: same BioSample, exact strain/isolate labels after case/whitespace normalization, or explicit semicolon-separated strain aliases. Punctuation and derivative identifiers are retained. BY4741 and BY4742 are not collapsed, nor are distinct PE-2 haploid derivatives inferred to be the same strain.
- A representative is selected by: S288C reference anchor; Full representation; not extremely sparse in protein annotation; assembly level; lower N fraction; larger FASTA-sequence N50; later release date; accession as a deterministic tie breaker. Large genome/protein counts are not treated as evidence of higher quality.
- S288C partial and alternative records collapse to GCA_000146045.2. The explicit `CEN.PK113-7D; CBS 8340` alias is preserved in the audit.
- **121 retained records lack an explicit resolved strain/isolate attribute.** They retain separate BioSample/accession identities and titles. N/A labels are never merged with each other. Their BioSample titles often contain useful isolate labels, but titles are not used as unverified equivalence assertions.
- GCA_014905575.1 has conflicting assembly/BioSample strain labels; this is recorded and not used to bridge different samples.
- These are metadata-defined strain/sample units, not a guarantee of sequence-verified unique biological clones. Undocumented aliases can remain. Removed assembly metadata and every merge edge are retained for review.

Outputs: `outputs/orthology_20261030/strain_deduplication.csv`, `strain_deduplication_edges.json`, `selected_assemblies.csv`, `cohort_summary.json`.

## Complete source and gene inputs

All 279 source packages were downloaded: genome FASTA, GFF3, nucleotide CDS FASTA and protein FASTA. Each compressed file was checked against its NCBI-published MD5; source URLs and SHA-256 are retained in per-assembly manifests under `data/orthology/source/`.

- **1,541,431 original protein records**, retained without source modification.
- **1,541,009 orthology input coding units**, using the longest translated isoform per GFF-mapped coding locus.
- Zero unmapped protein units in this snapshot; zero missing-protein cases among the explicitly tagged protein-coding genes checked by the input audit.
- Distinct loci with identical peptides are retained as distinct genes. Sequences are not collapsed across genomes or gene loci.
- Only trailing `*` formatting is removed from representative FASTA; transformations and internal-stop flags are recorded. All original isoforms remain in the original FASTA and in the gene-to-protein correspondence table.
- Every exported sequence ID has a reversible link to assembly, gene ID, original protein IDs, coordinates, strand and isoform policy.
- This is **all available protein-coding genes**, not an assertion of orthology for noncoding RNA genes or genes absent from the submitted annotation. Unannotated DNA cannot enter protein orthology until models are recovered.
- GCA_903819185.2, with only 16 protein records, remains in the final cohort. It is flagged and excluded only from core initialization, never from the final comparison. Its apparent absences must not be interpreted as thousands of biological gene losses.

Source and intermediate directories are excluded from Git. Current source plus prepared FASTA and mappings use several GB. Input summary: `outputs/orthology_20261030/all_gene_input_summary.csv`.

## Orthology workflow and actual run state

The installed production workflow is **OrthoFinder 3.2**, with DIAMOND 2.2.8, MCL, FAMSA, FastTree and ASTRAL-Pro from a project-local dependency environment. The installed Python/package versions are recorded in `orthofinder_environment.json`.

OrthoFinder's official guidance recommends the scalable core/assign workflow for >100 species: https://github.com/OrthoFinder/OrthoFinder#advanced-usage---scaling-to-thousands-of-species . Here each proteome is one retained within-species assembly, not a separate biological species.

1. Use all coding loci of 64 initial core genomes to infer groups, gene trees and a rooted strain tree.
2. Add all coding loci of the other 215 genomes with official `--assign --core`.
3. Allow the official workflow to process genes not assigned to existing profiles; they are not silently discarded. Run the full assignment workflow, not an early `-og` stop.
4. Export final phylogenomic orthogroup membership, copy-number and presence/absence matrices.
5. Verify every input coding unit appears exactly once in the assigned or unassigned outputs, both globally and per assembly. Fail rather than publishing a complete-looking table if any unit is unaccounted for.

The initial core includes the S288C reference and uses greedy farthest selection on bottom-256 hashes of exact protein sequences. This is a transparent diversity-initialization heuristic, **not a phylogenetic distance estimate**. It does not limit the final cohort or select a small subset of genes. See `core_selection.json` for every chosen core and additional accession.

**The full 279-genome run was launched on 2026-10-30 at approximately 14:36 KST. It is a long-running computation; no full-cohort orthogroup counts are claimed in this document.** The live state is authoritative:

- `outputs/orthology_20261030/orthofinder_status.json`: current stage, timestamps, PID, cohort size and disk availability.
- `core_inference.log`: initial 64-genome inference.
- `all_genome_assignment.log`: created when full-cohort expansion starts.
- `scalable_run_state.json`: actual commands and completed-stage result paths.
- `runs/orthofinder32_core_20261030/`: upstream results and checkpoints.

The supervisor automatically advances from core inference through full assignment and gene-accounting export. No follow-up command is required while it remains running. It uses six search threads and one analysis thread, retaining at least 6 GiB free disk by stopping with checkpoints if necessary. Idle-sleep inhibition was started for the supervisor lifetime; closing/powering off the computer still interrupts execution.

## Validation performed before the long run

- Real two-genome installation/output test: S288C + RM11-1a, all 11,379 coding units accounted for in 5,120 MCL groups plus 344 unassigned rows. These numbers are **only the two-genome test**, not the 279-genome result.
- Unit/integration suite: 19 tests passed, covering existing sequence-evidence guards, metadata classification, strain identity boundaries, representative selection, retention of the sparse proteome, and full input accounting.
- Original OrthoFinder 3.1.5 was used only during installation testing; production uses 3.2, whose stage/restart behavior and result names are recorded separately. Old smoke directories are not incorporated as production biological results.

## Final outputs generated after successful completion

`outputs/orthology_20261030/orthogroups/`:

- `gene_membership.tsv`: locus/protein/coordinate links for every input gene, with assignment status.
- `copy_number.tsv`: orthogroup × all 279 genomes, counts of coding loci.
- `presence_absence.tsv`: observed annotation membership, not validated biological absence.
- `presence_absence_scope_aware.tsv`: the same view with N/A where relevant nuclear annotation or supplied organelle scope is absent; original counts are preserved.
- `orthogroup_summary.csv`: prevalence, single/multicopy counts and core/accessory summaries.
- `summary.json`: completion and whole-cohort gene-accounting checks.

Orthogroups can contain paralogs and are not equivalent to strict one-to-one ortholog pairs. Pairwise orthologs and trees are retained in the upstream OrthoFinder result. A zero matrix cell means no submitted gene model assigned to that group; subsequent genome searches, neighborhood checks and assembly-boundary inspection are needed before a biological-loss claim.

## Reproduction and resume

```sh
python3 scripts/prepare_orthology_cohort.py --fetch-biosamples
python3 scripts/prepare_all_gene_inputs.py --download
python3 scripts/select_orthofinder_core.py --core-size 64
python3 scripts/run_scalable_orthology.py --threads 6 --analysis-threads 1
```

Do not rerun core selection or regenerate proteomes underneath an active run. Inputs are checksum-verified before launch, but external mutation during a run is not supported.

After an actual interruption, inspect status/logs, stop any surviving worker processes before resuming, then:

```sh
python3 scripts/run_scalable_orthology.py --resume --threads 6 --analysis-threads 1
```

The scalable runner above is the production entry point. No Meta model or biological-loss classifier has yet interpreted these cohort-wide groups.
