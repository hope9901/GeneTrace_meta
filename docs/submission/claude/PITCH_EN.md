> Historical draft, superseded by the verified execution workflow and final video. See ../FRESH_EXECUTION.md and ../EXECUTION_COMPARISON.md. Old PENDING markers and proposed UI controls are not the current product state.

# GeneTrace — Pitch Copy (EN)

Source of truth: `docs/submission/FACTS.md` (2026-11-03). Every number below is taken from it or from the reports it cites. Anything marked **[PENDING LIVE META VERIFICATION]** must stay that way until a real Meta model run with valid evidence citations has been recorded.

---

## One-liner

**GeneTrace is an evidence-first research agent that checks whether a specific stretch of a public genome can be trusted before a biologist builds on it.**

---

## Short description (50 words)

> GeneTrace helps biologists check whether a suspicious gene region in a public genome is trustworthy before follow-up work. Instead of judging from genome-wide statistics, it gathers local evidence (reads, neighboring genes, external contamination reports) and shows evidence IDs, exact commands, versions, file hashes, remaining uncertainty and next checks to run.

*(50 words)*

---

## Long description (250 words)

> Biologists often pick a gene from a public genome assembly and start an experiment, a comparison or a paper. Genome-wide quality scores may look excellent while one gene is broken by an assembly artifact, a sample mix-up or contamination. Checking that properly means reads, alignments and provenance metadata: work that sits between biologists and computational analysts.
>
> GeneTrace bridges that gap. A user asks about one region; the agent assembles local evidence and reports it with evidence IDs, the exact commands, tool versions and file hashes behind each result, what is still uncertain, and what to check next. It does not call an error from genome-wide statistics alone, and it never treats a minority allele or indel as an automatic error.
>
> We built it on public NCBI and ENA data: 279 deduplicated yeast strains and 100 quality-filtered E. coli genomes selected from 5,512 candidates. Three real cases show the approach. In nimT, reads disagree with the assembly, but neighboring genes disagree too, so the cause stays open. In another sample, unusually regular reads trigger a provenance warning instead of being counted as support. Against NCBI FCS contamination calls in three flagged assemblies, every eligible actioned contig was flagged, alongside 94 extra flags we cannot yet explain.
>
> A Meta model plans the investigation: it may only choose allow-listed evidence tools, and its claims must cite evidence IDs that code verifies. The underlying checks are deterministic and run without the model. [PENDING LIVE META VERIFICATION]
>
> These are curated cases, not an accuracy benchmark.

*(250 words — recount after any edit)*

---

## Section building blocks

Use these when the submission form has separate fields.

### Problem / users
- **Users:** bench and computational biologists who reuse public genome assemblies (here: *S. cerevisiae* and *E. coli*) as the starting point for follow-up research.
- **Problem:** genome-level QC (completeness, contamination scores) does not answer "is *this* gene region reliable?" A single local artifact can send a follow-up study in the wrong direction. Answering the local question needs reads, alignment, neighboring-gene context and provenance — skills and time many biologists don't have.

### Product
- A user names a gene/region; GeneTrace gathers local evidence and presents it as an evidence trail.
- Every finding carries an **evidence ID**, the **actual command, tool version and file hash**, the **remaining uncertainty**, and **suggested next checks**.
- Neighbor/control comparison is built in: a disagreement is only interesting if the surrounding region behaves differently.
- Raw evidence (with source SHA256) can be exported as JSON.
- The demo replays **stored, real analyses** of public data; it does not rerun alignments live.

### Meta's role
- A Meta model acts as the investigation planner and explainer: it selects from an **allow-list of evidence tools** and must cite **evidence IDs**, which code checks against the stored evidence before display.
- Core evidence checks are deterministic and do **not** use the model; the model cannot invent a measurement.
- Status: code path implemented; live call **[PENDING LIVE META VERIFICATION]**. Do not name a model ID or claim a successful run until verified.

### Technical differentiators
- **Local, not genome-wide judgments.** Read depth, mismatch rate and indel disagreement per gene, always compared with flanking genes.
- **No automatic error calls.** Minority alleles and indels are reported as evidence, not verdicts.
- **Provenance awareness.** Read sets with implausibly regular structure are flagged for provenance review and excluded from "independent sequencing support".
- **Pre-registered external comparison.** Contamination detector thresholds were fixed before comparing with NCBI FCS; FCS coordinates were never used as detector input.
- **Reproducibility.** Genome/GFF/CDS/protein/assembly-report files (500 for the E. coli panel) checked against publisher MD5; where no publisher MD5 exists, a SHA256 snapshot is recorded.

### Validation and limits
- **Data:** 279 yeast strains after operational strain/BioSample deduplication (not sequence-verified clone deduplication). 100 E. coli genomes: 2,161 of 5,512 Complete RefSeq candidates passed CheckM filters (completeness ≥99%, contamination ≤1%, etc.); MG1655 plus 99 chosen by an explicit quality ranking and deduplication. Not a claim of species-wide representativeness.
- **Real-read cases:** 63 of the 100 E. coli genomes have WGS runs (220) under the exact BioSample; whether a given run was used for the specific assembly version is unconfirmed.
- **FCS comparison:** three high-contamination positive cases, contig-level concordance only. Only 131 of 3,095 FCS-actioned contigs were ≥2 kb and therefore testable; all 131 were flagged. 94 contigs without FCS action were also flagged — specificity is not established, and "no FCS action" is not a confirmed negative.
- **Gene-family context:** 9,233 fast protein homology groups (MMseqs2, 90% identity / 80% bidirectional coverage) — candidate groups, not orthogroups. Formal OrthoFinder runs for yeast and E. coli are still in progress.
- **Tests:** existing automated tests are development checks, not real-world diagnostic accuracy.
- **Not claimed:** accuracy, detection of all errors, confirmed gene loss, confirmed assembly errors, or self-validation via FCS.

### Attribution line (fill before submitting)
Built during the hackathon: investigation pipeline, data audit, evidence-linked agent and UI. Third-party: public NCBI/ENA data; NCBI FCS reports; CheckM metadata; minimap2; MMseqs2; OrthoFinder; Pfam; [add remaining tools/versions from README once Codex finalizes it].
