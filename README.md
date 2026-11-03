# GeneTrace

**Investigate genome reliability before biological interpretation.**

GeneTrace is a research assistant that connects same-species population context, structural annotation, read–assembly disagreement and data provenance. It helps researchers decide what to verify before spending time on an apparent gene loss, expansion or disrupted model. It preserves uncertainty rather than treating every indel or minority allele as an error.

Built for the Meta Global AI Developer Hackathon, Indie Dev Track. The investigation workflow, evidence guards, analysis scripts and interfaces were built during the event. NCBI/ENA data and established third-party tools are credited separately below.

## Demo video

[Watch or download the 2m54s English video](media/GeneTrace_EN.mp4) · [English subtitles](media/GeneTrace_EN.srt) · [Narration](media/NARRATION_EN.md). Actual execution captures are time-condensed; explanatory cards are interleaved.

## Run the small submission demo

Python 3.9+ standard library for the web app and saved evidence review. Fresh execution additionally requires minimap2 and samtools. No full cohort download or completed orthogroup job is needed.

```sh
python3 -m genetrace.demo --port 8766
```

Open **http://127.0.0.1:8766/**. Choose one of three cases and start an investigation. Inspect the evidence IDs, measurements, source SHA256 and exported JSON.

**Two explicit paths:** saved evidence review for all three cases, and **fresh local execution for nimT**. Enable **Run fresh nimT analysis** to parse GFF, realign a public read subset with minimap2, analyze the candidate region, and expose neighboring-control checks when coherent disagreement is observed. Meta selects among the currently allowed tools and reassesses returned evidence. A fixed minimum-check policy bounds both modes. The guarded structured judgment is deterministic; Meta adds selection rationales and cited interpretations. A configured key enables AI-guided review; otherwise use model-free Evidence explorer.

| Case | What it demonstrates | What it does not establish |
|---|---|---|
| JEONG-1266 / nimT | Read–assembly disagreement plus neighboring controls | A confirmed local assembly error |
| GF3-3 / read provenance | Perfect matching reads can still require source review | Proof that the reads are simulated |
| Three NCBI contaminated E. coli controls | Comparison with external FCS action contigs | General contamination diagnostic accuracy |

The exported snapshots are in `demo/evidence/`; their SHA256 values are pinned in `demo/manifest.json`. Source paths identify the original local analysis artifacts and are accompanied by their snapshot hashes. These large artifacts are not required to read the small demo.

## Run a fresh nimT investigation

Install minimap2 and samtools through your package manager (for example `brew install minimap2 samtools`, or `conda install -c bioconda minimap2 samtools`) and keep both on `PATH`. The web app also recognizes the original local `tools/readqc-env/bin` installation. Start the server, select the nimT case, enable **Run fresh nimT analysis**, and choose AI-guided review or Evidence explorer.

The included execution case is **756 paired FASTQ records selected from public ENA run SRR2762258**, plus the original NCBI assembly and structural annotation for JEONG-1266 / GCF_001558995.2. Reads were selected by names overlapping a previously aligned window; original sequence and quality records were retained. New execution maps those selected reads to the **whole assembly**. This reference-dependent subset cannot establish unbiased whole-run depth, absence, independent sequencing validation, or general diagnostic accuracy. Source whole-file MD5, extraction method, selected counts and SHA256 are pinned in `demo/execution/nimt/manifest.json`. [Case protocol](docs/submission/FRESH_EXECUTION.md).

`runs/execution/<job-id>/` contains new alignment and analysis artifacts. Exports identify saved population context separately from fresh execution, with commands, tool versions, timings, input hashes and output hashes. Missing tools, altered inputs and API failures fail visibly; there is no silent model-free fallback.

[Paired fixed-policy vs Meta case study](docs/submission/EXECUTION_COMPARISON.md): both modes executed four fresh tools and produced equal deterministic observations. One Meta run used five validated responses. This is a single-case demonstration, not evidence of accuracy or productivity improvement.

## Live Meta integration

Create `.env` locally using `.env.example`:

```dotenv
MODEL_API_KEY=your_local_key
META_MODEL=exact_model_id_from_your_Meta_dashboard
META_API_BASE=https://api.meta.ai/v1
```

Never commit the key, show it in a recording, or send it to chat/Claude. Use the exact model available to your account. Refresh the demo after configuring it, select **AI-guided review**, then run a case. The report records the model ID, successful model-call count, selected tools, usage and evidence citations. API or citation-validation failures are displayed as failures; they are not silently converted into model-free success.

**Live verification completed on 2026-11-03:** `muse-spark-1.3` successfully performed four API calls on the nimT case, inspected four evidence items and produced four evidence-cited provisional interpretations. Provider-returned model IDs and citation membership were checked. See `demo/live_meta_verification.json` and `demo/verified_meta_report.json`. This is a saved verification record; opening it does not make a new API call or rerun genomic analysis. A configured key or mocked test alone is not live integration proof.

## Full research workspace

```sh
python3 -m genetrace.web --port 8765
```

This workspace uses the local research dataset and external binaries. It supports S. cerevisiae 279 assemblies and quality-selected E. coli 100 assemblies, with species-specific reference genes, actual BLAST searches, annotation/sequence checks, selected Pfam evidence, analysis status and JSON exports. The background panel checks process identity, stage, heartbeat and disk space every 15 seconds. A core search file percentage is not the full pipeline completion percentage.

The workspace expects data manifests and source files acquired by the research scripts. Use the small submission demo above for a fresh checkout. See [research setup](docs/DEVELOPMENT.md), [agent architecture](docs/AGENT_WORKSPACE.md), [quality-selected bacterial panel](docs/ECOLI_QUALITY100.md), and [restart notes](docs/RESUME_AFTER_SHUTDOWN.md).

## Completed evidence and limitations

- **Yeast:** annotated NCBI assemblies 435 → 279 operational strain/sample representatives. This is metadata deduplication, not sequence-verified clone deduplication. YJU3 was compared across all 279; ambiguous or uncallable loci are excluded from allele-frequency denominators.
- **Bacteria:** all 5,512 Complete RefSeq / non-atypical candidates were surveyed. After explicit quality eligibility and strain/sample deduplication, MG1655 plus 99 ranked assemblies were selected. Their 500 source files were checked against NCBI MD5 manifests. This quality panel is not representative of all E. coli diversity.
- **rpoB:** all 100 loci had single-locus and two-flank support; 50 variable reference repeat blocks were found. These blocks are not normalized VCF calls.
- **Domains:** six selected Pfam profiles were scanned against all 487,647 bacterial protein records. Protein counts and domain-instance counts are separate. This is not full Pfam annotation or orthogroup inference.
- **Candidate context:** fast MMseqs homolog groups and a provisional 30-protein tree were used for candidate triage. They do not establish reconciled orthogroups or recombination-corrected phylogeny.
- **Reads:** exact BioSample WGS links were found for 63 assemblies (220 runs). Two paired-run datasets were checksum verified and mapped. Exact run-to-assembly-version usage remains unresolved.
- **nimT:** a coherent 1bp deletion disagreement was observed in 21/23 alignment records, with both strands represented. Nearby genes also disagree, so sample/version/mixture possibilities remain. Paired records are not independent molecules.
- **Read provenance:** one dataset had regular 10bp starts, fixed 450bp template lengths and near-constant quality with exact local matches. It is flagged for provenance review, not independently confirmed physical sequencing and not proven simulation.
- **FCS controls:** 4/4, 1/1 and 126/126 eligible (≥2kb) external FCS action contigs were flagged in three high-contamination controls. Many shorter contigs were out of scope, and 0, 31 and 63 non-action contigs were additionally flagged. This is contig-level concordance, not sensitivity/specificity or general error diagnosis. [Protocol and limitations](docs/CONTAMINATION_VALIDATION.md).
- **Synthetic validation:** expected regional patterns were recovered on held-out simulations. Simulation recovery is not established real-world diagnostic performance.

Full OrthoFinder analyses for yeast 279 and bacteria 100 were paused at the user’s request; checkpoints and completed files are retained locally. Bacterial coding-locus inputs contain 487,637 entries after choosing the longest protein model per gene ID. Do not describe these long-running analyses as complete until cohort and gene accounting pass.

## Verification

For a fresh checkout, run the self-contained demo guards:

```sh
python3 -m unittest discover -s tests -p test_demo.py -v
```

The full local research test suite requires research datasets and optional analysis dependencies. Skipped tests do not validate the missing components. See `docs/submission/FACTS.md` for the submission fact sheet and current pending items.

## Reproducibility and attribution

Source MD5/SHA256, input paths, versions, commands and stderr are recorded by the analysis scripts. Some auxiliary NCBI FCS reports have no publisher MD5; that absence and the downloaded snapshot SHA256 are recorded explicitly. Saved case summaries are in `demo/evidence`. The separate `demo/execution/nimt` pack includes the public reference genome, structural annotation and a selected paired FASTQ subset for real reanalysis; it does not include all cohort genomes or full raw runs.

Public data: [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/), [ENA](https://www.ebi.ac.uk/ena/browser/), [NCBI FCS reports](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/quality/contamination/fcs-contamination/), [Pfam/InterPro](https://www.ebi.ac.uk/interpro/).

External analysis tools: NCBI BLAST+, OrthoFinder, DIAMOND, MMseqs2, HMMER, MAFFT, FastTree, minimap2 and SAMtools. GeneTrace integrates their outputs; it does not claim to have invented these algorithms. Meta is used for bounded evidence planning and interpretation when live mode is actually executed. [Meta API reference](https://dev.meta.ai/docs/api-reference).

Repository: [hope9901/GeneTrace_meta](https://github.com/hope9901/GeneTrace_meta). Live Meta evidence review and fresh execution were verified locally. Video upload and hackathon submission remain separate steps. No public server deployment is implied. Project code has no open-source license grant yet; copyright remains with the authors. Public data and third-party tools retain their own terms. See NOTICE.md.
