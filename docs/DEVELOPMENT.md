> Historical seed-control setup. For the current self-contained submission demo, start with README.md and `python3 -m genetrace.demo --port 8766`. The roadmap statements below describe the earlier baseline.

# Running the evidence baseline

Python 3.9+ standard library only. NCBI BLAST+ is a separate executable dependency.

## Setup

1. Download BLAST+ for your operating system from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ . Extract under `tools/` or set `BLAST_BIN` to its `bin` directory. The first local run used BLAST+ 2.17.0 ARM64 macOS.
2. Download the reference inputs and build controls:

```sh
python3 scripts/download_data.py
python3 scripts/prepare_demo.py
```

3. Investigate a control:

```sh
python3 -m genetrace --dataset data/demo/annotation_removed --output runs/annotation-removed-001
```

The output directory must be new. Outputs include copied input files, checksum manifest, tool versions, complete command arguments and stdout/stderr, parsed alignments, evidence JSON, report JSON and Markdown. BLAST commands reference the copied inputs. After moving a bundle, replace the original absolute input paths with the new bundle's `inputs/` paths to rerun.

4. Validate scientific guards and four real-search scenarios:

```sh
python3 -m unittest discover -s tests -v
```

Integration tests skip if data or BLAST executables are missing; skipped checks do not validate the pipeline.

## Current behavior and interpretation

| Control | Expected result |
|---|---|
| Intact | annotated_homolog |
| ACT1 annotation/protein removed, DNA unchanged | annotation_gap_candidate |
| ACT1 annotation/protein removed, expected locus masked with N | assembly_limited |
| ACT1 annotation/protein removed, expected locus replaced with A | inconclusive |

This is a **seeded local-locus analysis**, not a de novo discovery or comparison of different genomes. The seed location is known from the public reference. The nucleotide input is cropped to 15 kb on each side of ACT1; the proteome search uses the full source proteome, with ACT1 removed in perturbed cases. Related actin proteins can remain and produce strong hits. These do not establish ACT1 orthology or family-wide presence/absence.

The control generation retains source records unchanged and writes transformed files separately. Source accessions, URLs, file checksums and crop coordinate offsets are recorded. Dataset snapshots under `data/source` and `data/demo`, local executables and generated runs are excluded from Git. For a public release, distribute a pinned checksum manifest and confirm applicable source data policies before packaging data. NCBI policy: https://www.ncbi.nlm.nih.gov/home/about/policies/ .

Single-HSP query coverage is used for triage; multiple exon hits are not chained. Candidate detection does not establish a full gene model, function, or orthology. The classifier never emits `loss_supported` in this version. GFF attributes are retained in evidence; the demo assumes NCBI GFF3 inputs, not arbitrary user GFF dialects.

## Next implementation steps

- Add Meta tool calling with bounded steps and persisted call/response trace, keeping baseline and model-driven runs visibly distinct.
- Add a second independently assembled genome, seed localization with flanking homologs, and spliced alignment where needed.
- Build the review interface and export a relocatable evidence bundle.
- Add pinned public data acquisition and three-minute video workflow.

No Meta API calls, hosting, public repository creation or publication are performed by the baseline.
