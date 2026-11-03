# Fresh nimT investigation

## Inputs and scope

- Assembly: GCF_001558995.2, E. coli JEONG-1266; whole reference genome and structural GFF from NCBI RefSeq.
- Read run: SRR2762258, ENA Illumina MiSeq paired-end WGS; BioSample SAMN04191546. Assembly metadata records PacBio/FALCON. Exact run-to-assembly-version use is unverified.
- Execution subset: 756 pairs copied unchanged from publisher-MD5-verified original FASTQ records. Read names were selected by an earlier whole-run alignment window around nimT and its neighbors. The subset is reference-dependent; unmapped/divergent evidence may be missed. It is not independent validation, unbiased whole-run depth, an absence test, or whole-run QC.
- `demo/execution/nimt/manifest.json` pins source MD5, selected counts, extraction window and distributed input SHA256.
- `scripts/prepare_execution_case.py` reproduces extraction when the full local raw inputs and BAM are available. Those full inputs are not shipped.

## State-dependent execution

1. Saved 100-genome population context is loaded and labeled as saved.
2. Fresh structural annotation parsing is the first required tool.
3. Fresh minimap2 paired short-read alignment to the whole assembly, SAMtools sorting/indexing and candidate regional analysis follow.
4. Candidate coherent disagreement exposes a new neighboring-control analysis tool. A read-structure/source-linkage audit is also required.
5. Meta chooses among currently available checks and sees each returned result before the next choice. The shared minimum-check policy prevents early termination without required evidence. Model commands, paths and arguments are never executed.
6. The deterministic structured judgment is updated after each tool. Final Meta interpretations must cite existing evidence IDs. No confirmed causal diagnosis or biological loss is assigned.

Only this fixed case is enabled for fresh execution. Other cases remain saved evidence review. This is a bounded working prototype, not an arbitrary uploaded-genome analysis service.

## Comparison

The fixed baseline uses the same input pack, tools and conditional minimum-check policy, with a fixed order and no model. The Meta arm adds tool-selection rationale and cited synthesis. Both can find the neighbor discrepancy. The one-case comparison does not establish improved accuracy, expert time savings or unique AI discovery. See EXECUTION_COMPARISON.md.

## Reproducibility

Each tool writes a JSON artifact with input hashes, analysis source hash, commands, versions when called, return codes and timings. Evidence cites the tool artifact SHA256. Generated alignments are stored outside source under `runs/execution/<job-id>`.

`python3 -m unittest discover -s tests -p test_execution.py -v` checks input tampering rejection, fixed tool allowlists, conditional follow-up, incomplete judgments and unresolved cause guards without requiring external binaries.
