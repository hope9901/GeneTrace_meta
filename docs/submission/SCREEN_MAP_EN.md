> Update 2026-11-03: fresh nimT execution is implemented and verified. See [FRESH_EXECUTION.md](FRESH_EXECUTION.md) and [EXECUTION_COMPARISON.md](EXECUTION_COMPARISON.md). DIAMOND/OrthoFinder jobs are paused. The final video is generated from actual execution captures plus explanatory cards; older snapshot-only notes below describe the earlier stage.

# English recording screen map

The submission site is http://127.0.0.1:8766/ . The research workspace at port 8765 remains a separate interface.

- Case selection: **Choose a case** → **A disrupted gene, or a mismatched dataset?** (nimT).
- Live model path: **Review mode → AI-guided review → Start investigation**. This makes new model calls over stored evidence; no alignments run live.
- Existing report: **Saved investigation** means a previously completed review, not a newly executed call. **Run details** holds the model ID, validated response count and timestamp.
- Evidence: **Population & neighborhood**, **Sequence & structural annotation**, **Read provenance & integrity**, **Candidate & local controls**.
- Measurements: **Mean depth**, **Mismatch**, **Coherent indels**, **Local controls**.
- Audit: **Inspect source evidence** opens measurements and source hashes. **Export evidence JSON** saves the report.
- Interpretation: **Interpretation & next checks**; evidence-ID links identify cited items.
- Other cases: **Perfect matches. Uncertain provenance.** and **Testing against known contamination**.

Removed from the primary flow: API setup/configuration messages, local .env instructions, raw function names, repeated hash boilerplate, and technical call counters. Snapshot/live distinction and scientific limitations remain visible. Original saved model responses are preserved; new demo responses are generated in English.
