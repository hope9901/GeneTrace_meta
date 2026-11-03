> Historical draft, superseded by the verified execution workflow and final video. See ../FRESH_EXECUTION.md and ../EXECUTION_COMPARISON.md. Old PENDING markers and proposed UI controls are not the current product state.

# GeneTrace — 180-second demo video plan

- Demo: `python3 -m genetrace.demo --port 8766` → http://127.0.0.1:8766/ (per `SUBMISSION_CHECKLIST.md`).
- Screen/click column describes the *intended* shot. Claude has not seen the UI (`web/`, `genetrace/` were out of scope), so button/panel names in `[brackets]` must be matched to the real UI before recording.
- Narration total: **349 words** (45/47/96/61/48/52 per segment, ≈116 wpm; the 0:50–1:35 segment is densest at ≈128 wpm). Recount after edits; target 320–380.
- Rule for every shot: stored analyses are labelled as stored; only the Meta call is live. Never show `.env` or a key.

| Time | Screen | Clicks / actions | Narration (EN) |
|---|---|---|---|
| 0:00–0:25 | Title card "GeneTrace", then GeneTrace home with the E. coli 100-genome panel visible. | None (static). Optional slow zoom on panel summary. | Before a biologist builds an experiment on a gene, they need to know whether that stretch of genome can be trusted. Genome-wide quality scores can look excellent while one gene is broken. GeneTrace investigates that local question with evidence, before anyone invests in follow-up work. |
| 0:25–0:50 | Meta configuration status (model ID shown, key hidden), then the question box. Investigation trace appears: tool selections + evidence IDs. | Show `[Meta status]` badge → type question about nimT in GCF_001558995.2 → `[Investigate]` (live Meta mode). Keep real latency or mark cut with on-screen "time skipped". | I ask about nimT in E. coli assembly GCF_001558995.2, which is annotated as frameshifted. A Meta model chooses from an allow-list of evidence tools, and every claim must cite an evidence ID that code verifies. The read results it retrieves were computed earlier; nothing is realigned live. **[PENDING LIVE META VERIFICATION]** |
| 0:50–1:35 | nimT evidence card: run SRR2762258, depth, mismatch, indel disagreement, strand split. Then neighbor comparison view: nimR and yeaO side by side. Uncertainty + next-checks panel. | Click evidence ID for nimT read check → open `[neighbors/controls]` → scroll to `[uncertainty / next checks]`. Briefly hover command/version/hash line. | The stored evidence uses public run SRR2762258. Across nimT, mean depth is about 31, and about 3.1 percent of aligned bases disagree with the assembly. At one single-base indel, 21 of 23 alignment records disagree, on both strands; those are records, not independent molecules. That looks like a local error. But GeneTrace checks the neighbors: nimR and yeaO also disagree, at about 2 percent. So we do not call an assembly error. The run may not match this assembly version, or the sample may be mixed. The cause stays unresolved, and the next checks are listed. |
| 1:35–2:05 | Provenance case: GCF_013112335.1 / SRR11451201. Four gene rows (ybiU, tonB, ybiV, ybiT) with depth 25, mismatch 0; read-start spacing, template length, quality summary; `provenance review required = true` flag. | Open case → expand `[read structure]` → point to the provenance flag. | Perfect support can also be a warning. For GCF_013112335.1, run SRR11451201 shows exactly 25-fold depth and zero mismatches across four genes, reads starting every 10 bases on the same strand, a fixed 450-base template length, and nearly one quality character. ENA lists it as Illumina. GeneTrace does not call it synthetic, but it will not count it as independent sequencing support. |
| 2:05–2:30 | FCS comparison table (three GCA assemblies) with **both** columns: eligible ≥2 kb denominator and "extra flags without FCS action". | Open `[contamination validation]` → highlight the 2,922 → 1 eligible row → highlight the 31 and 63 extra flags. | For external checking, we took three E. coli assemblies NCBI flags as contaminated, fixed our thresholds first, then compared with NCBI FCS calls. Only 131 FCS-actioned contigs were long enough to test; all were flagged. We also flagged 94 contigs FCS left alone, so specificity is not established. |
| 2:30–3:00 | Export: raw evidence JSON with source SHA256 visible. Then end card: value, Meta role, limits. | Click `[Export evidence JSON]` → show file with hashes → cut to end card. | GeneTrace gives biologists an evidence trail: evidence IDs, exact commands, tool versions, file hashes, uncertainty, and next checks. Meta's model plans and explains, but only from cited evidence; the checks themselves are deterministic. These are a few curated public-data cases, not an accuracy benchmark, and our formal ortholog analysis is still running. |

## End card text (on screen, not narrated)

- Built at the hackathon: investigation pipeline, data audit, evidence-linked agent, UI.
- Uses public NCBI/ENA data, NCBI FCS reports, and third-party tools (list in README).
- Meta model: [MODEL ID] — **[PENDING LIVE META VERIFICATION]**
- Repo: [URL] · Not an accuracy benchmark.

## Recording checklist (shots that must not be missed)

1. Meta status visible with model ID and **no key**; a real live call, not fallback/no-model mode. If the call fails, fix and re-record — do not edit around it.
2. Tool-selection trace and at least one clicked evidence ID resolving to stored evidence.
3. A "stored result" label (or narration) on read/alignment evidence.
4. nimT **and** nimR/yeaO on screen together.
5. Provenance flag on SRR11451201.
6. FCS table showing eligible denominator and extra-flag column, not just the "flagged" column.
7. One command/version/hash line and the JSON export.
8. If time is skipped while waiting for the model, show an on-screen "time skipped" marker.

## If live Meta verification is not done by recording time

Replace the 0:25–0:50 narration's last two sentences with: "In this recording the model step is shown from code, not a live call; the read results were computed earlier." and remove the model ID from the end card. Do not imply a successful call.
