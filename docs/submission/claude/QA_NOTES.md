> Historical draft, superseded by the verified execution workflow and final video. See ../FRESH_EXECUTION.md and ../EXECUTION_COMPARISON.md. Old PENDING markers and proposed UI controls are not the current product state.

# QA notes — submission copy and video

Checked against `docs/submission/FACTS.md`, `RELEASE_PLAN.md`, `SUBMISSION_CHECKLIST.md`, `docs/CONTAMINATION_VALIDATION.md`, `outputs/ecoli_followup_20261102/report.md`, `outputs/ecoli100_candidate_audit_20261102/report.md`. Code, UI and raw outputs were not read (out of scope), so UI-dependent items are unverified.

## P0 — blocks submission

1. **Meta live call not verified.** FACTS: `.env` setup and real API call unverified; `submission.json` has `live_Meta_verification: "pending"`. PITCH_EN.md and VIDEO_180S.md keep `[PENDING LIVE META VERIFICATION]`. Before submitting: record one successful run with model ID, call count, tool trace and valid evidence-ID citations, then replace the tag. Until then, no sentence may say the model "ran", "investigated" or "found" anything. RELEASE_PLAN treats this as a mandatory requirement.
2. **Model ID unknown.** Copy says "a Meta model". Insert the exact model ID only after verification; never guess a product name.
3. **Repository / video / demo URLs are empty** in `submission.json`. End card and form fields have placeholders.
4. **Attribution list incomplete.** PITCH_EN.md lists tools visible in the reports (minimap2, MMseqs2, OrthoFinder, Pfam, CheckM metadata, NCBI FCS). Versions, licenses and any other tools must come from Codex's final README. "Built during the hackathon" vs pre-existing work must be confirmed by the owner.

## P1 — numbers that need confirmation or careful wording

5. **nimT 21/23 indel figure is only in FACTS.** The follow-up report lists `coherent_indel_disagreement` for nimT but gives no 21/23 or 15 forward / 6 reverse counts. Confirm against the raw analysis manifest before narrating it.
6. **Direction of the nimT indel disagreement is ambiguous.** FACTS says "1bp deletion disagreement". It's unclear whether the assembly carries the deletion and reads don't, or the reverse. The script says "at one single-base indel, 21 of 23 alignment records disagree" to stay neutral. Confirm and, if possible, state the direction. This matters because nimT is annotated frameshifted (candidate audit).
7. **Provenance case details are only in FACTS.** The 10 bp start step, fixed template length 450 and near-single quality character are not quantified in the follow-up report (which only says "regular start spacing, fixed template length"). Strain name "GF3-3" is also only in FACTS. Confirm against the manifest before recording; the script doesn't use the strain name.
8. **FCS denominators: always show both sides.** Correct framing: 131 of 3,095 FCS-actioned contigs were ≥2 kb and testable, and all 131 were flagged. 94 extra flags (0 + 31 + 63) among 556 non-actioned contigs (13 + 458 + 85, from the follow-up report; FACTS omits this column). Do **not** say "matched all FCS contigs": for GCA_027889185.1 only 1 of 2,922 actioned contigs was testable.
9. **Mismatch rounding.** Report: nimT 0.030853 (3.09%; script says "about 3.1 percent"), nimR 0.020657, yeaO 0.019601 ("about 2 percent"). Depth nimT 31.370. Consistent.
10. **Accession prefixes.** Read cases use GCF_ (RefSeq); FCS cases use GCA_ (GenBank). Don't mix them on screen or in captions.
11. **nimT appears in two genomes.** The candidate audit also lists nimT for GCF_031588155.1. The read case is GCF_001558995.2 only; don't generalize.
12. **"5,512 candidates → 100" wording.** 2,161 passed filters; MG1655 + 99 chosen by explicit ranking and deduplication. Avoid "best 100 E. coli genomes" or "representative of the species".
13. **Yeast 279.** Operational strain/BioSample deduplication, not sequence-verified clone deduplication. The pitch says "deduplicated yeast strains"; the long form explains it in the limits section. Shorten further if needed but keep the qualifier somewhere.
14. **Run-to-assembly link unconfirmed.** 63/100 genomes have 220 WGS runs under the exact BioSample, but whether a given run produced the specific assembly version is unconfirmed. The nimT narration uses this as a possible explanation; don't present the run as "the assembly's reads".

## P1 — overclaim risks (wording to avoid)

15. "Every eligible actioned contig was flagged" must always appear with the 94 extra flags and "specificity not established". Never say FCS "validates" GeneTrace's accuracy.
16. Never call the 9,233 MMseqs2 groups or the 3,027 single-copy groups "orthogroups". Formal OrthoFinder is still running (yeast core inference; E. coli additional 84 assignment). The 915 "missing combinations" are not gene losses; tonB is an explicit counterexample (homologous DNA and annotation present).
17. The provisional tree (30 proteins, 7,664 columns) has no recombination correction; the copy doesn't mention it. Keep it out of the video.
18. "80 tests pass" must not appear as an accuracy figure. Record the final test/skip count separately just before submission.
19. HN1 PacBio (ERR4194695) is in FACTS but not in the script. If shown, don't infer low quality from `!` quality strings, and note that ZMW independence and run-version use are unconfirmed.
20. Pfam search covered 6 profiles, not all of Pfam.

## P2 — Meta contribution clarity

21. Copy describes Meta's role as: choose allow-listed tools, explain, cite evidence IDs that code validates. Core checks run without the model. Make sure the UI actually shows (a) which tools the model chose and (b) citation validation. Otherwise judges can't see a Meta contribution beyond text generation.
22. If the demo has a no-model fallback mode, the recording must show it is off. A fallback run must not be presented as a Meta run.
23. Hackathon judges may ask "what fails without Meta?". Prepare one sentence: without the model, users browse the evidence manually; the model turns a natural-language question into a tool plan and a cited explanation.

## P2 — recording gaps

24. UI element names in VIDEO_180S.md are placeholders (`[Investigate]`, `[neighbors/controls]`, `[Export evidence JSON]`, etc.). Map them to the real UI.
25. Visibly label stored vs live content (on-screen badge or caption), not only in narration.
26. If model latency is cut, add an on-screen "time skipped" marker (checklist: don't misrepresent latency).
27. Confirm the FCS table in the UI shows the eligible ≥2 kb column and the extra-flag column. If it shows only "flagged", fix the UI (Codex) or add an overlay.
28. Hide `.env`, terminal history and environment variables. Check browser autofill and tab titles.
29. Optional port 8765 OrthoFinder status shot: label it "in progress". Don't show partial counts as results.
30. Deadline: internal target is 2026-11-04 morning (Asia/Seoul). The official 11/4 21:00 EST deadline is unverified; the owner must check the platform's current requirements and video length limit.
