# GeneTrace video narration

## GeneTrace

A gene appears to be missing. Is that biology, or a data problem? GeneTrace investigates genome reliability before researchers build experiments on an uncertain result.

## A missing label is not a missing gene.

Sequencing produces DNA fragments called reads. Assembly reconstructs a genome, and annotation marks gene locations. Think of reconstructing a book and marking its chapters. A missing chapter label does not prove the text is gone.

## How Meta powers the agent.

Meta's Muse Spark one point three receives the research question, current evidence, and allowed checks. It returns a tool choice, reason, and cited hypotheses. GeneTrace validates the response, executes the check, and sends new evidence back. A final Meta call synthesizes the findings.

## Public genomes. Selected read checks.

We use one hundred quality selected, complete Escherichia coli genomes from NCBI RefSeq, with gene annotations. MG1655 anchors the comparison panel. For this demonstration, seven hundred fifty six read pairs were selected from a public Illumina run. This reference selected subset is not whole run quality control or independent validation.

## Meta selects a check. Code executes it.

We select nimT, a gene with frameshift related annotation, and enable fresh analysis. Meta requests structural annotation inspection, then read alignment. The application actually executes minimap2 and regional checks on fixed public inputs, and returns newly computed evidence.

## A new result opens a follow-up.

The candidate disagrees with selected reads. That observation exposes a neighboring control tool. Meta requests it to investigate whether the discrepancy is local. The application calculates nimR and yeaO from the new alignment, then updates the judgment. These are executed checks, not just saved result lookups.

## A discrepancy needs context.

The neighboring genes also disagree. A strictly local assembly error is therefore insufficient as an explanation. The structured result lists possible biological, assembly, read, or dataset linkage causes, with limitations and next verification. Exact run to assembly version linkage remains unverified; the cause is unresolved.

## What does Meta add?

We ran the same case with a fixed policy and with Meta. Both executed four fresh checks and produced equal measurements. Four Meta planning calls and one synthesis added selection reasons and cited interpretation. Both share minimum check safeguards, so this single case does not establish improved accuracy or unique AI discovery.

## Bounded reasoning. Traceable execution.

Population context is saved; regional checks are newly executed on a selected subset. Commands, versions, hashes and cited evidence support review. This is a bounded research prototype.

## Evidence before interpretation.

GeneTrace connects biological context, real analysis execution, and Meta guided interpretation, helping researchers investigate uncertain genome regions before deciding what to trust and what to verify next.