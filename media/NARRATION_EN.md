# GeneTrace video narration

## GeneTrace

A gene appears to be missing. Is that biology, or a data problem? GeneTrace investigates genome reliability before researchers build experiments on an uncertain result.

## A missing label is not a missing gene.

Sequencing produces DNA fragments called reads. Assembly reconstructs a genome, and annotation marks gene locations. Think of reconstructing a book and marking its chapters. A missing chapter label does not prove the text is gone.

## Automate the investigation.

Experts need data collection, computational analysis, and biological interpretation together. Our goal is to automate that investigation. A rare difference can be natural, so the agent must also know when the cause is unresolved.

## Public genomes. Selected read checks.

We use one hundred quality selected, complete Escherichia coli genomes from NCBI RefSeq, with gene annotations. MG1655 anchors the comparison panel. For this demonstration, seven hundred fifty six read pairs were selected from a public Illumina run. This reference selected subset is not whole run quality control or independent validation.

## Meta selects a check. Code executes it.

We select nimT, a gene with frameshift related annotation, and enable fresh analysis. Meta requests structural annotation inspection, then read alignment. The application actually executes minimap2 and regional checks on fixed public inputs, and returns newly computed evidence.

## A new result opens a follow-up.

The candidate disagrees with selected reads. That observation exposes a neighboring control tool. Meta requests it to investigate whether the discrepancy is local. The application calculates nimR and yeaO from the new alignment, then updates the judgment. These are executed checks, not just saved result lookups.

## A discrepancy needs context.

The neighboring genes also disagree. A strictly local assembly error is therefore insufficient as an explanation. The structured result lists possible biological, assembly, read, or dataset linkage causes, with limitations and next verification. Exact run to assembly version linkage remains unverified; the cause is unresolved.

## What does Meta add?

We ran the same case with a fixed policy and with Meta. Both executed four fresh checks and produced equal measurements. Meta added evidence dependent selection reasons and cited synthesis. Both share minimum check safeguards, so this single case does not establish improved accuracy or unique AI discovery.

## Show what the evidence can support.

Population context is still a saved analysis. Fresh execution is enabled for one bounded case. Commands, versions, hashes, and evidence exports make the investigation reviewable. This is a working prototype, not a validated diagnostic system.

## Evidence before interpretation.

GeneTrace connects biological context, real analysis execution, and Meta guided interpretation, helping researchers investigate uncertain genome regions before deciding what to trust and what to verify next.