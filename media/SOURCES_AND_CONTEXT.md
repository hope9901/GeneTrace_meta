# Sources and biological context for the submission cut

## What the gene functions do and do not establish

- **nimT**: the shipped GCF_001558995.2 GFF describes a 2-nitroimidazole transporter, marked `pseudo=true` and `Note=frameshifted`. Reference E. coli biology links NimT to compound transport. This does not establish that the JEONG-1266 product is functional or that the strain has a measured drug-response phenotype. [UniProt P76242](https://www.uniprot.org/uniprotkb/P76242/entry); [experimental E. coli transporter study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9405208/).
- **nimR**: reference annotation identifies a DNA-binding transcriptional regulator/repressor. Regulation of the nimT locus has been reported in E. coli K-12. The video explains this as control of gene activity, rather than implying regulation was measured in JEONG-1266. [NCBI Gene 946590](https://www.ncbi.nlm.nih.gov/gene/946590).
- **yeaO**: the shipped annotation describes a DUF488 protein; the reference record lists an uncharacterized conserved protein. The film explicitly keeps its precise function unclear. [NCBI Gene 945675](https://www.ncbi.nlm.nih.gov/gene/945675).

These two neighbors were selected as local sequence controls, not because they were proven essential, unaffected, or independent. Their disagreement helps investigate the spatial scope of a signal. It neither excludes a coincident nimT assembly error nor identifies the cause. Biological variation, technical effects and run/version mismatch remain possible.

## BUSCO direction

BUSCO estimates completeness using lineage-specific expectations for conserved single-copy orthologs. These markers are not equivalent to a set of all experimentally essential genes. The proposed extension would prioritize apparently missing or fragmented markers for targeted sequence/annotation/read investigation. It is not implemented in this prototype. A missing BUSCO marker does not independently prove sequencing failure. [Official user guide](https://busco.ezlab.org/busco_userguide).

## Actual demonstration

Recorded report: f64d1b2f2a5c4761b1d023d3bd125796. Source: demo/verified_execution_report.json. Four fresh tools, five validated Meta responses. The video uses enlarged crops of the recorded investigation UI, not a new live run during video playback.

- Population context: saved 100-genome analysis.
- Fresh regional analysis: 756 reference-selected original Illumina paired records; no unbiased whole-run depth or absence claim.
- Coherent read disagreement at nimT exposes a neighboring-control tool in application code. Meta chooses among allowed follow-ups; required steps are not unrestricted model choices.
- Model selection reasons, actual tool completion, neighboring measurements and the provisional judgment are shown in matching order.
- The fixed-policy arm reproduced measurements; the demonstration does not establish an AI accuracy advantage.

## Non-specialist language choices

Reads, assembly and annotation are explained once as pieces, reconstruction and labels. "Mismatch" is explained on screen as read/genome base disagreement. Unexplained orthogroup, synteny, HMM, CDS and frameshift terminology is removed from the spoken script. NimT is described as "labeled as disrupted" rather than asking the audience to know frameshift biology. Reference gene functions are separated from this strain's unmeasured biological activity. BUSCO is introduced as a future marker-based completeness tool.

The MP4 remains local and must not be restored to GitHub without a new user request.
