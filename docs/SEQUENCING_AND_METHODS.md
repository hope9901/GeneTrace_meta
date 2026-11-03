# Sequencing, assembly and structural annotation metadata

Snapshot: 2026-10-30. This extends the existing 435-assembly S. cerevisiae cohort. Previously calculated size, GC, proteins and assembly statistics are preserved. The sequencing boxplots use the original 434 Full-representation assemblies, excluding the one Partial assembly.

## NCBI sources and missing values

Each accession was inspected using:

1. Its complete `assembly_report.txt`: sequencing technology, assembly method and coverage.
2. Its complete GFF3: gene, transcript/mRNA, CDS and exon source-column labels, header metadata and explicit ab-initio-prediction inference attributes. GFF compressed-file MD5 was verified against the prior NCBI file manifest.
3. The first GenBank-record header before `FEATURES`, including structured Genome Assembly/Annotation Data and explicit structural-tool comments. This is a bounded header inspection, **not a complete audit of every GenBank record**. These partial responses are not described as full-file MD5 verified.

Where the report is missing a value, explicit GenBank structured metadata is used as a fallback. Raw source labels, source URLs, normalized platform combinations, reported assembly-tool names, structural-tool names/versions and evidence excerpts are retained separately.

`N/A` means the value was not explicit in the checked sources. It does not establish that the technology or tool was not used, or that no information exists elsewhere in NCBI. A failed download is an error, not N/A; all 435 source collections completed without unresolved errors. Raw source responses are cached under `data/inventory/s_cerevisiae_methods_2026-10-30/`.

NCBI GFF source semantics: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/ . Structural annotation pipeline context: https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process/ .

## Sequencing normalization

| Reported sequencing group | Assemblies, all 435 |
|---|---:|
| Illumina + Nanopore | 272 |
| Illumina | 101 |
| PacBio, mode unspecified | 17 |
| Nanopore | 10 |
| 454 | 10 |
| Illumina + PacBio, mode unspecified | 8 |
| MGI/DNBSEQ | 7 |
| N/A | 7 |
| Illumina + 454 | 2 |
| Sanger/ABI capillary | 1 |

Group counts are mutually exclusive and sum to 435. The plot has N/A n=6 because the excluded Partial assembly has unreported sequencing technology. Hybrid assemblies are retained as combinations rather than arbitrarily assigned to one platform.

- Solexa and the literal `llumina` spelling variant are normalized to Illumina; original text remains available.
- Nanopore aliases include explicitly reported Oxford Nanopore, ONT, MinION, PromethION and GridION.
- MGI/DNBSEQ requires an explicit MGI or DNBSEQ label.
- Hi-C is a library/contact protocol and has its own explicit-report field; it is not inferred from an Illumina instrument or a scaffolding-tool name.
- HiFi/CCS requires explicit HiFi/CCS/circular-consensus wording in sequencing metadata. PacBio/Sequel/RS II/Revio or an assembly-tool name such as Hifiasm is not sufficient.
- **No HiFi/CCS or Hi-C use was explicitly reported in the checked sources for this cohort.** The corresponding per-record fields are N/A, not a claim of confirmed non-use.
- Linked SRA assays were not automatically attributed to the assembly. A matching BioSample may contain experiments unrelated to the reads used to assemble that genome.

## Assembly methods and tools

Assembly method is available for 427 records and N/A for 8. Original labels preserve versions and method descriptions. Normalized tool names include assemblers as well as any scaffolders/polishers explicitly listed in that method. Tool counts can overlap and do not establish execution order or the primary assembler.

Most commonly reported names: Canu 218, ABySS 97, Velvet 94, SMARTdenovo 92 and Pilon 78. Pilon is reported within assembly methods and must not be interpreted as a standalone de novo assembler. Other names include Flye, MIRA, SPAdes, CLC Genomics Workbench, Newbler, MaSuRCA, Hifiasm, HGAP, Unicycler, Medaka, RaGOO, Racon, DBG2OLC, Sprai, Phrap and Maia.

## Structural annotation only

| Explicit tool/pipeline | Records | Reported version |
|---|---:|---|
| AUGUSTUS | 5 | 3.3.3 |
| YGAP | 1 | N/A |
| Funannotate | 1 | 1.5.1-1e7a090 |
| N/A | 428 | N/A |

The five AUGUSTUS records explicitly say `Annotated by AUGUSTUS ver.3.3.3` in their first GenBank header. YGAP is explicitly named in an Annotation Pipeline field. Funannotate is explicitly named in an annotation comment; the pipeline is recorded, but its underlying structural predictors are **not inferred**. Funannotate can perform multiple annotation tasks; no functional subtools are imported into the structural-tool list.

Annotation provider names and generic GFF sources such as Genbank, EMBL, DDBJ, tpg and SGD are not converted to software names. Gene products, protein function descriptions, GO terms, Pfam, BLAST, InterPro and eggNOG results are not used to assign structural tools. Since public source metadata is sparse, N/A is expected for many records despite the existence of GFF/CDS/protein files.

## Deliverables and reproduction

In `outputs/s_cerevisiae_methods_20261030/`:

- `s_cerevisiae_genome_metadata_with_methods.xlsx`: original Overview and Assemblies worksheets retained, extra method columns added to the existing table, and Methods comparison worksheet added.
- `genome_metadata_with_methods.csv`: all 435 enriched records.
- `structural_annotation_tools.csv`: structural tools, versions, evidence and NCBI URLs, including N/A rows.
- `sequencing_group_summary.csv`: sequencing-group medians and quartiles for the 434 plotted assemblies.
- `assembly_tool_counts.csv`: overlapping explicit tool-name counts.
- `sequencing_technology_boxplots.png` / `.pdf`: descriptive comparisons, all outliers retained.

```sh
python3 scripts/collect_methods_metadata.py
.venv/bin/python scripts/analyze_methods_metadata.py
.venv/bin/python -m unittest discover -s tests -p 'test_method_classification.py' -v
```

The source collector uses cached successful responses. Use a fresh cache directory or update the snapshot workflow explicitly for a new-date inventory; do not silently mix dates. The workbook uses the bundled artifact runtime and its extension builder. Existing per-assembly metric values and the original Overview were not recomputed or replaced. Plots and summary statistics are static dated outputs, not dynamic workbook formulas.

The comparisons are descriptive. Sequencing platform is confounded with submitter, annotation method, genome representation and sampling. The plots do not isolate a causal effect of sequencing technology.
