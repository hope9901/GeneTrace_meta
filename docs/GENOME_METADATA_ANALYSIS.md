# S. cerevisiae assembly metadata and boxplots

Snapshot: 2026-10-30. The cohort is the previously audited 435 GenBank assembly accessions with GFF3, nucleotide CDS FASTA and protein FASTA listed by NCBI. The paired S288C RefSeq record is not counted again. These are assemblies, not necessarily independent strains.

## Outputs

Final deliverables are in `outputs/s_cerevisiae_20261030/`:

- `s_cerevisiae_genome_metadata.xlsx`: Overview (summary, static boxplots and definitions) and Assemblies (435 rows, filterable table, typed values, fixed headers).
- `genome_metadata.csv`: all 435 records, one field per column, UTF-8 with BOM for Excel.
- `assembly_metadata_boxplots.png` and `.pdf`: four panels showing FASTA size, recomputed GC%, actual protein FASTA records, and NCBI contig count by assembly level.
- `assembly_level_summary.csv`: n, minimum, Q1, median, Q3, maximum, mean and missing count for each plotted metric and group.
- `analysis_notes.md`: compact results and definitions.

All 435 assemblies are retained in the table. The plots use 434 `Full` representation assemblies; partial GCA_000091065.2 remains explicitly marked as excluded. `Full` representation differs from `Complete Genome` assembly level. There is no outlier exclusion or significance test.

## Measurement definitions

- **Genome size:** total sequence length across all records in the downloaded genomic FASTA, including ambiguous bases. Both bp and Mb (1 Mb = 1,000,000 bp) are provided.
- **NCBI primary genome size:** the original `genome_size` field, retained in its own column. NCBI defines this as total top-level sequence length in the primary assembly. The full downloaded FASTA may include additional units, so this is not always identical to the FASTA total. The difference is recorded without labeling it a data error.
- **GC%:** recomputed from all genomic FASTA sequences as `100 × (G+C)/(A+C+G+T)`. N and other ambiguity codes are excluded from the denominator. This describes the provided FASTA, including any organellar or additional-unit sequence it contains.
- **NCBI rounded GC%:** retained separately. NCBI's summary rounds this field to the nearest 0.5%; it was not used for the GC boxplot.
- **Number of proteins:** actual number of records in `protein.faa.gz`, not an assumed equivalent of protein-coding gene count. Repeated or isoform sequences are not collapsed. The original protein-coding and total gene counts remain separate.
- **Contig/scaffold count:** NCBI primary assembly metadata. Whole-FASTA record count and sequence N50 are separate fields; FASTA sequence N50 is not relabeled contig N50 because records can be scaffolds.

Source definitions: https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt . Source metadata: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt . Every assembly's NCBI page, genome/protein URLs, and compressed-file SHA-256 values are included in the table.

## Validation and observations

- All 435 protein FASTA files and all 435 genome FASTA files were downloaded, checked against NCBI's published compressed-file MD5, decompressed and counted. No failed downloads were treated as zeros.
- No missing genome size, GC%, protein count or assembly-level values in the final table.
- Primary plot counts: 18 Complete Genome, 113 Chromosome, 95 Scaffold, 208 Contig (total 434).
- Overall medians: 12.065082 Mb FASTA span, 38.351610% GC, 5,690 protein records, 18 contigs.
- Protein record count differs from protein-coding gene count in 103 of the 435 assemblies.
- Full representation GCA_903819185.2 has 16 protein records despite an approximately 11.76 Mb genomic FASTA. Files being present does not establish comprehensive annotation; investigate this record before gene-family absence inference.
- Large assemblies up to 57.099456 Mb and protein counts up to 27,828 are retained. No biological cause is inferred. All source `assembly_type` fields read `haploid`; these labels are not independent ploidy measurements.
- Summary statistics were calculated with NumPy. Boxes span Q1–Q3, median lines are shown, and whiskers use 1.5×IQR. Every assembly is also shown as a jittered point using a fixed random seed (4932). The contig-count axis is logarithmic.

## Reproduction

The inventory and per-accession count/checksum cache live under `data/inventory/` (excluded from Git). Protein/genome FASTA bytes were processed in memory; the count cache retains provenance rather than the full sequence files. Existing successful cache records are reused only if their expected MD5 matches the audited source manifest.

```sh
python3 scripts/collect_genome_metadata.py --output data/inventory/s_cerevisiae_metadata_2026-10-30 --recompute-gc
```

Plotting requires the analysis dependencies in `requirements-analysis.txt` (the development run used Python 3.12):

```sh
.venv/bin/python scripts/plot_genome_metadata.py --input data/inventory/s_cerevisiae_metadata_2026-10-30/genome_metadata.json --output outputs/s_cerevisiae_20261030
```

The workbook was authored with the bundled `@oai/artifact-tool` runtime. Summaries are a dated, externally calculated snapshot and the plots are embedded static images. Editing workbook data does not refresh them; rerun the analysis and workbook builder.
