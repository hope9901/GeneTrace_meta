# NCBI S. cerevisiae annotation inventory

Checked on 2026-10-30 (KST). Scope: current/latest live assembly records with `species_taxid=4932`, including strain taxids classified under S. cerevisiae. Historical/replaced/suppressed versions and other species/hybrids outside this species taxid are not included.

## Result

**435 distinct assemblies list all three required files**, after removing the GenBank/RefSeq paired duplicate. **434 represent full genomes**; the remaining assembly is partial. These are assembly counts, not deduplicated biological strain counts.

| Assembly level | All eligible assemblies | Full representation only |
|---|---:|---:|
| Complete Genome | 18 | 18 |
| Chromosome | 114 | 113 |
| Scaffold | 95 | 95 |
| Contig | 208 | 208 |
| Total | 435 | 434 |

`Full` genome representation means the submission represents the full genome; it is distinct from the `Complete Genome` assembly level. A Contig assembly may still have Full representation.

The partial assembly is GCA_000091065.2 (S288C). The one paired duplicate is GCA_000146045.2 / GCF_000146045.2 (S288C). Keeping only assemblies other than the S288C reference yields 433 full-representation assemblies; these still need strain/provenance checks before being treated as independent isolates.

## Method

1. Downloaded NCBI current GenBank and RefSeq fungi assembly summaries:
   - https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt
   - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
2. Filtered to species_taxid 4932 and latest version status: 1,769 GenBank records and 1 RefSeq record.
3. Retrieved each of the **1,770 assembly-specific md5checksums.txt manifests**, with retries and persisted responses. No retrieval errors remained.
4. Required exact matching basenames for all three assembly-level files:
   - `<assembly>_genomic.gff.gz` — GFF3 annotation
   - `<assembly>_cds_from_genomic.fna.gz` — nucleotide CDS FASTA
   - `<assembly>_protein.faa.gz` — protein FASTA
5. Found 435 eligible GenBank records plus 1 eligible RefSeq record. Deduplicated their paired assembly to obtain 435 unique assemblies.

This verifies NCBI file-manifest availability, not biological annotation accuracy, completeness, consistency, or the validity of every file's content. Genome/annotation contents were not downloaded for the inventory. The count should not be interpreted as 435 experimentally validated gene annotations.

NCBI's file descriptions: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/ .

## Outputs and reproduction

- `NCBI_S_CEREVISIAE_ANNOTATED_2026-10-30.tsv`: 435 unique assemblies, using GenBank accession as the primary identifier and keeping the paired RefSeq identifier where present.
- `NCBI_S_CEREVISIAE_COUNTS_2026-10-30.json`: counts, snapshot hashes and time.
- `../data/inventory/s_cerevisiae_2026-10-30/`: all filtered metadata, raw file manifests and full audit (excluded from Git).
- `../scripts/count_ncbi_annotations.py`: reproducible audit, requiring downloaded fungi summary files. Run with `--genbank`, `--refseq`, and `--output`; use a fresh output directory for a fresh remote audit because existing per-accession manifest responses are cached.

The taxonomy and annotation availability can change. This is a dated snapshot, not a permanent count.
