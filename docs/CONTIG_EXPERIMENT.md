# ACT1 recovery in two NCBI contig-level S. cerevisiae strains

Experiment date: 2026-10-30. Actual local NCBI BLAST+ 2.17.0 runs, no model calls.

## Question and result

Does the S288C ACT1 protein recover a homologous sequence in other strains assembled at Contig level? **Yes, for both tested assemblies.** Each yields a 372-aa HSP at 100% identity and 99.2% query coverage (query positions 4–375 of 375 aa). The missing three query residues are not evidence of gene loss; this is a single-HSP translated search, not a spliced gene reconstruction.

This experiment supports cross-strain sequence recovery. It does **not** establish naturally occurring annotation errors or recovery of a naturally contig-split ACT1 gene.

## Source assemblies

| Strain | NCBI accession | Level | FASTA sequences | Total bp | Computed sequence N50 |
|---|---|---|---:|---:|---:|
| CLIB215 | GCA_019394775.1 | Contig | 244 | 11,295,577 | 83,712 bp |
| KACC 30049 | GCA_060154965.1 | Contig | 21 | 12,073,384 | 1,090,256 bp |

NCBI assembly reports were downloaded and checked for organism, accession, and `Assembly level: Contig`. Statistics above are computed from the downloaded FASTA; all records are included. Both FASTA files have zero N bases. Assembly level alone does not measure fragmentation: KACC 30049 is far more contiguous than CLIB215.

- CLIB215: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_019394775.1/ . BioProject PRJNA412468, Illumina MiSeq, Newbler 2.7. **NCBI explicitly records reference-guided assembly against GCA_000146045.2**, so this is not wholly independent of the S288C reference assembly.
- KACC 30049: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_060154965.1/ . BioProject PRJNA1468922, PacBio Revio, Hifiasm 0.25.0. The assembly report does not list a reference-guided assembly field.

Neither FTP file listing contains genomic GFF or protein FASTA; downloaded feature counts record zero protein-coding gene features. That is an annotation availability issue, not evidence that these genomes contain zero genes. No CDS or protein annotations were fabricated or copied from S288C.

## Search results

Query: S288C ACT1, NP_116614.1. TBLASTN searches each **entire assembly**, not a reference-derived crop. Search arguments: E-value retrieval threshold 1e-5, SEG enabled, max target sequences 100, explicit tabular alignment columns. Qualifying HSP thresholds match the existing baseline: identity ≥40%, single-HSP query coverage ≥60%, E-value ≤1e-10.

| Strain | Target contig | HSP coordinates, 1-based inclusive | Strand | Identity | Query coverage | Bit score |
|---|---|---|---|---:|---:|---:|
| CLIB215 | PCEV01000060.1 | 44,982–46,097 | + | 100% | 99.2% | 775 |
| KACC 30049 | JBYNJL010000013.1 | 38,225–39,340 | − | 100% | 99.2% | 775 |

Both alignments report E-value `0.0` in BLAST's output formatting; this should not be interpreted as a mathematically exact zero probability. The highest-scoring match is an ACT1 homolog candidate. Orthology, full gene structure, and function are not independently established.

## Controlled perturbations

For each strain, the strongest observed HSP plus 500 bp on each side defines the seeded locus. This location is derived from the original search, so the follow-up experiment is explicitly seeded, not a blind loss benchmark. Original source files are preserved. The full assembly remains the search subject in every scenario.

| Scenario | CLIB215 | KACC 30049 |
|---|---|---|
| Original assembly, no sequence edits | homolog_detected_annotation_unavailable | homolog_detected_annotation_unavailable |
| Seeded region replaced with N | assembly_limited | assembly_limited |
| Seeded region replaced with A | inconclusive | inconclusive |

Two qualifying actin-related hits elsewhere in each assembly remain after both perturbations. They are reported separately; their presence is not used to claim that ACT1 still exists at the seeded locus. This is a meaningful difference from searching only a short cropped locus.

Because genuine GFF/protein inputs are absent, **the prior annotated-versus-annotation-removed test cannot be repeated on these native submissions**. Empty annotation files would falsely turn every detected gene into an annotation-gap candidate. The experiment instead reports annotation unavailable and directly evaluates sequence recovery and conservative handling of masked/missing sequence.

## Reproduce and inspect

Run from the project root with Python 3.9+ and BLAST+ installed:

```sh
python3 scripts/download_data.py
python3 scripts/prepare_demo.py
python3 scripts/download_strains.py
python3 scripts/test_contig_strains.py --output runs/new-contig-experiment
python3 -m unittest discover -s tests -v
```

`docs/contig-sources.json` pins original SHA-256 checksums and URLs. Downloads fail if the snapshot changes. Existing run directories are not overwritten. Recorded first results are under `runs/contig-strains-001/`: each strain/scenario includes the exact genome searched, query, tool version, arguments, raw and parsed alignments, and local assessment. `comparison.json` combines results. Outputs under `runs` and source downloads are excluded from Git.

Tests include fresh execution of the original four S288C controls plus validation of stored two-strain BLAST evidence. The stored-result test expects the first run path and skips if absent; it does not rerun that experiment by itself.

## Implication for GeneTrace

Sequence recovery transfers to these two other strains, including the more fragmented CLIB215 assembly. The next scientifically distinct challenge is a locus that is genuinely split at contig ends, followed by flanking-gene localization and chained/spliced alignment. Native annotation error detection still requires an annotated draft assembly or an explicitly labeled annotation pipeline.
