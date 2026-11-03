# Raw-read linked dataset rebuild — 2026-10-30

279 representative assemblies retained; all 1,116 genome/GFF/CDS/protein files are locally present. Exact BioSample lookup completed for 278 samples with 0 query failures; one reference assembly has no BioSample.

- Same-BioSample WGS available: 216 assemblies, 429 unique WGS runs.
- Any same-BioSample read library: 218 assemblies, 461 unique runs.
- Two assemblies have other libraries only; 60 BioSamples returned no runs. This is not proof that raw reads do not exist.
- FASTQ links total 712.51 GB across all library types, deduplicated by URL.
- Assembly-use provenance verified: 0. FASTQ downloaded: 0. Local read alignment/QC completed: 0.

`assemblies.csv` retains the previous assembly/annotation metadata and adds local file paths, read lookup status and provenance fields. `read_runs.csv` and `read_files.csv` separate run-level metadata from file-level downloads and checksums. `library_platform_summary.csv` separates WGS from RNA-seq and other library types. A PacBio platform or Sequel instrument does not establish HiFi status.

`project_strain_candidates.csv` separately records 8 title-matched read candidates in PRJEB38713, including different-BioSample records. These are NOT counted as exact-BioSample links and have no verified assembly-use provenance. Project-level fallback has only been checked for this project; other unmatched assemblies need further submission/publication investigation.

`investigation_queue.csv` preserves five diagnostic examples for follow-up. The queue does not claim sequencing errors or biological gene loss.

A pilot download of JXXY16.1 PacBio ERR4194697 was attempted but failed because ftp.sra.ebi.ac.uk could not be resolved. No FASTQ was saved or marked as downloaded. Its BioSample differs from the assembly BioSample, so it remains a candidate even if downloaded later.

Previous boxplots and completed results are archived under `back_data/20261030_before_read_qc/`. All 67 archived file SHA256 values were rechecked. Compatibility symlinks preserve prior report references. Active orthology directories remain in place; its archive snapshot is not a live result.

Reproduce: `python3 scripts/rebuild_read_dataset.py` then `python3 scripts/summarize_read_dataset.py`.

Metadata source: [ENA file reports](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html), joined to NCBI assembly BioSample accessions. Responses are retained under `data/read_dataset_20261030/ena_reports/`.
