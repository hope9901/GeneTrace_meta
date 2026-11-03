#!/usr/bin/env python3
"""Validate and describe the read-linked cohort without inferring assembly provenance."""
import csv,hashlib,json,pathlib,collections,re
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
 out=ROOT/'outputs/read_dataset_20261030';summary=json.loads((out/'summary.json').read_text());rows=json.loads((out/'assemblies.json').read_text())
 runs=list(csv.DictReader((out/'read_runs.csv').open()));files=list(csv.DictReader((out/'read_files.csv').open()))
 assert len(rows)==279 and len({r['assembly_accession'] for r in rows})==279
 for r in rows:
  for kind in ['genome','gff','cds','protein']:assert (ROOT/r[kind+'_local_path']).is_file()
 for r in runs:assert r['assembly_biosample'] in (r['sample_accession'],r['secondary_sample_accession'])
 for f in files:assert re.fullmatch('[0-9a-f]{32}',f['md5']) and int(f['bytes'])>0
 archive=ROOT/'back_data/20261030_before_read_qc';m=json.loads((archive/'archive_manifest.json').read_text())
 for path,digest in m['files'].items():assert hashlib.sha256((archive/path).read_bytes()).hexdigest()==digest,path
 counts=collections.Counter((r['library_strategy'],r['library_source'],r['instrument_platform']) for r in runs)
 with (out/'library_platform_summary.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['library_strategy','library_source','platform','run_count']);w.writerows((*k,v) for k,v in sorted(counts.items()))
 probe=ROOT/'data/read_dataset_20261030/project_candidates/PRJEB38713.tsv'
 candidates=[]
 if probe.exists():
  rr=list(csv.DictReader(probe.open(),delimiter='\t'))
  for r in rows:
   if r['bioproject']!='PRJEB38713':continue
   label=r['resolved_strain']
   if not label or label=='N/A':continue
   for read in rr:
    if re.search(r'(?<![A-Za-z0-9])'+re.escape(label)+r'(?![A-Za-z0-9])',read['sample_title'],re.I):
     candidates.append(dict(assembly_accession=r['assembly_accession'],strain=label,assembly_biosample=r['biosample'],read_biosample=read['sample_accession'],run_accession=read['run_accession'],platform=read['instrument_platform'],library_strategy=read['library_strategy'],sample_title=read['sample_title'],evidence='same_project_exact_strain_label_in_title',biosample_match=r['biosample']==read['sample_accession'],assembly_use_verified=False,fastq_ftp=read['fastq_ftp'],fastq_bytes=read['fastq_bytes'],fastq_md5=read['fastq_md5']))
  if candidates:
   with (out/'project_strain_candidates.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(candidates[0]));w.writeheader();w.writerows(candidates)
 # A small actionable queue retains diagnostic examples even where direct read links are absent.
 priority={'GCA_903819185.2':'nuclear annotation scope mismatch; investigate ACS1/PSY4','GCA_903819125.2':'YJU3 annotation-gap candidate','GCA_000190255.1':'SWC3 frameshift model; sequence versus assembly unresolved','GCA_000326005.1':'partial/disrupted models; unresolved','GCA_000146045.2':'reference search calibration; raw-read provenance needs separate validation'}
 queue=[]
 for r in rows:
  acc=r['assembly_accession']
  if acc in priority:queue.append(dict(assembly_accession=acc,strain=r['resolved_strain'],reason=priority[acc],raw_read_status=r['raw_read_status'],sample_linked_wgs=r['raw_run_accessions'],project_candidate_runs=';'.join(x['run_accession'] for x in candidates if x['assembly_accession']==acc) or 'N/A',assembly_use_verified=False,next_action='verify publication/submission run provenance before regional read QC'))
 with (out/'investigation_queue.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(queue[0]));w.writeheader();w.writerows(queue)
 report=f'''# Raw-read linked dataset rebuild — 2026-10-30

279 representative assemblies retained; all 1,116 genome/GFF/CDS/protein files are locally present. Exact BioSample lookup completed for 278 samples with {summary['query_failures']} query failures; one reference assembly has no BioSample.

- Same-BioSample WGS available: {summary['assemblies_with_sample_linked_wgs']} assemblies, {summary['unique_wgs_runs']} unique WGS runs.
- Any same-BioSample read library: {summary['assemblies_with_any_sample_linked_reads']} assemblies, {summary['unique_runs']} unique runs.
- Two assemblies have other libraries only; 60 BioSamples returned no runs. This is not proof that raw reads do not exist.
- FASTQ links total {summary['fastq_bytes_unique']/1e9:.2f} GB across all library types, deduplicated by URL.
- Assembly-use provenance verified: 0. FASTQ downloaded: 0. Local read alignment/QC completed: 0.

`assemblies.csv` retains the previous assembly/annotation metadata and adds local file paths, read lookup status and provenance fields. `read_runs.csv` and `read_files.csv` separate run-level metadata from file-level downloads and checksums. `library_platform_summary.csv` separates WGS from RNA-seq and other library types. A PacBio platform or Sequel instrument does not establish HiFi status.

`project_strain_candidates.csv` separately records {len(candidates)} title-matched read candidates in PRJEB38713, including different-BioSample records. These are NOT counted as exact-BioSample links and have no verified assembly-use provenance. Project-level fallback has only been checked for this project; other unmatched assemblies need further submission/publication investigation.

`investigation_queue.csv` preserves five diagnostic examples for follow-up. The queue does not claim sequencing errors or biological gene loss.

A pilot download of JXXY16.1 PacBio ERR4194697 was attempted but failed because ftp.sra.ebi.ac.uk could not be resolved. No FASTQ was saved or marked as downloaded. Its BioSample differs from the assembly BioSample, so it remains a candidate even if downloaded later.

Previous boxplots and completed results are archived under `back_data/20261030_before_read_qc/`. All {len(m['files'])} archived file SHA256 values were rechecked. Compatibility symlinks preserve prior report references. Active orthology directories remain in place; its archive snapshot is not a live result.

Reproduce: `python3 scripts/rebuild_read_dataset.py` then `python3 scripts/summarize_read_dataset.py`.

Metadata source: [ENA file reports](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html), joined to NCBI assembly BioSample accessions. Responses are retained under `data/read_dataset_20261030/ena_reports/`.
'''
 (out/'dataset_report.md').write_text(report)
 (ROOT/'docs/READ_QC_DATASET.md').write_text(report)
 print('Validated 279 assemblies, 1,116 local source files, run joins and archive hashes. Project-label candidates:',len(candidates))
if __name__=='__main__':main()
