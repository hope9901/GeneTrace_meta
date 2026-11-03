#!/usr/bin/env python3
"""Keep publication-corroborated links separate from exact-BioSample joins."""
import csv,hashlib,json,pathlib,re,sys
root=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from genetrace.io import read_fasta,read_gff

def main():
 data=root/'data/read_dataset_20261030';out=root/'outputs/read_dataset_20261030'
 cohort=json.loads((out/'assemblies.json').read_text());long=list(csv.DictReader((data/'project_candidates/PRJEB38713.tsv').open(),delimiter='\t'));short=list(csv.DictReader((data/'provenance/PRJNA396809_aliases.tsv').open(),delimiter='\t'))
 # Table S1 and S5, visually checked PDF page 2 and page 6; no run accession is given in the paper.
 expected={'Y55':(11791742,80847),'BJ4':(11888108,86123),'HLJ1':(11702338,78472),'HN1':(11820137,88103),'SX2':(11918758,80393),'EM14S01-3B':(11803938,86728),'JXXY16.1':(11692380,68736),'XXYS1.4':(11654158,85245)}
 links=[];metrics=[]
 for row in cohort:
  strain=row['resolved_strain']
  if row['bioproject']!='PRJEB38713' or strain not in expected:continue
  acc=row['assembly_accession'];seqs=read_fasta(root/'data/orthology/source'/acc/'genome.fna.gz')
  # Explicit organelle labels only; unclassified records are never assumed absent.
  features=read_gff(root/'data/orthology/source'/acc/'annotation.gff.gz')
  mito={f['seqid'] for f in features if f['type']=='region' and any(f['attrs'].get(k,'').lower().startswith('mitochondrion') for k in ('organelle','genome'))}
  observed_mito=sum(len(seq) for key,seq in seqs.items() if key in mito);observed_nuclear=sum(map(len,seqs.values()))-observed_mito
  matching=(observed_nuclear,observed_mito)==expected[strain]
  metrics.append(dict(assembly_accession=acc,strain=strain,paper_nuclear_bp=expected[strain][0],observed_non_mito_bp=observed_nuclear,paper_mito_bp=expected[strain][1],observed_mito_bp=observed_mito,exact_length_match=matching))
  selected=[('PacBio',r,'project_and_exact_strain_in_title') for r in long if re.search(r'(?<![A-Za-z0-9])'+re.escape(strain)+r'(?![A-Za-z0-9])',r['sample_title'],re.I)]
  selected += [('Illumina',r,'paper_reused_project_and_exact_sample_alias') for r in short if r['sample_alias']==strain and strain not in ('EM14S01-3B','Y55')]
  for platform,r,basis in selected:
   links.append(dict(assembly_accession=acc,strain=strain,assembly_biosample=row['biosample'],read_biosample=r['sample_accession'],run_accession=r['run_accession'],platform=platform,instrument=r['instrument_model'],library_strategy=r['library_strategy'],read_project='PRJEB38713' if platform=='PacBio' else 'PRJNA396809',exact_biosample_match=row['biosample']==r['sample_accession'],paper_assembly_lengths_match=matching,evidence_basis=basis+'; publication Table S5',link_status='publication_corroborated_run_candidate',exact_run_use_in_assembly_version_verified=False,fastq_ftp=r['fastq_ftp'],fastq_md5=r['fastq_md5'],fastq_bytes=r['fastq_bytes'],paper_doi='10.1093/gbe/evab001',paper_source='https://www.pure.ed.ac.uk/ws/portalfiles/portal/480257147/evab001_supplementary_data.pdf'))
 for filename,rows in [('publication_supported_read_links.csv',links),('publication_assembly_metric_checks.csv',metrics)]:
  with (out/filename).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 for row in cohort:
  linked=[x for x in links if x['assembly_accession']==row['assembly_accession']]
  row['publication_supported_run_count']=len(linked)
  row['publication_supported_run_candidates']=';'.join(x['run_accession'] for x in linked) or 'N/A'
  row['publication_supported_platforms']=';'.join(sorted({x['platform'] for x in linked})) or 'N/A'
  row['publication_exact_run_version_use_verified']=False
 (out/'assemblies.json').write_text(json.dumps(cohort,indent=2)+'\n')
 with (out/'assemblies.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(cohort[0]));writer.writeheader();writer.writerows(cohort)
 provenance={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (data/'provenance').glob('*') if p.is_file()}
 (out/'publication_provenance.json').write_text(json.dumps(dict(sources_sha256=provenance,table_s5_pdf_page=6,table_s1_pdf_page=2,manually_checked=True,links=len(links),assemblies=len(metrics),assembly_metric_exact_matches=sum(r['exact_length_match'] for r in metrics),limitations=['Table S5 gives strain and project, not run accessions or exact assembly version.','Run candidate selection uses exact strain alias/title, not substring matching.','Reference S288C reads in this paper do not establish raw-read provenance for GCA_000146045.2.','Y55 and EM14S01-3B reused short-read projects have not yet been surveyed.']),indent=2)+'\n')
 print(len(links),'publication-corroborated run candidates;',sum(x['exact_length_match'] for x in metrics),'of',len(metrics),'assembly length matches')
if __name__=='__main__':main()
