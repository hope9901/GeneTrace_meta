#!/usr/bin/env python3
"""Scan complete proteomes and bacterial six-frame ORFs with a pinned Pfam panel."""
import argparse,collections,csv,datetime,json,pathlib,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff,write_fasta,write_json,sha256
from genetrace.domains import scan,six_frame_orfs,domain_coordinates,genomic_assessment,matching_cds
OUT=ROOT/'outputs/domain_bacteria_20261101';RUN=ROOT/'runs/domain_bacteria_20261101'

def table(path,rows):
 if not rows:return
 with path.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--all-yeast',action='store_true');a=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True)
 bacteria=json.loads((ROOT/'data/bacteria/cohort.json').read_text());yeast=json.loads((ROOT/'outputs/orthology_20261030/selected_assemblies.json').read_text());models=json.loads((ROOT/'data/domain_models/manifest.json').read_text());panel=[m['accession'].split('.')[0] for m in models]
 pilot={'GCA_000146045.2','GCA_903819125.2','GCA_903819185.2','GCA_000190255.1','GCA_000326005.1'}
 targets=[dict(assembly_accession=r['assembly_accession'],strain=r['strain'],species='Escherichia coli',source=ROOT/r['source_directory'],bacterial=True) for r in bacteria]
 targets += [dict(assembly_accession=r['assembly_accession'],strain=r['resolved_strain'],species='Saccharomyces cerevisiae',source=ROOT/'data/orthology/source'/r['assembly_accession'],bacterial=False) for r in yeast if a.all_yeast or r['assembly_accession'] in pilot]
 summaries=[];members=[];genomic=[];controls=[];all_arch=[]
 for index,t in enumerate(targets):
  acc=t['assembly_accession'];directory=RUN/acc;directory.mkdir(exist_ok=True)
  write_json(OUT/'status.json',dict(stage='scanning',completed=index,total=len(targets),current=acc,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_yeast=a.all_yeast))
  proteins=read_fasta(t['source']/'protein.faa.gz');proteinfile=directory/'proteins.faa'
  if not proteinfile.exists():write_fasta(proteinfile,proteins)
  hits=scan(proteinfile,directory/'proteome');by_protein=collections.defaultdict(list)
  for h in hits:members.append(dict(assembly_accession=acc,strain=t['strain'],species=t['species'],**h));by_protein[h['protein_id']].append(h)
  architectures=[]
  for pid,hh in by_protein.items():
   ordered=sorted(hh,key=lambda h:(h['aa_start'],h['aa_end'],h['profile']));architectures.append(dict(assembly_accession=acc,protein_id=pid,panel_architecture=' > '.join(h['profile'] for h in ordered),domain_instances=len(ordered),protein_length=len(proteins[pid]),panel_only=True))
  all_arch.extend(architectures)
  for pfam in panel:
   hh=[h for h in hits if h['profile']==pfam];full=[h for h in hh if h['coverage_class']=='near_complete'];summaries.append(dict(assembly_accession=acc,strain=t['strain'],species=t['species'],profile=pfam,proteins_with_domain=len({h['protein_id'] for h in hh}),domain_instances=len(hh),near_complete_instances=len(full),partial_instances=len(hh)-len(full),proteome_records=len(proteins),annotation_scope_warning=(acc=='GCA_903819185.2'),orthology_inferred=False))
  if t['bacterial']:
   sequences=read_fasta(t['source']/'genome.fna.gz');features=read_gff(t['source']/'annotation.gff.gz');orffile=directory/'six_frame_orfs.faa';mapfile=directory/'orf_mapping.json'
   if not orffile.exists():
    orfs,mapping=six_frame_orfs(sequences);write_fasta(orffile,orfs);write_json(mapfile,mapping)
   mapping=json.loads(mapfile.read_text());dh=scan(orffile,directory/'genome_domain')
   assessed=[genomic_assessment(h,domain_coordinates(mapping[h['protein_id']],h['aa_start'],h['aa_end']),features,hits) for h in dh]
   genomic.extend(dict(assembly_accession=acc,strain=t['strain'],**r) for r in assessed);write_json(directory/'genomic_domain_assessments.json',assessed)
   # Structural annotation control: remove a known locus from both proteins and GFF.
   if t['strain']=='MG1655':
    chosen=next(h for h in hits if h['profile']=='PF00005' and h['coverage_class']=='near_complete');pid=chosen['protein_id'];cds=[f for f in features if f['type']=='CDS' and f['attrs'].get('protein_id')==pid]
    if not cds:raise ValueError('Control protein lacks genomic CDS mapping')
    locus=cds[0];control=directory/'annotation_removed_control';control.mkdir(exist_ok=True)
    write_fasta(control/'proteins.faa',{k:v for k,v in proteins.items() if k!=pid});removed_ids={x for f in cds for x in f['attrs'].get('Parent','').split(',') if x}
    keep=[f for f in features if f['attrs'].get('protein_id')!=pid and f['attrs'].get('ID') not in removed_ids]
    with (control/'annotation.gff3').open('w') as f:
     f.write('##gff-version 3\n');f.writelines('\t'.join(x['fields'])+'\n' for x in keep)
    modified=scan(control/'proteins.faa',control/'proteome');local=[h for h in dh if h['profile']=='PF00005' and any(x['attrs'].get('protein_id')==pid for x in matching_cds(domain_coordinates(mapping[h['protein_id']],h['aa_start'],h['aa_end']),features))]
    recovered=[genomic_assessment(h,domain_coordinates(mapping[h['protein_id']],h['aa_start'],h['aa_end']),keep,modified) for h in local]
    outcome=dict(scenario='synthetic_annotation_removal',assembly=acc,removed_protein_id=pid,profile='PF00005',original_protein_has_domain=True,removed_protein_has_domain=any(h['protein_id']==pid for h in modified),genomic_domain_recovered=bool(recovered),recovered_as_annotation_gap=any(r['classification']=='unannotated_genomic_domain_candidate' for r in recovered),source_genome_sha256=sha256(t['source']/'genome.fna.gz'),biological_loss_inferred=False)
    assert outcome['recovered_as_annotation_gap'] and not outcome['removed_protein_has_domain'];controls.append(outcome);write_json(control/'report.json',outcome)
  table(OUT/'domain_counts.csv',summaries);table(OUT/'domain_members.csv',members);table(OUT/'panel_architectures.csv',all_arch);table(OUT/'genomic_domain_candidates.csv',genomic);write_json(OUT/'controls.json',controls)
  print(f'{index+1}/{len(targets)} {acc}: {len(hits)} domain instances',flush=True)
 write_json(OUT/'summary.json',dict(assemblies=len(targets),bacterial_assemblies=len(bacteria),yeast_assemblies=len(targets)-len(bacteria),profiles=models,counts=summaries,controls=controls,all_yeast=a.all_yeast,notes=['All annotated proteins in each included assembly were searched, not a few reference genes.','Six-frame genomic profile searches run only for bacteria, translation table 11 internal codons; alternative starts are not gene predictions.','Only six selected Pfam profiles were searched, not the full Pfam database.','Domains shared by proteins do not imply orthology, identical function or gene copy number.','Genome-level ORF hits are candidates and do not establish an intact gene.','HMM coverage>=0.7 is an extra completeness label; partial GA-passing domains remain in the tables.','Yeast JXXY16.1 has organelle-only annotation; low proteome domain counts cannot establish nuclear gene loss.']))
 write_json(OUT/'status.json',dict(stage='complete',completed=len(targets),total=len(targets),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
if __name__=='__main__':main()
