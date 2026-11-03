"""Validate all-gene coverage and export cohort-wide orthogroup matrices."""
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import write_json
from genetrace.orthology_evidence import scope_aware_presence

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--results',required=True)
    parser.add_argument('--output',default=str(ROOT/'outputs/orthology_20261030/orthogroups'))
    parser.add_argument('--cohort',default=str(ROOT/'outputs/orthology_20261030/selected_assemblies.json'))
    args=parser.parse_args()
    results=Path(args.results)
    output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
    selected=json.loads(Path(args.cohort).read_text())
    accessions=[r['assembly_accession'] for r in selected]
    og_dir=results/'Orthogroups'
    candidates=[og_dir/'Orthogroups.tsv',og_dir/'MCL_Orthogroups.tsv']
    og_file=next((f for f in candidates if f.exists()),None)
    if og_file is None:
        raise FileNotFoundError('No completed OrthoFinder orthogroup membership file')
    # Only interpret stage output here; -og is not gene-tree orthology inference.
    is_mcl=og_file.name.startswith('MCL_')
    unassigned=og_dir/('MCL_Orthogroups_UnassignedGenes.tsv' if is_mcl else 'Orthogroups_UnassignedGenes.tsv')
    if not unassigned.exists(): unassigned=None
    maps={}
    audit_dir=ROOT/'outputs/orthology_20261030/quality_control'
    scope_file=audit_dir/'annotation_scope_audit.json'
    labels_file=audit_dir/'sequence_compartments.json'
    scopes={r['assembly_accession']:r for r in json.loads(scope_file.read_text())} if scope_file.exists() else {}
    compartment_labels=json.loads(labels_file.read_text()) if labels_file.exists() else {}
    expected_by_assembly=Counter()
    for accession in accessions:
        with (ROOT/'data/orthology/gene_maps'/(accession+'.tsv')).open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                if row['sequence_id'] in maps:
                    raise ValueError('Duplicate exported sequence ID')
                maps[row['sequence_id']]=row
                row['sequence_compartment']=compartment_labels.get(accession,{}).get(row['seqid'],'unknown')
                expected_by_assembly[accession]+=1
    seen=set();counts=Counter();stats=[]
    member_fields=['orthogroup','assignment_type','sequence_id','assembly_accession','gene_id','ncbi_protein_id',
                   'all_isoform_protein_ids','mapping_type','gene_name','locus_tag','seqid','start','end','strand','sequence_compartment']
    with (output/'gene_membership.tsv').open('w',newline='') as member_handle, \
         (output/'copy_number.tsv').open('w',newline='') as copy_handle, \
         (output/'presence_absence.tsv').open('w',newline='') as presence_handle, \
         (output/'presence_absence_scope_aware.tsv').open('w',newline='') as scoped_handle:
        member_writer=csv.DictWriter(member_handle,fieldnames=member_fields,delimiter='\t',extrasaction='ignore')
        member_writer.writeheader()
        copy_writer=csv.writer(copy_handle,delimiter='\t');copy_writer.writerow(['Orthogroup']+accessions)
        presence_writer=csv.writer(presence_handle,delimiter='\t');presence_writer.writerow(['Orthogroup']+accessions)
        scoped_writer=csv.writer(scoped_handle,delimiter='\t');scoped_writer.writerow(['Orthogroup']+accessions)
        def consume(filename,assignment_type):
            if filename is None:
                return
            with filename.open() as handle:
                reader=csv.DictReader(handle,delimiter='\t')
                species=[k for k in reader.fieldnames if k!='Orthogroup']
                if set(species)!=set(accessions):
                    raise ValueError('OrthoFinder species columns do not match the selected cohort')
                for record in reader:
                    group=record['Orthogroup']
                    row_counts=[]
                    compartments=set()
                    for accession in accessions:
                        tokens=[t.strip() for t in record[accession].split(',') if t.strip()]
                        row_counts.append(len(tokens))
                        for token in tokens:
                            # OrthoFinder normally prefixes species unless -X was used.
                            if token not in maps:
                                prefix=accession+'_'
                                if token.startswith(prefix): token=token[len(prefix):]
                            if token not in maps or maps[token]['assembly_accession']!=accession:
                                raise ValueError('Unknown or mismatched gene ID: '+token)
                            if token in seen:
                                raise ValueError('Gene assigned more than once: '+token)
                            seen.add(token);counts[accession]+=1
                            compartments.add(maps[token]['sequence_compartment'])
                            member_writer.writerow(dict(maps[token],orthogroup=group,assignment_type=assignment_type))
                    copy_writer.writerow([group]+row_counts)
                    presence_writer.writerow([group]+[int(n>0) for n in row_counts])
                    scoped=[scope_aware_presence(n,compartments,scopes.get(accession,{})) for accession,n in zip(accessions,row_counts)]
                    scoped_writer.writerow([group]+scoped)
                    present=sum(n>0 for n in row_counts)
                    stats.append(dict(orthogroup=group,assignment_type=assignment_type,
                                      total_coding_units=sum(row_counts),genomes_present=present,
                                      genomes_absent=len(accessions)-present,prevalence=present/len(accessions),
                                      single_copy_genomes=sum(n==1 for n in row_counts),
                                      multicopy_genomes=sum(n>1 for n in row_counts),max_copies=max(row_counts),
                                      present_in_all=present==len(accessions),
                                      single_copy_in_all=all(n==1 for n in row_counts)))
                    stats[-1].update(sequence_compartments=';'.join(sorted(compartments)),
                                     scope_unavailable_genomes=sum(x=='N/A' for x in scoped))
        consume(og_file,'MCL_orthogroup' if is_mcl else 'phylogenomic_orthogroup')
        consume(unassigned,'unassigned_by_OrthoFinder')
    missing=set(maps)-seen
    if missing:
        write_json(output/'unaccounted_gene_ids.json',sorted(missing))
        raise ValueError('%d input coding units not accounted for in group/unassigned outputs'%len(missing))
    if counts!=expected_by_assembly:
        raise ValueError('Per-assembly gene accounting failed')
    with (output/'orthogroup_summary.csv').open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(stats[0]));writer.writeheader();writer.writerows(stats)
    summary=dict(status='complete',assemblies=len(accessions),coding_units=len(maps),accounted_coding_units=len(seen),
                 orthogroups=sum(r['assignment_type']!='unassigned_by_OrthoFinder' for r in stats),
                 unassigned_rows=sum(r['assignment_type']=='unassigned_by_OrthoFinder' for r in stats),
                 ubiquitous_groups=sum(r['present_in_all'] for r in stats),
                 single_copy_in_all_groups=sum(r['single_copy_in_all'] for r in stats),
                 biological_loss_inferred=False,
                 method='OrthoFinder 3.2, DIAMOND, MCL' if is_mcl else 'OrthoFinder 3.2, DIAMOND, MCL, FAMSA, FastTree, ASTRAL-Pro; scalable core+assign',
                 assignment_stage='MCL clusters only' if is_mcl else 'Final phylogenomic orthogroups',
                 scope_aware_matrix_applied=bool(scopes and compartment_labels),
                 sparse_annotation_assemblies=[r['assembly_accession'] for r in selected if r.get('number_of_proteins',0)<1000],
                 sparse_annotation_caveat='All retained strains included; observed zero cells are not validated biological gene absence')
    write_json(output/'summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
