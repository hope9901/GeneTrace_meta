"""Combine calibrated searches conservatively and summarize diagnostic evidence."""
import csv
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/genetrace-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/genetrace-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from genetrace.io import read_gff
OUT=ROOT/'outputs/loss_artifact_tests_20261030'
RUN=ROOT/'runs/loss_artifact_tests_20261030'


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    default=json.loads((OUT/'natural_default_reports.json').read_text()) if (OUT/'natural_default_reports.json').exists() else json.loads((OUT/'natural_case_reports.json').read_text())
    sensitivity=json.loads((OUT/'natural_sensitivity_reports.json').read_text())
    default_callable={r['reference_query_id'] for r in default[0]['results'] if r['protein_homolog_found'] and r['genome_homolog_found']}
    sensitive_callable={r['reference_query_id'] for r in sensitivity[0]['results'] if r['protein_homolog_found'] and r['genome_homolog_found']}
    if len(sensitive_callable)!=64:raise ValueError('Sensitivity panel is not calibrated against the reference')
    by_default={r['assembly_accession']:r for r in default}
    scope={r['assembly_accession']:r for r in json.loads((OUT/'annotation_scope_audit.json').read_text())}
    all_rows=[];summary=[]
    for report in sensitivity:
        accession=report['assembly_accession']; defaults={r['reference_query_id']:r for r in by_default[accession]['results']}
        features=read_gff(ROOT/'data/orthology/source'/accession/'annotation.gff.gz')
        evidence=json.loads((RUN/'sensitivity_no_mask'/accession/'evidence.json').read_text())
        genome_raw=next(e['result']['hits'] for e in evidence if e['tool']=='tblastn_parsed')
        rows=[]
        for current in report['results']:
            row=dict(current);query=row['reference_query_id'];old=defaults[query]
            use_default=query in default_callable
            protein=row['protein_homolog_found'] or (use_default and old['protein_homolog_found'])
            genome=row['genome_homolog_found'] or (use_default and old['genome_homolog_found'])
            if not row['genome_homolog_found'] and genome:
                for key in ('best_genome_seqid','best_genome_start','best_genome_end','strand','genome_identity','genome_query_coverage','overlapping_cds'):
                    row[key]=old[key]
            partial=[h for h in genome_raw if h['qseqid']==query and h['pident']>=70 and h['query_coverage']>=.5 and h['evalue']<=1e-20]
            best_partial=max(partial,key=lambda h:h['bitscore']) if partial else None
            row.update(protein_homolog_found=bool(protein),genome_homolog_found=bool(genome),
                       default_reference_callable=use_default,
                       settings_disagree=(old['protein_homolog_found']!=current['protein_homolog_found'] or old['genome_homolog_found']!=current['genome_homolog_found']),
                       partial_or_divergent_genome_evidence=bool(best_partial) and not genome,
                       weaker_hit_identity=best_partial['pident'] if best_partial and not genome else '',
                       weaker_hit_query_coverage=best_partial['query_coverage'] if best_partial and not genome else '',
                       weaker_hit_seqid=best_partial['sseqid'] if best_partial and not genome else '',
                       weaker_hit_start=best_partial['start'] if best_partial and not genome else '',
                       weaker_hit_end=best_partial['end'] if best_partial and not genome else '',
                       weaker_hit_strand=best_partial['strand'] if best_partial and not genome else '',
                       calibrated_criteria='>=90% identity, >=90% single-HSP coverage, E<=1e-20; union of reference-callable settings')
            context=[]
            if genome:
                context=[f for f in features if f['type']=='gene' and f['seqid']==row['best_genome_seqid'] and f['strand']==row['strand']
                         and max(0,min(f['end'],row['best_genome_end'])-max(f['start'],row['best_genome_start'])+1)
                         /(row['best_genome_end']-row['best_genome_start']+1)>=.5]
            notes=' | '.join(f['attrs'].get('Note','') for f in context if f['attrs'].get('Note'))
            biotypes=';'.join(sorted({f['attrs'].get('gene_biotype','unspecified') for f in context}))
            disrupted=bool(re.search(r'frameshift|interrupted|premature|stop codon|pseudogene',notes+' '+biotypes,re.I)
                           or any(f['attrs'].get('pseudo')=='true' for f in context))
            variations=[f for f in features if genome and f['type'] in ('sequence_alteration','variation')
                        and f['seqid']==row['best_genome_seqid'] and f['start']<=row['best_genome_end'] and f['end']>=row['best_genome_start']]
            row.update(overlapping_gene_ids=';'.join(f['attrs'].get('ID','') for f in context),
                       overlapping_gene_biotypes=biotypes,overlapping_gene_notes=notes,
                       recorded_disruption_or_pseudogene=disrupted,overlapping_sequence_alterations=len(variations),
                       sequence_alteration_notes=';'.join(sorted({f['attrs'].get('Note','') for f in variations})))
            if genome and not protein:
                if disrupted:row['classification']='disrupted_or_pseudogene_model_unresolved'
                elif row['overlapping_cds']:row['classification']='genomic_homolog_with_protein_model_discrepancy'
                elif context:row['classification']='annotated_gene_without_protein_model_unresolved'
                else:row['classification']='unannotated_genomic_homolog'
            elif genome and protein:row['classification']='protein_and_genome_homolog_detected'
            elif row['partial_or_divergent_genome_evidence']:row['classification']='partial_or_divergent_homolog_unresolved'
            elif protein:row['classification']='protein_homolog_only_genome_search_unresolved'
            else:row['classification']='no_qualifying_homolog_unresolved'
            rows.append(row)
        all_rows.extend(rows)
        summary.append(dict(assembly_accession=accession,strain=report['strain'],assembly_level=report['assembly_level'],
                  panel_queries=64,source_proteins=report['proteins_in_submission'],
                  protein_positive=sum(r['protein_homolog_found'] for r in rows),
                  full_genome_positive=sum(r['genome_homolog_found'] for r in rows),
                  additional_partial_or_divergent=sum(r['partial_or_divergent_genome_evidence'] for r in rows),
                  unannotated_genomic_homologs=sum(r['classification']=='unannotated_genomic_homolog' for r in rows),
                  protein_model_discrepancies=sum(r['classification']=='genomic_homolog_with_protein_model_discrepancy' for r in rows),
                  disrupted_model_candidates=sum(r['classification']=='disrupted_or_pseudogene_model_unresolved' for r in rows),
                  settings_sensitive_queries=sum(r['settings_disagree'] for r in rows),
                  nuclear_coding_loci=scope[accession]['nuclear_coding_loci'],
                  organelle_coding_loci=scope[accession]['organelle_coding_loci'],
                  nuclear_bp=scope[accession]['nuclear_or_unclassified_bp'],
                  biological_loss_calls=0))
    write_csv(OUT/'adjudicated_locus_results.csv',all_rows)
    write_csv(OUT/'adjudicated_case_summary.csv',summary)
    (OUT/'adjudicated_cases.json').write_text(json.dumps(dict(cases=summary,loci=all_rows),indent=2)+'\n')
    modes=dict(default_reference_callable=len(default_callable),sensitivity_reference_callable=len(sensitive_callable),
               baseline_uninformative_queries=[r for r in default[0]['results'] if r['reference_query_id'] not in default_callable],
               evidence_scope='Positive homolog sequence evidence; no functional gene-loss inference or cohort prevalence estimate')
    (OUT/'search_calibration_summary.json').write_text(json.dumps(modes,indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(11.5,5.8));y=np.arange(len(summary));height=.32
    p=[r['protein_positive'] for r in summary];g=[r['full_genome_positive'] for r in summary];partial=[r['additional_partial_or_divergent'] for r in summary]
    ax.barh(y-height/2,p,height,color='#438A83',label='Submitted protein FASTA: strong homolog')
    ax.barh(y+height/2,g,height,color='#355C7D',label='Genomic DNA: strong homolog')
    ax.barh(y+height/2,partial,height,left=g,color='#D1D7DE',hatch='///',label='Genomic DNA: additional partial/divergent evidence')
    for i,(pp,gg,part) in enumerate(zip(p,g,partial)):
        ax.text(pp+.65,i-height/2,str(pp),va='center',fontsize=10)
        ax.text(gg+part+.65,i+height/2,str(gg)+((' + '+str(part)+' partial') if part else ''),va='center',fontsize=10)
    ax.set_yticks(y,[r['strain']+'\n'+r['assembly_level'] for r in summary]);ax.invert_yaxis()
    ax.set_xlim(0,82);ax.set_xticks(range(0,65,8));ax.set_xlabel('Reference queries with homolog evidence (64 tested per assembly)')
    ax.grid(axis='x',color='#E1E5E9');ax.set_axisbelow(True)
    ax.set_title('Apparent protein absence can coexist with retained genomic sequence',loc='left',fontsize=15,fontweight='bold',pad=18)
    ax.legend(loc='upper left',bbox_to_anchor=(0,-.15),frameon=False,fontsize=9)
    fig.text(.22,.025,'Same fixed 64-locus panel across five unmodified assemblies; calibrated on S288C.\n'
             'Selected diagnostic cases. Homology is not proof of orthology, intact function, or biological loss.',fontsize=9,linespacing=1.5)
    fig.subplots_adjust(left=.22,right=.98,top=.88,bottom=.32)
    fig.savefig(OUT/'protein_vs_genome_evidence.png',dpi=190,facecolor='white')
    fig.savefig(OUT/'protein_vs_genome_evidence.pdf',facecolor='white');plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
