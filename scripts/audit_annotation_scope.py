"""Audit compartment scope and assembly-level confounding, without touching OF inputs."""
import csv
import gzip
import json
import os
import re
import sys
from collections import Counter,defaultdict
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/genetrace-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/genetrace-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import read_gff,write_json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/loss_artifact_tests_20261030'
SOURCE=ROOT/'data/orthology/source'
LEVELS=['Complete Genome','Chromosome','Scaffold','Contig']
COLORS=['#355C7D','#438A83','#B6873D','#8064A2']


def compartment(text):
    text=text.lower()
    if re.search(r'mitochondri|plastid|chloroplast|chromosome[=: ]+(?:mt|mito)\b',text):return 'organelle'
    if 'plasmid' in text:return 'plasmid'
    return 'nuclear_or_unclassified'


def csv_write(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def audit(row):
    accession=row['assembly_accession']; source=SOURCE/accession
    features=read_gff(source/'annotation.gff.gz')
    labels={}
    for f in features:
        if f['type']=='region':
            kind=compartment(' '.join(k+'='+v for k,v in f['attrs'].items()))
            if kind!='nuclear_or_unclassified':labels[f['seqid']]=kind
    counts=defaultdict(Counter); current=None
    with gzip.open(source/'genome.fna.gz','rt') as handle:
        for line in handle:
            if line.startswith('>'):
                seqid=line[1:].split()[0]; kind=labels.get(seqid,compartment(line))
                labels[seqid]=kind; current=seqid; counts[kind]['records']+=1
            else:
                seq=line.strip().upper()
                if seq and current is None:raise ValueError('Invalid genome FASTA')
                if not seq:continue
                count=counts[labels[current]]
                count['bases']+=len(seq);count['gc']+=seq.count('G')+seq.count('C')
                count['acgt']+=sum(seq.count(b) for b in 'ACGT')
    coding=Counter()
    with (ROOT/'data/orthology/gene_maps'/(accession+'.tsv')).open() as handle:
        for gene in csv.DictReader(handle,delimiter='\t'):
            coding[labels.get(gene['seqid'],'missing_genome_sequence')]+=1
    nuclear=counts['nuclear_or_unclassified']; other=counts['organelle']+counts['plasmid']
    all_gc=sum(c['gc'] for c in counts.values());all_acgt=sum(c['acgt'] for c in counts.values())
    if not nuclear['acgt']:raise ValueError('No retained putative nuclear sequence: '+accession)
    result=dict(assembly_accession=accession,resolved_strain=row['resolved_strain'],assembly_level=row['assembly_level'],
                annotation_provider=row['annotation_provider'],sequencing_group=row['sequencing_group'],
                bioproject=row.get('bioproject') or 'N/A',number_of_proteins=row['number_of_proteins'],
                genome_fasta_bp=sum(c['bases'] for c in counts.values()),
                nuclear_or_unclassified_bp=nuclear['bases'],explicit_organelle_bp=counts['organelle']['bases'],
                explicit_plasmid_bp=counts['plasmid']['bases'],gc_all_percent=100*all_gc/all_acgt,
                gc_organelle_excluded_percent=100*nuclear['gc']/nuclear['acgt'],
                nuclear_coding_loci=coding['nuclear_or_unclassified'],organelle_coding_loci=coding['organelle'],
                plasmid_coding_loci=coding['plasmid'],
                coding_loci_on_missing_genome_sequence=coding['missing_genome_sequence'],
                nuclear_coding_loci_per_mb=coding['nuclear_or_unclassified']/(nuclear['bases']/1e6),
                annotation_scope_mismatch=coding['nuclear_or_unclassified']==0 and coding['organelle']>0 and nuclear['bases']>1e6,
                gc_scope_shift_percentpoints=100*all_gc/all_acgt-100*nuclear['gc']/nuclear['acgt'])
    if result['genome_fasta_bp']!=row['genome_size_bp'] or abs(result['gc_all_percent']-row['gc_percent'])>1e-9:
        raise ValueError('Audit changed the original whole-FASTA measurement: '+accession)
    return result,labels


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=json.loads((ROOT/'outputs/orthology_20261030/selected_assemblies.json').read_text())
    audit_rows=[];all_labels={}
    for i,row in enumerate(rows,1):
        result,labels=audit(row);audit_rows.append(result);all_labels[row['assembly_accession']]=labels
        if i%40==0:print('Scope audits: %d/%d'%(i,len(rows)),flush=True)
    write_json(OUT/'annotation_scope_audit.json',audit_rows)
    write_json(OUT/'sequence_compartments.json',all_labels)
    csv_write(OUT/'annotation_scope_audit.csv',audit_rows)
    qc=ROOT/'outputs/orthology_20261030/quality_control'
    qc.mkdir(parents=True,exist_ok=True)
    write_json(qc/'annotation_scope_audit.json',audit_rows)
    write_json(qc/'sequence_compartments.json',all_labels)
    csv_write(qc/'annotation_scope_audit.csv',audit_rows)
    old=json.loads((ROOT/'outputs/s_cerevisiae_20261030/analysis.json').read_text())['rows']
    previous=[r for r in old if r['included_in_primary_plots']]
    level_summary=[]
    for level in LEVELS:
        group=[r for r in audit_rows if r['assembly_level']==level]
        before=[r for r in previous if r['assembly_level']==level]
        level_summary.append(dict(assembly_level=level,n_before_dedup=len(before),n_after_dedup=len(group),
               protein_median_before_dedup=float(np.median([r['number_of_proteins'] for r in before])),
               protein_median_after_dedup=float(np.median([r['number_of_proteins'] for r in group])),
               gc_all_median=float(np.median([r['gc_all_percent'] for r in group])),
               gc_organelle_excluded_median=float(np.median([r['gc_organelle_excluded_percent'] for r in group])),
               organelle_included_assemblies=sum(r['explicit_organelle_bp']>0 for r in group),
               organelle_only_annotation_assemblies=sum(r['annotation_scope_mismatch'] for r in group)))
    csv_write(OUT/'assembly_level_scope_summary.csv',level_summary)
    strata=defaultdict(list)
    for r in audit_rows:strata[(r['annotation_provider'],r['sequencing_group'])].append(r)
    stratum_rows=[]
    for (provider,technology),group in strata.items():
        for level in sorted({r['assembly_level'] for r in group}):
            subset=[r for r in group if r['assembly_level']==level]
            stratum_rows.append(dict(annotation_provider=provider,sequencing_group=technology,assembly_level=level,n=len(subset),
                                    levels_in_stratum=len({r['assembly_level'] for r in group}),
                                    median_proteins=float(np.median([r['number_of_proteins'] for r in subset]))))
    csv_write(OUT/'provider_sequencing_strata.csv',stratum_rows)
    write_json(OUT/'audit_summary.json',dict(assemblies=len(audit_rows),level_summary=level_summary,
               scope_mismatch_accessions=[r['assembly_accession'] for r in audit_rows if r['annotation_scope_mismatch']],
               assemblies_with_explicit_organelles=sum(r['explicit_organelle_bp']>0 for r in audit_rows),
               provider_sequencing_strata=len(strata),
               multi_level_strata=sum(len({r['assembly_level'] for r in group})>1 for group in strata.values()),
               assumptions='Unlabeled contigs retained as nuclear_or_unclassified; only explicitly labeled organelle/plasmid records excluded',
               causal_effect_of_assembly_level_identified=False))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(13,5.8));rng=np.random.default_rng(4932)
    for ax,key,title in [(axes[0],'gc_all_percent','GC% across all supplied FASTA records'),
                         (axes[1],'gc_organelle_excluded_percent','GC% after excluding labeled organelles/plasmids')]:
        values=[[r[key] for r in audit_rows if r['assembly_level']==level] for level in LEVELS]
        boxes=ax.boxplot(values,patch_artist=True,showfliers=False,medianprops={'color':'#182634','linewidth':1.5})
        for j,(v,box,color) in enumerate(zip(values,boxes['boxes'],COLORS),1):
            box.set_facecolor(color);box.set_alpha(.35)
            ax.scatter(j+rng.uniform(-.15,.15,len(v)),v,s=13,color=color,alpha=.45,edgecolors='none')
        ax.set_xticks(range(1,5),[l.replace('Complete Genome','Complete\nGenome')+'\n(n='+str(len(v))+')' for l,v in zip(LEVELS,values)])
        ax.set_title(title,loc='left',fontsize=12,fontweight='bold',pad=12)
        ax.set_ylabel('GC (%)');ax.set_xlabel('Assembly level');ax.grid(axis='y',color='#E1E5E9');ax.set_axisbelow(True)
    limits=[ax.get_ylim() for ax in axes];lo=min(x[0] for x in limits);hi=max(x[1] for x in limits)
    for ax in axes:ax.set_ylim(lo,hi)
    fig.suptitle('A within-data scope check, not a biological-loss test',x=.075,ha='left',fontsize=16,fontweight='bold',y=.985)
    fig.text(.075,.92,'Same 279 strain/sample representatives in both panels. No outlier exclusion.',fontsize=10)
    fig.text(.075,.025,'GC = 100 × (G+C)/(A+C+G+T). Only explicitly labeled non-nuclear records are removed.\n'
             'Differences by assembly level remain confounded with submitter, sequencing and annotation.',fontsize=9,linespacing=1.5)
    fig.subplots_adjust(left=.075,right=.98,bottom=.26,top=.80,wspace=.22)
    fig.savefig(OUT/'gc_scope_comparison.png',dpi=190,facecolor='white')
    fig.savefig(OUT/'gc_scope_comparison.pdf',facecolor='white');plt.close(fig)
    print(json.dumps(level_summary,indent=2),flush=True)


if __name__=='__main__':main()
