"""Within-BioSample assembly-representation contrasts from audited metadata."""
import csv
import itertools
import json
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/genetrace-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/genetrace-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/loss_artifact_tests_20261030'


def main():
    all_rows=json.loads((ROOT/'outputs/s_cerevisiae_methods_20261030/methods_analysis.json').read_text())['rows']
    samples=defaultdict(list)
    for row in all_rows:
        if row.get('biosample'):samples[row['biosample']].append(row)
    pairs=[]
    for sample,group in sorted(samples.items()):
        by_level=defaultdict(list)
        for row in group:by_level[row['assembly_level']].append(row)
        for left,right in itertools.combinations(sorted(by_level),2):
            # Only one contrast per level pair/sample; duplicate deposits must not
            # double the apparent sample size. Prefer later release deterministically.
            a=max(by_level[left],key=lambda r:(r['release_date'],r['assembly_accession']))
            b=max(by_level[right],key=lambda r:(r['release_date'],r['assembly_accession']))
            pairs.append(dict(biosample=sample,level_a=left,level_b=right,
                              assembly_a=a['assembly_accession'],assembly_b=b['assembly_accession'],
                              proteins_a=a['number_of_proteins'],proteins_b=b['number_of_proteins'],
                              delta_proteins_b_minus_a=b['number_of_proteins']-a['number_of_proteins'],
                              size_mb_a=a['genome_size_mb'],size_mb_b=b['genome_size_mb'],
                              delta_size_mb_b_minus_a=b['genome_size_mb']-a['genome_size_mb'],
                              gc_percent_a=a['gc_percent'],gc_percent_b=b['gc_percent'],
                              same_reported_sequencing=a['sequencing_group']==b['sequencing_group'],
                              same_annotation_provider=a['annotation_provider']==b['annotation_provider'],
                              sequencing_group_a=a['sequencing_group'],sequencing_group_b=b['sequencing_group'],
                              annotation_provider_a=a['annotation_provider'],annotation_provider_b=b['annotation_provider']))
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'same_biosample_cross_level_pairs.csv').open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(pairs[0]));writer.writeheader();writer.writerows(pairs)
    focus=[r for r in pairs if r['level_a']=='Contig' and r['level_b']=='Scaffold']
    summary=dict(unique_biosamples=len({r['biosample'] for r in pairs}),level_pair_contrasts=len(pairs),
                 contig_scaffold_pairs=len(focus),
                 scaffold_more_proteins=sum(r['delta_proteins_b_minus_a']>0 for r in focus),
                 equal_proteins=sum(r['delta_proteins_b_minus_a']==0 for r in focus),
                 scaffold_fewer_proteins=sum(r['delta_proteins_b_minus_a']<0 for r in focus),
                 median_protein_difference=float(np.median([r['delta_proteins_b_minus_a'] for r in focus])),
                 min_protein_difference=min(r['delta_proteins_b_minus_a'] for r in focus),
                 max_protein_difference=max(r['delta_proteins_b_minus_a'] for r in focus),
                 median_span_difference_mb=float(np.median([r['delta_size_mb_b_minus_a'] for r in focus])),
                 same_sequencing_pairs=sum(r['same_reported_sequencing'] for r in focus),
                 same_provider_pairs=sum(r['same_annotation_provider'] for r in focus),
                 biological_sample_identity='Same NCBI BioSample, not independently sequence-verified clone identity',
                 conclusion='Observed protein counts vary between assembly representations of the same registered biological sample.',
                 limitations=['Observational paired metadata; not randomized assembly methods.',
                              'Differences may reflect completeness, phasing/haplotype representation or annotation; no single mechanism proven.',
                              'Not a claim that Scaffold level is always better, or that all group differences are technical.'])
    (OUT/'same_biosample_pair_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.5,5.6))
    for ax,key,title in [(axes[0],'proteins','Submitted protein FASTA records'),(axes[1],'size_mb','Genome FASTA span (Mb)')]:
        for pair in focus:
            ax.plot([0,1],[pair[key+'_a'],pair[key+'_b']],color='#B8C2CC',linewidth=1,alpha=.75,zorder=1)
        ax.scatter(np.zeros(len(focus)),[r[key+'_a'] for r in focus],color='#8064A2',s=26,zorder=2,label='Contig')
        ax.scatter(np.ones(len(focus)),[r[key+'_b'] for r in focus],color='#B6873D',s=26,zorder=2,label='Scaffold')
        ax.set_xticks([0,1],['Contig','Scaffold']);ax.set_xlim(-.15,1.15)
        ax.set_ylabel(title);ax.set_title(title,loc='left',fontweight='bold',fontsize=12)
        ax.grid(axis='y',color='#E1E5E9');ax.set_axisbelow(True)
    fig.suptitle('Same registered sample, different assembly representation',x=.09,ha='left',fontsize=16,fontweight='bold',y=.98)
    fig.text(.09,.91,str(len(focus))+' paired BioSamples; duplicate deposits collapsed within each sample/level',fontsize=10)
    fig.text(.09,.025,'Lines connect the same BioSample. These observations do not establish an independent causal effect of assembly level.\n'
             'All pairs have the same reported sequencing group and annotation provider; representation/annotation can still differ.',fontsize=9,linespacing=1.5)
    fig.subplots_adjust(left=.09,right=.98,top=.82,bottom=.19,wspace=.29)
    fig.savefig(OUT/'same_sample_assembly_pairs.png',dpi=190,facecolor='white')
    fig.savefig(OUT/'same_sample_assembly_pairs.pdf',facecolor='white');plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
