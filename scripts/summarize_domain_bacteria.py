#!/usr/bin/env python3
"""Adjudicate cached genomic domain hits and describe completed domain scans."""
import collections,csv,json,pathlib,sys,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/genetrace-matplotlib');os.environ.setdefault('XDG_CACHE_HOME','/tmp/genetrace-cache')
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_gff,write_json
from genetrace.domains import domain_coordinates,genomic_assessment
OUT=ROOT/'outputs/domain_bacteria_20261101';RUN=ROOT/'runs/domain_bacteria_20261101'
def main():
 status=json.loads((OUT/'status.json').read_text())
 if status['stage']!='complete':raise RuntimeError('Wait for the full-cohort scan before final summary')
 summary=json.loads((OUT/'summary.json').read_text());sources=json.loads((ROOT/'data/bacteria/cohort.json').read_text());rows=list(csv.DictReader((OUT/'domain_counts.csv').open()));adjudicated=[]
 for source in sources:
  acc=source['assembly_accession'];d=RUN/acc;features=read_gff(ROOT/source['source_directory']/'annotation.gff.gz');mapping=json.loads((d/'orf_mapping.json').read_text());hits=json.loads((d/'genome_domain/hits.json').read_text());ph=json.loads((d/'proteome/hits.json').read_text())
  assessments=[genomic_assessment(h,domain_coordinates(mapping[h['protein_id']],h['aa_start'],h['aa_end']),features,ph) for h in hits]
  write_json(d/'genomic_domain_assessments.json',assessments);adjudicated.extend(dict(assembly_accession=acc,strain=source['strain'],**r) for r in assessments)
 with (OUT/'genomic_domain_candidates.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(adjudicated[0]));w.writeheader();w.writerows(adjudicated)
 counts=collections.Counter((r['strain'],r['classification']) for r in adjudicated);write_json(OUT/'genomic_adjudication_summary.json',[dict(strain=k[0],classification=k[1],domain_instances=n) for k,n in sorted(counts.items())])
 # Profile counts are distinct proteins, not domain instances or orthogroups.
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;import numpy as np
 names=[s['strain'] for s in sources];profiles=[m['accession'].split('.')[0] for m in summary['profiles']];data=np.array([[next(int(r['proteins_with_domain']) for r in rows if r['strain']==name and r['profile']==p) for p in profiles] for name in names])
 fig,ax=plt.subplots(figsize=(10,4),layout='constrained');im=ax.imshow(data,cmap='YlGn',aspect='auto',vmin=0);ax.set_xticks(range(len(profiles)),[p+'\n'+m['name'] for p,m in zip(profiles,summary['profiles'])],fontsize=9);ax.set_yticks(range(len(names)),names)
 for i in range(len(names)):
  for j in range(len(profiles)):ax.text(j,i,str(data[i,j]),ha='center',va='center',color='white' if data[i,j]>data.max()*.65 else '#16392e')
 ax.set_title('E. coli: proteins carrying each selected Pfam profile');fig.colorbar(im,ax=ax,label='Distinct annotated proteins');fig.savefig(OUT/'bacterial_domain_counts.png',dpi=180);plt.close(fig)
 text=f'''# Domain profiles and bacterial transfer test — 2026-11-01

## Scope

All annotated proteins from {summary['yeast_assemblies']} representative S. cerevisiae assemblies and four E. coli laboratory model strains were scanned with six selected Pfam HMM profiles: PF00005 (ABC_tran), PF00072 (Response_reg), PF00512 (HisKA), PF02518 (HATPase_c), PF00069 (Pkinase), PF00076 (RRM_1). This is a selected-panel test, not a whole-Pfam annotation. Current profile versions and SHA256 are in `data/domain_models/manifest.json`.

Bacterial genome, GFF, CDS and protein files were downloaded from NCBI and checked against the NCBI MD5 manifests. Strains are MG1655, W3110, REL606 and BL21(DE3), all complete assemblies. Their source reports and exact accessions are in `data/bacteria/cohort.json`.

## Method

HMMER 3.4 `hmmsearch --cpu 1 --cut_ga` searches each complete proteome. Both sequence and domain scores must meet profile-specific Pfam gathering thresholds. HMM coverage >=70% labels a near-complete domain; shorter GA-passing hits remain as partial evidence. Per-profile count tables distinguish proteins carrying the domain from domain instances; repeated domains in one protein are not extra genes. Architecture strings describe only these six profiles and are not full protein architectures.

For bacteria, all six reading frames are translated, split at stop codons and ORFs of at least 30 amino acids are searched with the same profiles. Internal codon translations follow table 11; these ORFs are not gene predictions and alternative start codons, complete CDS structure and function are not inferred. Hits are projected to 1-based genomic coordinates and compared with same-strand CDS. Pseudogene, frameshift and interrupted model records are inspected before suggesting annotation omission. Whole-genome profile searches have not been run on intron-containing yeast genomes.

## Findings and controls

MG1655 has 77 proteins carrying PF00005, with 94 domain instances. This illustrates why domain instance count and gene count must remain separate. A genomic PF00005 hit outside protein-coding annotation overlaps the annotated pseudogene **ycjV**: classify it as `pseudogene_or_disrupted_model_unresolved`, not an omitted functional gene.

The synthetic positive control removes MG1655 protein NP_416685.1 and its structural annotation from copied input files. Its PF00005 hit disappears from the protein search while the unchanged genome still contains the domain, recovering `unannotated_genomic_domain_candidate`. Original source files remain unchanged. This establishes recovery of one known annotation omission, not sensitivity to natural errors generally.

The existing whole-protein algorithm was also run on rpoB in all four strains: a full-length DNA hit at 100% amino-acid identity was detected in each, along with protein/structural annotation evidence. This is a conserved-gene positive control, not proof that every bacterial gene or artifact type is handled correctly.

JXXY16.1 has zero selected-profile hits in the supplied 16-protein organelle-only annotation. This does not establish nuclear domain or gene loss. Pkinase/RRM_1 profile absence in E. coli does not establish absence of every kinase or RNA-binding protein family; these profiles are specific families.

## Outputs

- `domain_counts.csv`: per-assembly, per-profile distinct protein / domain / partial counts.
- `domain_members.csv`: domain coordinates, scores, HMM coverage and protein IDs.
- `panel_architectures.csv`: order and repeats of the selected profiles on proteins.
- `genomic_domain_candidates.csv`: bacterial DNA coordinates and structural-annotation adjudication.
- `controls.json`: synthetic annotation-omission control outcome.
- `bacterial_domain_counts.png`: bacterial profile count comparison.

Raw HMMER tables, command/version/hash manifests, six-frame maps and control inputs are under `runs/domain_bacteria_20261101/`. No bacterial raw reads were downloaded; no sequencing error, biological loss or functional gene status is confirmed by these tests. Natural differences between the four assemblies require locus-specific read/synteny evidence and provenance checks.

References: [InterPro API](https://interpro-documentation.readthedocs.io/en/latest/api.html), [HMMER documentation](https://hmmer-web-docs.readthedocs.io/en/latest/searches.html), [NCBI MG1655 assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/).
'''
 (OUT/'report.md').write_text(text);(ROOT/'docs/DOMAIN_AND_BACTERIA.md').write_text(text);print('Adjudicated',len(adjudicated),'genomic domain instances; final report written')
if __name__=='__main__':main()
