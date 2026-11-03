#!/usr/bin/env python3
"""Find genuine NCBI FCS-positive E. coli controls; keep labels out of inference."""
import csv,gzip,json,pathlib,subprocess,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1];D=ROOT/'data/contamination_validation_20261102'
def main():
 candidates=[]
 with gzip.open(D/'fcs_summary_genbank.txt.gz','rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   fraction=float(r['fraction_contaminated']);length=int(r['total_contaminant_length'])
   if r['asserted_division'] in ['prok:g-proteobacteria','prok:proteobacteria'] and .05<=fraction<=.25 and 4500000<=length/fraction<=6500000:candidates.append(r)
 candidates.sort(key=lambda r:float(r['fraction_contaminated']),reverse=True)
 # Examine a deterministic large batch, not an assumption of taxonomic identity.
 selected=candidates[:300];reports=[]
 for index in range(0,len(selected),40):
  batch=selected[index:index+40];accessions=','.join(r['#assembly_accession'] for r in batch);p=D/f'gamma_batch_{index//40:02d}.json';url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+accessions+'/dataset_report?page_size=1000'
  if not p.exists():subprocess.run(['curl','-fLsS','--retry','2','--max-time','120',url,'-o',str(p)],check=True)
  reports.extend(json.loads(p.read_text()).get('reports',[]));print('metadata',len(reports),flush=True)
 lookup={r['#assembly_accession']:r for r in selected};ecoli=[]
 for r in reports:
  if r.get('organism',{}).get('organism_name','').startswith('Escherichia coli') and r['accession'] in lookup:ecoli.append(dict(report=r,fcs=lookup[r['accession']]))
 (D/'confirmed_ecoli_fcs_positive.json').write_text(json.dumps(dict(screened_fcs_candidates=len(candidates),metadata_examined=len(reports),ecoli_controls=ecoli,selection='Top 300 FCS-contaminated Proteobacteria by contaminated fraction; species verified by NCBI metadata; not a representative sample',labels_for_validation_only=True),indent=2)+'\n');print('confirmed E. coli',len(ecoli),flush=True)
if __name__=='__main__':main()
