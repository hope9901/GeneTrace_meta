#!/usr/bin/env python3
"""Targeted DNA rescue and flanking homolog evidence for missing protein clusters."""
import collections,csv,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.io import read_fasta,read_gff,write_fasta,sha256
from genetrace.engine import find_binary,parse_hits,FIELDS
from genetrace.cohort_comparison import anchor_support
GROUP=ROOT/'outputs/ecoli100_gene_groups_20261102';OUT=ROOT/'outputs/ecoli100_candidate_audit_20261102';RUN=ROOT/'runs/ecoli100_candidate_audit_20261102'
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True);metadata={r['assembly_accession']:r for r in json.loads((ROOT/'data/bacteria/cohort.json').read_text())};missing=list(csv.DictReader((GROUP/'missing_candidates.csv').open()));runs=list(csv.DictReader((ROOT/'outputs/ecoli100_reads_20261102/runs.csv').open()));small={}
 for r in runs:
  if not r['fastq_bytes']:continue
  size=sum(int(x) for x in r['fastq_bytes'].split(';'))
  if size<500000000 and (r['assembly_accession'] not in small or size<small[r['assembly_accession']][0]):small[r['assembly_accession']]=(size,r)
 ranked=sorted([r for r in missing if r['assembly_accession'] in small],key=lambda r:(small[r['assembly_accession']][0],-int(r['group_prevalence']),r['reference_locus']));selected=[];perasm=collections.Counter()
 for r in ranked:
  if perasm[r['assembly_accession']]>=3:continue
  selected.append(r);perasm[r['assembly_accession']]+=1
  if len(selected)==9:break
 ref='GCF_000005845.2';proteins=read_fasta(ROOT/metadata[ref]['source_directory']/'protein.faa.gz');gm=list(csv.DictReader((ROOT/'data/bacteria/gene_maps'/f'{ref}.tsv').open(),delimiter='\t'));ordered=sorted(gm,key=lambda r:(r['seqid'],int(r['start'])));results=[]
 for i,c in enumerate(selected):
  acc=c['assembly_accession'];g=next(g for g in ordered if g['locus_tag']==c['reference_locus']);idx=ordered.index(g)
  if idx==0 or idx==len(ordered)-1:continue
  flanks={'left':ordered[idx-1],'right':ordered[idx+1]};d=RUN/(acc+'_'+g['locus_tag']);d.mkdir(exist_ok=True);query=d/'query.faa';write_fasta(query,{'center':proteins[g['ncbi_protein_id']],**{k:proteins[v['ncbi_protein_id']] for k,v in flanks.items()}});genome=read_fasta(ROOT/metadata[acc]['source_directory']/'genome.fna.gz');subject=d/'genome.fna';write_fasta(subject,genome);cmd=[find_binary('tblastn'),'-query',str(query),'-subject',str(subject),'-seg','no','-comp_based_stats','0','-num_threads','1','-evalue','1e-5','-outfmt','6 '+FIELDS,'-max_target_seqs','100'];p=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=180);(d/'hits.tsv').write_text(p.stdout);hits=parse_hits(p.stdout);centers=sorted([h for h in hits if h['qseqid']=='center'],key=lambda h:h['bitscore'],reverse=True);h=centers[0] if centers else None
  result=dict(**c,run_candidate=small[acc][1],command=cmd,query_sha256=sha256(query),source_genome_sha256=sha256(ROOT/metadata[acc]['source_directory']/'genome.fna.gz'),best_genome_hit=h,classification='no_hit_unresolved',assembly_read_use_verified=False,raw_read_validation_complete=False)
  if h:
   support={role:any(hit['qseqid']==role and anchor_support(hit,h,g,flank,role) for hit in hits) for role,flank in flanks.items()};features=read_gff(ROOT/metadata[acc]['source_directory']/'annotation.gff.gz');cds=[f for f in features if f['type']=='CDS' and f['seqid']==h['sseqid'] and f['strand']==h['strand'] and f['start']<=h['end'] and f['end']>=h['start']];result.update(flank_support=support,both_flanks=all(support.values()),overlapping_cds=[dict(start=f['start'],end=f['end'],attrs=f['attrs']) for f in cds],classification='DNA_homolog_present_missing_protein_cluster_unresolved' if h['query_coverage']>=.9 and h['pident']>=90 else 'partial_or_divergent_DNA_homolog_unresolved')
  dump(d/'audit.json',result);results.append(result);dump(OUT/'results.json',results);dump(OUT/'status.json',dict(stage='auditing',completed=i+1,total=len(selected)));print(i+1,result['classification'],flush=True)
 dump(OUT/'status.json',dict(stage='complete',completed=len(results),total=len(selected),raw_read_validation_complete=False))
if __name__=='__main__':main()
