#!/usr/bin/env python3
"""Locate a gene across the entire same-species cohort; compare callable neighborhoods."""
import argparse,collections,concurrent.futures,csv,datetime,json,pathlib,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.agent import cohort,source_dir,reference_for,resolve_gene,genes,load
from genetrace.engine import find_binary,parse_hits,FIELDS
from genetrace.io import read_fasta,read_gff,write_fasta,write_json,sha256
from genetrace.cohort_comparison import platform_group,anchor_support,distinct_loci,parse_cs,repeat_allele,repeat_bounds,revcomp,summarize_marker,variable_reference_positions

def table(path,rows):
 if not rows:return
 temp=path.with_suffix('.tmp')
 with temp.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 temp.replace(path)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--assembly',default='GCA_903819125.2');ap.add_argument('--gene',default='YJU3');ap.add_argument('--workers',type=int,default=2);a=ap.parse_args()
 focal=a.assembly;ref=reference_for(focal);gene=resolve_gene(a.gene,ref);key=gene['locus_tag'] or gene['gene_name'];out=ROOT/'outputs/cohort_comparison_20261102'/key;run=ROOT/'runs/cohort_comparison_20261102'/key;out.mkdir(parents=True,exist_ok=True);run.mkdir(parents=True,exist_ok=True)
 selected=[r for acc,r in cohort().items() if reference_for(acc)==ref];ordered=sorted([g for g in genes(ref) if g['seqid']==gene['seqid']],key=lambda g:int(g['start']));index=next(i for i,g in enumerate(ordered) if g['gene_id']==gene['gene_id'])
 if index==0 or index==len(ordered)-1:raise ValueError('Two reference flanks are required')
 flanks={'left':ordered[index-1],'right':ordered[index+1]};proteins=read_fasta(source_dir(ref)/'protein.faa.gz');query=run/'queries.faa';write_fasta(query,{'center':proteins[gene['ncbi_protein_id']],**{k:proteins[v['ncbi_protein_id']] for k,v in flanks.items()}});query_hash=sha256(query);binary=find_binary('tblastn')
 started=time.monotonic();results=[]
 def investigate(row):
  acc=row['assembly_accession'];d=run/acc;d.mkdir(exist_ok=True);sig=dict(query_sha256=query_hash,genome_sha256=row.get('genome_fasta_sha256') or sha256(source_dir(acc)/'genome.fna.gz'))
  if (d/'result.json').exists():
   cached=load(d/'result.json')
   if cached.get('input_signature')==sig:return cached
  sequences=read_fasta(source_dir(acc)/'genome.fna.gz');genome=d/'genome.fna';write_fasta(genome,sequences)
  cmd=[binary,'-query',str(query),'-subject',str(genome),'-seg','no','-comp_based_stats','0','-num_threads','1','-evalue','1e-5','-outfmt','6 '+FIELDS,'-max_target_seqs','100']
  p=subprocess.run(cmd,capture_output=True,text=True,timeout=240);(d/'hits.tsv').write_text(p.stdout);(d/'stderr.txt').write_text(p.stderr)
  base=dict(assembly_accession=acc,strain=row['resolved_strain'],sequencing_group=platform_group(row),sequencing_metadata_raw=row.get('sequencing_group','N/A'),bioproject=row.get('bioproject','N/A'),submitter=row.get('submitter','N/A'),assembly_level=row['assembly_level'],input_signature=sig,command=cmd,biological_loss_inferred=False,error_inferred=False)
  if p.returncode:base.update(status='search_failed',strong_loci=0,both_flanks=False,eligible_for_frequency=False,error=p.stderr[-300:]);write_json(d/'result.json',base);return base
  hits=parse_hits(p.stdout);centers=distinct_loci([h for h in hits if h['qseqid']=='center' and h['pident']>=90 and h['query_coverage']>=.9 and h['evalue']<=1e-20]);base['strong_loci']=len(centers)
  if not centers:base.update(status='no_qualifying_hit_unresolved',both_flanks=False,eligible_for_frequency=False);write_json(d/'result.json',base);return base
  center=centers[0];supported={role:any(h['qseqid']==role and anchor_support(h,center,gene,flank,role) for h in hits) for role,flank in flanks.items()}
  both=all(supported.values());features=read_gff(source_dir(acc)/'annotation.gff.gz');overlap=[f for f in features if f['seqid']==center['sseqid'] and f['start']<=center['end'] and f['end']>=center['start'] and f['type'] in ['gene','pseudogene','CDS']]
  cds=[f for f in overlap if f['type']=='CDS' and f['strand']==center['strand']];sequence=sequences[center['sseqid']];lo=max(1,center['start']-500);hi=min(len(sequence),center['end']+500);fragment=sequence[lo-1:hi];fragment=fragment if center['strand']=='+' else revcomp(fragment)
  write_fasta(d/'locus.fna',{acc:fragment});base.update(status='single_synteny_supported_locus' if len(centers)==1 and both else 'multiple_loci_ambiguous' if len(centers)>1 else 'flank_support_incomplete',best_hit={k:center[k] for k in ['sseqid','start','end','strand','pident','query_coverage','bitscore']},both_flanks=both,flank_support=supported,eligible_for_frequency=len(centers)==1 and both,overlapping_cds=len(cds),recorded_disruption=any(f['type']=='pseudogene' or f['attrs'].get('pseudo')=='true' or any(w in str(f['attrs']).lower() for w in ['frameshift','interrupted','pseudogene']) for f in overlap),locus_start=lo,locus_end=hi,locus_strand=center['strand'],locus_sha256=sha256(d/'locus.fna'))
  write_json(d/'result.json',base);return base
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  futures={pool.submit(investigate,r):r for r in selected}
  for future in concurrent.futures.as_completed(futures):
   row=futures[future]
   try:result=future.result()
   except Exception as e:result=dict(assembly_accession=row['assembly_accession'],strain=row['resolved_strain'],status='search_failed',error=str(e),eligible_for_frequency=False,both_flanks=False,strong_loci=0)
   results.append(result);write_json(out/'status.json',dict(stage='locating_loci',focal_assembly=focal,gene=gene['gene_name'],reference_gene_id=key,completed=len(results),total=len(selected),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
   if len(results)%20==0:print('located',len(results),'/',len(selected),flush=True)
 results.sort(key=lambda r:r['assembly_accession']);write_json(out/'assemblies.json',results)
 columns=['assembly_accession','strain','status','strong_loci','both_flanks','eligible_for_frequency','sequencing_group','bioproject','submitter','assembly_level','overlapping_cds','recorded_disruption']
 table(out/'assemblies.csv',[{k:r.get(k,'N/A') for k in columns} for r in results])
 target=next(r for r in results if r['assembly_accession']==focal)
 if not target.get('eligible_for_frequency'):raise RuntimeError('Focal locus is not unambiguous and synteny supported')
 refseq=next(iter(read_fasta(run/focal/'locus.fna').values()));write_fasta(run/'focal_locus.fna',{'focal':refseq})
 subjects={r['assembly_accession']:next(iter(read_fasta(run/r['assembly_accession']/'locus.fna').values())) for r in results if r.get('eligible_for_frequency')};write_fasta(run/'candidate_loci.fna',subjects)
 mm=ROOT/'tools/readqc-env/bin/minimap2';cmd=[str(mm),'-cx','asm10','--cs=long','--secondary=no','-t','1',str(run/'focal_locus.fna'),str(run/'candidate_loci.fna')];p=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=120);(run/'locus_alignment.paf').write_text(p.stdout);(run/'minimap2.stderr').write_text(p.stderr)
 # Initial markers correspond to actual HN1 indel observations; novel targets can use a future marker selection step.
 marker_positions=[287652,287335,287295,287486] if focal=='GCA_903819125.2' and gene['gene_name']=='YJU3' else []
 comparisons=[];best={}
 for line in p.stdout.splitlines():
  f=line.split('\t');q=f[0];length=int(f[1]);strand=f[4];score=int(f[10]);cs=next(x[5:] for x in f[12:] if x.startswith('cs:Z:'))
  if q in best and best[q][0]>=score:continue
  oriented=subjects[q] if strand=='+' else revcomp(subjects[q]);qs=int(f[2]) if strand=='+' else length-int(f[3]);mapping,ins=parse_cs(cs,int(f[7]),qs,oriented);best[q]=(score,mapping,ins,float(f[9])/int(f[10]))
 # For targets without curated read markers, scan the entire located coding span.
 # Collapse changed positions in the same reference homopolymer into one block.
 marker_selection='previous_HN1_read_observations' if marker_positions else 'all_variable_reference_repeat_blocks_in_located_coding_span'
 if not marker_positions:
  coding_positions=[pos for pos in range(len(refseq)) if target['best_hit']['start']<=(target['locus_start']+pos if target['locus_strand']=='+' else target['locus_end']-pos)<=target['best_hit']['end']]
  positions=variable_reference_positions(refseq,[(v[1],v[2]) for acc,v in best.items() if acc!=focal],min(coding_positions),max(coding_positions)+1)
  marker_positions=[target['locus_start']+pos if target['locus_strand']=='+' else target['locus_end']-pos for pos in positions]
 markers=[];focal_technology=platform_group(cohort()[focal])
 for original_pos in marker_positions:
  pos=original_pos-target['locus_start'] if target['locus_strand']=='+' else target['locus_end']-original_pos
  lo,hi,base=repeat_bounds(refseq,pos);rows=[]
  for r in results:
   if r['assembly_accession']==focal:continue
   acc=r['assembly_accession'];allele=None
   if r.get('eligible_for_frequency') and acc in best:allele=repeat_allele(refseq,pos,best[acc][1],best[acc][2])
   info=dict(assembly_accession=acc,strain=r['strain'],position=original_pos,allele=allele['allele'] if allele else None,repeat_length=allele['repeat_length'] if allele else None,sequencing_group=r.get('sequencing_group','N/A'),bioproject=r.get('bioproject','N/A'),assembly_level=r.get('assembly_level','N/A'),locus_alignment_identity=best[acc][3] if acc in best else None);rows.append(info);comparisons.append(info)
  stats=summarize_marker(rows,refseq[lo:hi+1]);groups=[]
  for field in ['sequencing_group','bioproject','assembly_level']:
   for value in sorted({r[field] for r in rows}):groups.append(dict(field=field,value=value,**summarize_marker([r for r in rows if r[field]==value],refseq[lo:hi+1])))
  markers.append(dict(focal_reference_position=original_pos,repeat_base=base,focal_repeat_length=hi-lo+1,population=stats,same_technology=summarize_marker([r for r in rows if focal_technology!='N/A' and r['sequencing_group']==focal_technology],refseq[lo:hi+1]),different_known_technology=summarize_marker([r for r in rows if r['sequencing_group'] not in ['N/A',focal_technology]],refseq[lo:hi+1]),strata=groups))
 table(out/'marker_alleles.csv',comparisons)
 summary=dict(focal_assembly=focal,gene=gene['gene_name'],reference_gene_id=key,total_assemblies=len(results),other_assemblies=len(results)-1,unambiguous_synteny_supported=sum(r.get('eligible_for_frequency',False) for r in results),locus_status_counts=dict(collections.Counter(r['status'] for r in results)),markers=markers,marker_selection=marker_selection,focal_technology=focal_technology,elapsed_seconds=time.monotonic()-started,phylogenetic_adjustment='not_available; locus similarity is not genome phylogeny',raw_reads_examined='HN1 only' if focal=='GCA_903819125.2' else 'none for this cohort comparison',error_inferred=False,biological_loss_inferred=False,limitations=['Rare or private variants can be real; majority allele is not truth.','Strains and project batches are correlated; counts are not independent trials.','Technology labels do not establish run-to-assembly provenance or isolate instrument effects.','Markers use alignment-dependent repeat-length alleles, not normalized variant calls.','Ambiguous loci, incomplete synteny and uncallable bases are excluded from denominators, not counted as loss.','Orthogroup/phylogenetic confirmation remains pending; intron-rich genes may fail the single-HSP locator.'])
 write_json(out/'summary.json',summary);write_json(out/'alignment_manifest.json',dict(command=cmd,minimap2_version=subprocess.check_output([str(mm),'--version'],text=True).strip(),tblastn_version=subprocess.check_output([binary,'-version'],text=True).splitlines()[0],query_sha256=query_hash,focal_locus_sha256=sha256(run/'focal_locus.fna'),loci_sha256=sha256(run/'candidate_loci.fna')));write_json(out/'status.json',dict(stage='complete',completed=len(results),total=len(selected),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));print(json.dumps({k:summary[k] for k in ['total_assemblies','unambiguous_synteny_supported','locus_status_counts','elapsed_seconds']},indent=2),flush=True)
if __name__=='__main__':main()
