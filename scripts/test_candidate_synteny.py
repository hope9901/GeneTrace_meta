"""Check two flanking coding loci for selected diagnostic homolog candidates."""
import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import read_fasta,write_fasta,write_json
from genetrace.engine import EvidenceStore,find_binary,parse_hits,FIELDS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/loss_artifact_tests_20261030'
RUN=ROOT/'runs/loss_artifact_tests_20261030'
REFERENCE='GCA_000146045.2'


def compatible(hit,center,reference_center,reference_flank,role,max_distance=25000):
    if hit['sseqid']!=center['seqid'] or hit['pident']<70 or hit['query_coverage']<.5 or hit['evalue']>1e-20:return False
    same=reference_center['strand']==center['strand']
    expected=reference_flank['strand'] if same else ('-' if reference_flank['strand']=='+' else '+')
    if hit['strand']!=expected:return False
    # Each anchor must be a separate nonoverlapping locus on the expected side.
    before=(role=='left')==same
    if before:return 0<center['start']-hit['end']<=max_distance
    return 0<hit['start']-center['end']<=max_distance


def main():
    rows=json.loads((OUT/'adjudicated_cases.json').read_text())['loci']
    ref=list(csv.DictReader((ROOT/'data/orthology/gene_maps'/(REFERENCE+'.tsv')).open(),delimiter='\t'))
    genes={r['gene_id']:r for r in ref};ordered=defaultdict(list)
    for r in ref:ordered[r['seqid']].append(r)
    neighbors={}
    for seqid,chromosome in ordered.items():
        chromosome.sort(key=lambda r:int(r['start']))
        for i,gene in enumerate(chromosome):
            if 0<i<len(chromosome)-1:neighbors[gene['gene_id']]=(chromosome[i-1],chromosome[i+1])
    chosen=[]
    j=[r for r in rows if r['target_assembly']=='GCA_903819185.2' and r['genome_homolog_found']]
    seen=set()
    for row in j:
        seqid=genes[row['reference_gene_id']]['seqid']
        if seqid not in seen and len(seen)<3:chosen.append(row);seen.add(seqid)
    chosen.extend(r for r in rows if r['target_assembly']=='GCA_903819185.2' and r['partial_or_divergent_genome_evidence'])
    chosen.extend(r for r in rows if r['target_assembly']=='GCA_903819125.2' and r['genome_homolog_found'] and not r['protein_homolog_found'])
    for accession in ('GCA_000190255.1','GCA_000326005.1'):
        candidates=[r for r in rows if r['target_assembly']==accession and r['genome_homolog_found'] and not r['protein_homolog_found']]
        # One unannotated and one overlapping-model example where available.
        for overlap in (False,True):
            eligible=[r for r in candidates if bool(r['overlapping_cds'])==overlap]
            if eligible:chosen.append(eligible[0])
    by_target=defaultdict(list)
    for row in chosen:
        if row['reference_gene_id'] in neighbors:by_target[row['target_assembly']].append(row)
    proteins=read_fasta(ROOT/'data/orthology/source'/REFERENCE/'protein.faa.gz')
    results=[]
    for accession,cases in by_target.items():
        directory=RUN/'synteny'/accession
        if (directory/'report.json').exists():
            cached=json.loads((directory/'report.json').read_text())
            current={r['reference_gene_id']:r for r in rows if r['target_assembly']==accession}
            for candidate in cached:
                if candidate['reference_gene_id'] in current:
                    candidate['source_classification']=current[candidate['reference_gene_id']]['classification']
            results.extend(cached);continue
        store=EvidenceStore(directory)
        flanks={f['sequence_id']:f for row in cases for f in neighbors[row['reference_gene_id']]}
        write_fasta(directory/'flanks.faa',{key:proteins[f['ncbi_protein_id']] for key,f in flanks.items()})
        genome=RUN/'sensitivity_no_mask'/accession/'genome.fna'
        binary=find_binary('tblastn')
        args=[binary,'-query',str((directory/'flanks.faa').resolve()),'-subject',str(genome.resolve()),
              '-seg','no','-comp_based_stats','0','-num_threads','2','-max_target_seqs','100',
              '-evalue','1e-5','-outfmt','6 '+FIELDS]
        process=subprocess.run(args,capture_output=True,text=True,timeout=900)
        evidence=store.add('flank_tblastn',dict(stdout=process.stdout,stderr=process.stderr,exit_code=process.returncode),
                           status='ok' if process.returncode==0 else 'error',command=args)
        if process.returncode:raise RuntimeError('Flank search failed, not gene-loss evidence')
        hits=parse_hits(process.stdout)
        local=[]
        for row in cases:
            partial=not row['genome_homolog_found']
            center=dict(seqid=row['weaker_hit_seqid'] if partial else row['best_genome_seqid'],
                        start=row['weaker_hit_start'] if partial else row['best_genome_start'],
                        end=row['weaker_hit_end'] if partial else row['best_genome_end'],
                        strand=row['weaker_hit_strand'] if partial else row['strand'])
            gene=genes[row['reference_gene_id']];left,right=neighbors[row['reference_gene_id']]
            anchors={}
            for role,reference_flank in [('left',left),('right',right)]:
                eligible=[h for h in hits if h['qseqid']==reference_flank['sequence_id'] and compatible(h,center,gene,reference_flank,role)]
                anchors[role]=max(eligible,key=lambda h:h['bitscore']) if eligible else None
            result=dict(target_assembly=accession,reference_gene=row['reference_gene_name'],
                        reference_gene_id=row['reference_gene_id'],left_reference_gene=left['gene_name'],right_reference_gene=right['gene_name'],
                        target_seqid=center['seqid'],target_start=center['start'],target_end=center['end'],
                        partial_center_evidence=partial,left_anchor_supported=bool(anchors['left']),right_anchor_supported=bool(anchors['right']),
                        both_flanks_support=bool(anchors['left'] and anchors['right']),
                        evidence_id=evidence['id'],source_classification=row['classification'],
                        functional_gene_loss_inferred=False)
            local.append(result)
            store.add('neighborhood_assessment',dict(result=result,center=center,reference_center=gene,reference_flanks=[left,right],anchors=anchors),input_evidence=evidence['id'])
        results.extend(local);write_json(directory/'report.json',local);write_json(directory/'evidence.json',store.items)
        print(accession+': '+str(sum(r['both_flanks_support'] for r in local))+'/'+str(len(local))+' candidates supported by both flanks',flush=True)
    write_json(OUT/'synteny_candidate_tests.json',results)
    with (OUT/'synteny_candidate_tests.csv').open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)


if __name__=='__main__':main()
