"""Systematic 64-locus diagnostic panel on real, unmodified assemblies.

This pilot does not replace the all-gene OrthoFinder job or estimate genome-wide
loss prevalence. It tests whether apparent protein absences can reflect annotation
scope/model differences. Homology is not asserted to be orthology or function.
"""
import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter,defaultdict
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import read_fasta,read_gff,write_fasta,write_json,sha256
from genetrace.engine import EvidenceStore,find_binary,parse_hits,FIELDS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/loss_artifact_tests_20261030'
RUN=ROOT/'runs/loss_artifact_tests_20261030'
SOURCE=ROOT/'data/orthology/source'
REFERENCE='GCA_000146045.2'
CASES=[REFERENCE,'GCA_903819185.2','GCA_903819125.2','GCA_000190255.1','GCA_000326005.1']
CRITERIA=dict(min_identity=90.0,min_single_hsp_query_coverage=.90,max_evalue=1e-20)


def nonnuclear(text):return bool(re.search(r'mitochondri|plastid|chloroplast|plasmid',text,re.I))


def classify(protein_found,genome_found,cds_overlap,tool_ok=True):
    if not tool_ok:return 'search_failed_no_conclusion'
    if genome_found and not protein_found:
        return 'genomic_homolog_missing_from_proteome_with_CDS_overlap' if cds_overlap else 'unannotated_genomic_homolog'
    if genome_found and protein_found:return 'protein_and_genome_homolog_detected'
    if protein_found:return 'protein_homolog_only_genome_search_unresolved'
    return 'no_qualifying_homolog_in_panel_search_unresolved'


def panel():
    proteins=read_fasta(SOURCE/REFERENCE/'protein.faa.gz')
    features=read_gff(SOURCE/REFERENCE/'annotation.gff.gz')
    cds=Counter(f['attrs'].get('protein_id') for f in features if f['type']=='CDS')
    rows=list(csv.DictReader((ROOT/'data/orthology/gene_maps'/(REFERENCE+'.tsv')).open(),delimiter='\t'))
    grouped=defaultdict(list)
    for row in rows:
        protein=row['ncbi_protein_id'];sequence=proteins[protein]
        if cds[protein]==1 and 150<=len(sequence)<=800 and '*' not in sequence and 'X' not in sequence:
            grouped[row['seqid']].append(row)
    # The GenBank S288C reference input contains 16 nuclear chromosome records.
    genome=read_fasta(SOURCE/REFERENCE/'genome.fna.gz')
    if len(genome)!=16 or len(grouped)!=16:raise ValueError('Unexpected S288C nuclear chromosome set')
    chosen=[]
    for seqid,genes in sorted(grouped.items()):
        genes.sort(key=lambda r:int(r['start']))
        if len(genes)<4:raise ValueError('Too few eligible genes on '+seqid)
        for quantile in (.125,.375,.625,.875):chosen.append(genes[int((len(genes)-1)*quantile)])
    if len({r['sequence_id'] for r in chosen})!=64:raise ValueError('Panel is not 64 distinct loci')
    write_fasta(RUN/'panel.faa',{r['sequence_id']:proteins[r['ncbi_protein_id']] for r in chosen})
    write_json(OUT/'panel_manifest.json',dict(reference=REFERENCE,genes=chosen,criteria=CRITERIA,
                selection='Four coordinate-quantile single-CDS loci per nuclear chromosome, 150–800 aa, no stop/X; fixed before target searches',
                query_sha256=sha256(RUN/'panel.faa'),scope='64-locus diagnostic pilot, not all-gene loss inference',
                source_unchanged=True,orthofinder_inputs_unchanged=True))
    return {r['sequence_id']:r for r in chosen}


def run_search(store,name,subject,sensitivity=False):
    binary=find_binary(name)
    store.command(name+'_version',[binary,'-version'])
    args=[binary,'-query',str((RUN/'panel.faa').resolve()),'-subject',str(subject.resolve()),
          '-evalue','1e-5','-seg','no' if sensitivity else 'yes','-num_threads','2','-max_target_seqs','100','-outfmt','6 '+FIELDS]
    if sensitivity:args+=['-comp_based_stats','0']
    # Record command, stdout/stderr and status, retaining failures distinctly.
    import subprocess
    try:
        process=subprocess.run(args,capture_output=True,text=True,timeout=900,shell=False)
        record=store.add(name,dict(stdout=process.stdout,stderr=process.stderr,exit_code=process.returncode),
                         status='ok' if process.returncode==0 else 'error',command=args)
    except (OSError,subprocess.TimeoutExpired) as exc:
        record=store.add(name,dict(error=str(exc)),status='error',command=args)
    if record['status']!='ok':raise RuntimeError('Search failed, not biological absence: '+name)
    hits=parse_hits(record['result']['stdout'])
    for hit in hits:
        hit['passes_thresholds']=(hit['pident']>=90 and hit['query_coverage']>=.90 and hit['evalue']<=1e-20)
    store.add(name+'_parsed',dict(criteria=CRITERIA,hits=hits),input_evidence=record['id'])
    return hits


def run_case(accession,panel_rows,metadata,sensitivity=False):
    directory=(RUN/'sensitivity_no_mask'/accession) if sensitivity else RUN/accession
    if (directory/'diagnostic.json').exists():
        old=json.loads((directory/'diagnostic.json').read_text())
        if old['panel_sha256']!=sha256(RUN/'panel.faa'):raise ValueError('Cannot reuse results for a changed panel')
        return old
    store=EvidenceStore(directory)
    source=SOURCE/accession
    for original,target in [('genome.fna.gz','genome.fna'),('protein.faa.gz','protein.faa')]:
        with gzip.open(source/original,'rb') as read,(directory/target).open('wb') as write:
            import shutil
            shutil.copyfileobj(read,write)
    source_manifest=json.loads((source/'source_manifest.json').read_text())
    store.add('source_manifest',dict(assembly=accession,files=source_manifest,source_unchanged=True,panel_sha256=sha256(RUN/'panel.faa')))
    features=read_gff(source/'annotation.gff.gz')
    protein_hits=run_search(store,'blastp',directory/'protein.faa',sensitivity)
    genome_hits=run_search(store,'tblastn',directory/'genome.fna',sensitivity)
    protein_by_query=defaultdict(list);genome_by_query=defaultdict(list)
    for hit in protein_hits:
        if hit['passes_thresholds']:protein_by_query[hit['qseqid']].append(hit)
    for hit in genome_hits:
        if hit['passes_thresholds']:genome_by_query[hit['qseqid']].append(hit)
    results=[]
    for query,gene in panel_rows.items():
        proteins=protein_by_query[query];genomes=genome_by_query[query]
        best=max(genomes,key=lambda h:h['bitscore']) if genomes else None
        cds=[]
        if best:
            for f in features:
                if f['type']=='CDS' and f['seqid']==best['sseqid'] and f['strand']==best['strand']:
                    overlap=max(0,min(f['end'],best['end'])-max(f['start'],best['start'])+1)/(best['end']-best['start']+1)
                    if overlap>=.5:cds.append(f['attrs'].get('ID','unnamed_CDS'))
        result=dict(target_assembly=accession,reference_query_id=query,reference_gene_id=gene['gene_id'],
                    reference_gene_name=gene['gene_name'],reference_protein=gene['ncbi_protein_id'],
                    protein_homolog_found=bool(proteins),genome_homolog_found=bool(genomes),
                    protein_hit_count=len(proteins),genome_hit_count=len(genomes),
                    best_genome_seqid=best['sseqid'] if best else '',best_genome_start=best['start'] if best else '',
                    best_genome_end=best['end'] if best else '',strand=best['strand'] if best else '',
                    genome_identity=best['pident'] if best else '',genome_query_coverage=best['query_coverage'] if best else '',
                    overlapping_cds=';'.join(cds),classification=classify(bool(proteins),bool(genomes),bool(cds)),
                    biological_loss_inferred=False,
                    functional_status='not_assessed',orthology_status='not_yet_verified_by_full_orthogroups',
                    alignment_has_internal_stop=('*' in best['sseq']) if best else '')
        results.append(result)
    counts=Counter(r['classification'] for r in results)
    report=dict(assembly_accession=accession,strain=metadata[accession]['resolved_strain'],
                search_mode='SEG off; composition correction off' if sensitivity else 'SEG on; default composition correction',
                assembly_level=metadata[accession]['assembly_level'],proteins_in_submission=metadata[accession]['number_of_proteins'],
                queries=len(panel_rows),protein_positive=sum(r['protein_homolog_found'] for r in results),
                genome_positive=sum(r['genome_homolog_found'] for r in results),classifications=dict(counts),
                panel_sha256=sha256(RUN/'panel.faa'),results=results,
                limitations=['Selected cases are diagnostic extremes, not a prevalence estimate.',
                             'Strong homolog detection is not proof of orthology or intact function.',
                             'Single-HSP thresholds can miss divergent/spliced genes; negative results remain unresolved.',
                             'All original sequences and GFF are unmodified; no synthetic perturbation in these cases.'])
    store.add('panel_assessment',dict(counts=dict(counts),biological_loss_inferred=False))
    write_json(directory/'diagnostic.json',report);write_json(directory/'evidence.json',store.items)
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--sensitivity',action='store_true',help='Paired no-SEG/no-composition-correction sensitivity check; original thresholds unchanged')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);RUN.mkdir(parents=True,exist_ok=True)
    queries=panel()
    metadata={r['assembly_accession']:r for r in json.loads((ROOT/'outputs/orthology_20261030/selected_assemblies.json').read_text())}
    reports=[]
    for accession in CASES:
        print('Testing unmodified natural assembly '+accession,flush=True)
        report=run_case(accession,queries,metadata,args.sensitivity);reports.append(report)
        print(json.dumps({k:v for k,v in report.items() if k not in ('results','limitations')},indent=2),flush=True)
        write_json(OUT/('natural_sensitivity_reports.json' if args.sensitivity else 'natural_default_reports.json'),reports)
    flat=[r for report in reports for r in report['results']]
    reference=reports[0]
    callable_queries={r['reference_query_id'] for r in reference['results'] if r['protein_homolog_found'] and r['genome_homolog_found']}
    for row in flat:row['reference_control_callable']=row['reference_query_id'] in callable_queries
    prefix='natural_sensitivity' if args.sensitivity else 'natural_default'
    with (OUT/(prefix+'_locus_results.csv')).open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    write_json(OUT/(prefix+'_calibration.json'),dict(queries=64,reference_callable=len(callable_queries),
                    noncallable_queries=sorted(set(queries)-callable_queries),
                    masked_or_unresolved_reference_queries_are_not_loss_calls=True))
    if args.sensitivity and (reference['protein_positive']!=64 or reference['genome_positive']!=64):
        raise ValueError('Reference positive control did not recover all 64 loci')
    write_json(OUT/(prefix+'_case_summary.json'),[{k:v for k,v in report.items() if k not in ('results','limitations')} for report in reports])
    print('All five natural-case diagnostics complete; no biological loss calls made.',flush=True)


if __name__=='__main__':main()
