"""Choose a scalable core; every non-core genome remains in the final run.

Bottom-k sketches of exact protein sequences are only an initialization heuristic,
not an evolutionary tree, and do not remove genes or final-cohort genomes.
"""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import read_fasta, write_json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/orthology_20261030'
BASE=ROOT/'data/orthology'


def similarity(a,b,k=256):
    sample=sorted(a|b)[:k]
    return sum(x in a and x in b for x in sample)/len(sample)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--core-size',type=int,default=64)
    args=parser.parse_args()
    rows=json.loads((OUT/'all_gene_input_summary.json').read_text())
    if not (2<=args.core_size<len(rows)):
        raise ValueError('Core size must be between 2 and cohort size minus one')
    sketches={}
    for index,row in enumerate(rows,1):
        accession=row['assembly_accession']
        sequences=read_fasta(BASE/'proteomes'/(accession+'.faa'))
        hashes={int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),'big') for s in sequences.values()}
        sketches[accession]=set(sorted(hashes)[:256])
        if index%40==0:
            print('Full-proteome sketches computed: %d/%d'%(index,len(rows)),flush=True)
    eligible={r['assembly_accession'] for r in rows if not r['sparse_annotation_warning']}
    selected=['GCA_000146045.2']
    if selected[0] not in eligible:
        raise ValueError('S288C reference missing from core candidates')
    closest={a:similarity(sketches[a],sketches[selected[0]]) for a in eligible-set(selected)}
    decisions=[dict(assembly_accession=selected[0],core_order=1,reason='S288C reference anchor')]
    while len(selected)<args.core_size:
        next_id=min(closest,key=lambda a:(closest[a],a))
        decisions.append(dict(assembly_accession=next_id,core_order=len(selected)+1,
                              max_similarity_to_previous_core=closest[next_id],
                              reason='Greedy farthest exact-protein-sequence overlap sketch; not phylogenetic distance'))
        selected.append(next_id);closest.pop(next_id)
        for accession in closest:
            closest[accession]=max(closest[accession],similarity(sketches[accession],sketches[next_id]))
    core=set(selected)
    for name in ('core_proteomes','additional_proteomes'):
        directory=BASE/name
        directory.mkdir(parents=True,exist_ok=True)
        if list(directory.glob('*.faa')):
            raise FileExistsError('Core split already exists; do not silently change a running cohort')
    for row in rows:
        accession=row['assembly_accession']
        destination=BASE/('core_proteomes' if accession in core else 'additional_proteomes')/(accession+'.faa')
        # Symlinks avoid duplicate raw sequence storage; source files are immutable for this run.
        destination.symlink_to((BASE/'proteomes'/(accession+'.faa')).resolve())
    write_json(OUT/'core_selection.json',dict(core=selected,additional=sorted(set(sketches)-core),
               decisions=decisions,heuristic='bottom-256 exact protein sequence BLAKE2b hashes; greedy farthest',
               all_genes_used=True,all_retained_genomes_used=True,
               caveat='Core initialization is heuristic. Sparse proteome is excluded from core only, retained in additional cohort.'))
    print('Core %d + additional %d = %d retained genomes'%(len(core),len(rows)-len(core),len(rows)),flush=True)


if __name__=='__main__':
    main()
