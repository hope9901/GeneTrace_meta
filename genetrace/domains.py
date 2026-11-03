"""Pfam profile-based evidence; a domain hit is not an ortholog or functional gene."""
import csv,json,re,subprocess
from pathlib import Path
from .engine import ROOT
from .io import read_fasta,read_gff,sha256,write_fasta,write_json

# NCBI translation table 11 has the standard internal codon translations.
_BASES='TCAG';_AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODONS={a+b+c:aa for (a,b,c),aa in zip(((a,b,c) for a in _BASES for b in _BASES for c in _BASES),_AA)}

def six_frame_orfs(records,min_aa=30):
    peptides={};mapping={};i=0
    for seqid,dna in records.items():
        reverse=dna.translate(str.maketrans('ACGTRYMKBDHVN','TGCAYRKMVHDBN'))[::-1]
        for strand,seq in [('+',dna),('-',reverse)]:
            for offset in range(3):
                aa=''.join(CODONS.get(seq[p:p+3],'X') for p in range(offset,len(seq)-2,3));cursor=0
                for part in aa.split('*'):
                    if len(part)>=min_aa:
                        i+=1;name=f'orf_{i:07d}';s=offset+cursor*3+1;e=s+len(part)*3-1
                        if strand=='-':s,e=len(dna)-e+1,len(dna)-s+1
                        peptides[name]=part;mapping[name]=dict(seqid=seqid,start=s,end=e,strand=strand,frame=offset+1,length=len(part))
                    cursor+=len(part)+1
    return peptides,mapping

def domain_coordinates(orf,aa_start,aa_end):
    if not 1<=aa_start<=aa_end<=orf['length']:raise ValueError('Domain outside translated ORF')
    if orf['strand']=='+':start,end=orf['start']+(aa_start-1)*3,orf['start']+aa_end*3-1
    else:start,end=orf['end']-aa_end*3+1,orf['end']-(aa_start-1)*3
    return dict(seqid=orf['seqid'],start=start,end=end,strand=orf['strand'])

def parse_domtbl(text,models):
    ga={m['accession'].split('.')[0]:[float(x) for x in m['ga'].strip(';').split()] for m in models};hits=[]
    for line in text.splitlines():
        if not line or line.startswith('#'):continue
        f=line.split(maxsplit=22)
        if len(f)<22:raise ValueError('Malformed HMMER domain table')
        pfam=f[4].split('.')[0]
        if pfam not in ga:raise ValueError('Unrecognized profile')
        if float(f[7])<ga[pfam][0] or float(f[13])<ga[pfam][1]:continue
        coverage=(int(f[16])-int(f[15])+1)/int(f[5])
        hits.append(dict(protein_id=f[0],protein_length=int(f[2]),profile=pfam,profile_version=f[4],profile_name=f[3],profile_length=int(f[5]),sequence_score=float(f[7]),domain_score=float(f[13]),independent_evalue=float(f[12]),hmm_start=int(f[15]),hmm_end=int(f[16]),aa_start=int(f[17]),aa_end=int(f[18]),profile_coverage=coverage,coverage_class='near_complete' if coverage>=.7 else 'partial',gathering_threshold_passed=True))
    return hits

def scan(fasta,directory,models=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    models=models or json.loads((ROOT/'data/domain_models/manifest.json').read_text());panel=ROOT/'data/domain_models/panel.hmm'
    signature=dict(fasta_sha256=sha256(fasta),profile_panel_sha256=sha256(panel))
    old=directory/'scan_manifest.json'
    if old.exists() and json.loads(old.read_text()).get('inputs')==signature and (directory/'hits.json').exists():return json.loads((directory/'hits.json').read_text())
    binary=ROOT/'tools/readqc-env/bin/hmmsearch';cmd=[str(binary),'--cpu','1','--cut_ga','--noali','--domtblout',str(directory/'domains.domtbl'),str(panel),str(fasta)]
    with (directory/'hmmsearch.stdout').open('w') as out,(directory/'hmmsearch.stderr').open('w') as err:
        p=subprocess.run(cmd,stdout=out,stderr=err,timeout=900)
    if p.returncode:raise RuntimeError('HMMER failed; absence cannot be inferred')
    hits=parse_domtbl((directory/'domains.domtbl').read_text(),models)
    write_json(directory/'hits.json',hits);write_json(old,dict(command=cmd,version=subprocess.check_output([str(binary),'-h'],text=True).splitlines()[1],inputs=signature,thresholds='Pfam sequence AND domain gathering cutoffs; coverage>=0.7 marked near_complete, shorter hits retained separately'))
    return hits

def matching_cds(position,features):
    return [f for f in features if f['type']=='CDS' and f['seqid']==position['seqid'] and f['strand']==position['strand'] and max(0,min(f['end'],position['end'])-max(f['start'],position['start'])+1)/(position['end']-position['start']+1)>=.5]

def genomic_assessment(hit,position,features,protein_hits):
    cds=matching_cds(position,features);profile=hit['profile'];protein_ids={f['attrs'].get('protein_id') for f in cds}
    matching=[h for h in protein_hits if h['profile']==profile and h['protein_id'] in protein_ids]
    overlap=[f for f in features if f['seqid']==position['seqid'] and f['start']<=position['end'] and f['end']>=position['start'] and f['type'] in ['gene','pseudogene','CDS','sequence_alteration']]
    disrupted=any(f['type']=='pseudogene' or f['attrs'].get('pseudo')=='true' or re.search(r'frameshift|pseudogen|interrupted',json.dumps(f['attrs']),re.I) for f in overlap)
    if disrupted:label='pseudogene_or_disrupted_model_unresolved'
    elif cds and matching:label='protein_and_genome_domain_supported'
    elif cds:label='annotated_locus_domain_discrepancy'
    else:label='unannotated_genomic_domain_candidate'
    return dict(**position,profile=profile,domain_score=hit['domain_score'],profile_coverage=hit['profile_coverage'],coverage_class=hit['coverage_class'],overlapping_protein_ids=sorted(x for x in protein_ids if x),classification=label,recorded_disruption=disrupted,overlapping_gene_labels=sorted({f['attrs'].get('gene',f['attrs'].get('locus_tag','N/A')) for f in overlap if f['type'] in ['gene','pseudogene']}),biological_loss_inferred=False,functional_status='not_assessed')

def domain_snapshot(assembly,profile=None,protein_ids=None):
    """Read bounded existing results; do not report an unscanned assembly as a zero."""
    directory=ROOT/'outputs/domain_bacteria_20261101'
    expanded=ROOT/'outputs/ecoli100_domains_20261102'
    if assembly.startswith('GCF_') and (expanded/'status.json').exists():
        expanded_state=json.loads((expanded/'status.json').read_text())
        if expanded_state.get('stage')=='complete':directory=expanded
    try:status=json.loads((directory/'status.json').read_text())
    except (FileNotFoundError,json.JSONDecodeError):status={'stage':'not_available'}
    def rows(name):
        path=directory/name
        if not path.exists():return []
        with path.open(newline='') as f:return [r for r in csv.DictReader(f) if r.get('assembly_accession')==assembly and (not profile or r.get('profile')==profile)]
    counts=rows('domain_counts.csv');members=rows('domain_members.csv')
    if protein_ids is not None:members=[h for h in members if h.get('protein_id') in set(protein_ids)]
    return dict(assembly_accession=assembly,status=status,assembly_scanned=bool(counts),panel_only=True,counts=counts,members=members[:250],matching_domain_instances=len(members),display_limit=250,notes=['Only six selected Pfam profiles were searched. No panel hit does not establish domain or gene absence.','Domain instances and proteins carrying domains are separate counts; neither is an orthogroup copy-number estimate.'])
