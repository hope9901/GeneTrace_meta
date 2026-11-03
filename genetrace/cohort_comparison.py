"""Population context is not an error classifier. Frequencies require callable loci."""
import collections,re

COMPLEMENT=str.maketrans('ACGTRYMKBDHVN','TGCAYRKMVHDBN')
def revcomp(seq):return seq.translate(COMPLEMENT)[::-1]

def platform_group(row):
    publication=row.get('publication_supported_platforms','')
    text=publication if publication and publication!='N/A' else row.get('sequencing_group','')
    labels=[label for word,label in [('illumina','Illumina'),('nanopore','Nanopore'),('pacbio','PacBio'),('454','454'),('mgi','MGI'),('dnbseq','MGI'),('sanger','Sanger')] if word in text.lower()]
    return '+'.join(sorted(set(labels))) or 'N/A'

def anchor_support(hit,center,reference_center,reference_flank,role,max_distance=25000):
    if hit['sseqid']!=center['sseqid'] or hit['pident']<70 or hit['query_coverage']<.5 or hit['evalue']>1e-20:return False
    same=reference_center['strand']==center['strand'];expected=reference_flank['strand'] if same else ('-' if reference_flank['strand']=='+' else '+')
    if hit['strand']!=expected:return False
    before=(role=='left')==same
    return 0<center['start']-hit['end']<=max_distance if before else 0<hit['start']-center['end']<=max_distance

def distinct_loci(hits):
    result=[]
    for h in sorted(hits,key=lambda x:x['bitscore'],reverse=True):
        if not any(h['sseqid']==r['sseqid'] and h['strand']==r['strand'] and max(0,min(h['end'],r['end'])-max(h['start'],r['start'])+1)/min(h['end']-h['start']+1,r['end']-r['start']+1)>.5 for r in result):result.append(h)
    return result

def parse_cs(cs,target_start,query_start,query):
    """0-based coordinates; cs long tokens from minimap2, in target orientation."""
    pattern=r'=[A-Za-z]+|:\d+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}'
    tokens=re.findall(pattern,cs)
    if ''.join(tokens)!=cs:raise ValueError('Unsupported minimap2 cs string')
    t=target_start;q=query_start;mapping={};insertions=collections.defaultdict(str)
    for token in tokens:
        op=token[0]
        if op in '=:':
            n=len(token)-1 if op=='=' else int(token[1:])
            for i in range(n):mapping[t+i]=(query[q+i].upper(),q+i)
            t+=n;q+=n
        elif op=='*':mapping[t]=(token[2].upper(),q);t+=1;q+=1
        elif op=='+':insertions[t-1]+=token[1:].upper();q+=len(token)-1
        elif op=='-':
            for i in range(len(token)-1):mapping[t+i]=('-',None)
            t+=len(token)-1
        elif op=='~':t+=int(token[3:-2])
    return mapping,dict(insertions)

def repeat_bounds(reference,position):
    base=reference[position].upper();lo=position;hi=position
    while lo and reference[lo-1].upper()==base:lo-=1
    while hi+1<len(reference) and reference[hi+1].upper()==base:hi+=1
    return lo,hi,base

def repeat_allele(reference,position,mapping,insertions):
    lo,hi,base=repeat_bounds(reference,position)
    if lo==0 or hi+1>=len(reference):return None
    if any(p not in mapping for p in range(lo-1,hi+2)):return None
    if any(mapping[p][0] not in 'ACGT' for p in [lo-1,hi+1]):return None
    allele=insertions.get(lo-1,'')
    for p in range(lo,hi+1):
        b=mapping[p][0]
        if b not in 'ACGT-':return None
        allele+=('' if b=='-' else b)+insertions.get(p,'')
    return dict(reference_base=base,reference_repeat_length=hi-lo+1,allele=allele,repeat_length=len(allele) if all(c==base for c in allele) else None,alignment_dependent=True)

def summarize_marker(rows,focal_allele):
    callable_rows=[r for r in rows if r.get('allele') is not None];counts=collections.Counter(r['allele'] for r in callable_rows)
    return dict(callable_assemblies=len(callable_rows),uncallable_assemblies=len(rows)-len(callable_rows),allele_counts=dict(counts),focal_allele=focal_allele,focal_allele_count=counts[focal_allele],focal_allele_fraction=counts[focal_allele]/len(callable_rows) if callable_rows else None,error_inferred=False,biological_loss_inferred=False,interpretation='Frequency is descriptive population context, not proof of artifact or true variation.')

def variable_reference_positions(reference,alignments,span_start,span_end):
    """One marker per changed reference repeat block; callable assemblies only.

    Span is zero-based, half-open in the oriented focal locus. Not a VCF caller.
    """
    positions=[];blocks=set()
    for pos in range(span_start,span_end):
        lo,hi,_=repeat_bounds(reference,pos)
        if lo in blocks:continue
        for mapping,insertions in alignments:
            allele=repeat_allele(reference,pos,mapping,insertions)
            if allele is not None and allele['allele']!=reference[lo:hi+1]:
                blocks.add(lo);positions.append(pos);break
    return positions
