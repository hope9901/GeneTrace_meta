"""Local alignment observations, deliberately separated from causal diagnosis.

Coordinates are 1-based inclusive. SAM SEQ is already oriented to reference.
No variant calls or independent-molecule claims are made.
"""
import collections,re,statistics

POLICY=dict(min_mapq=20,min_depth=10,min_alt_records=5,coherent_fraction=.8,mixed_fraction=.2,min_strand_records=2,clip_bp=20,clip_radius=60,diffuse_mismatch_fraction=.08,min_boundary_depth=5,large_gap_bp=50)

def cigar_ops(cigar):
    ops=[(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar)]
    if not ops or ''.join(str(n)+op for n,op in ops)!=cigar or any(n<=0 for n,_ in ops):raise ValueError('Invalid CIGAR')
    return ops

def sam_records(text):
    for line in text.splitlines():
        if not line or line.startswith('@'):continue
        f=line.split('\t')
        if len(f)<11:raise ValueError('Malformed SAM record')
        yield dict(name=f[0],flag=int(f[1]),seqid=f[2],start=int(f[3]),mapq=int(f[4]),cigar=f[5],sequence=f[9].upper(),quality=f[10],tags=f[11:])

def analyze_sam(text,reference,seqid,start,end,policy=None,partial_input=False):
    policy=dict(POLICY,**(policy or {}))
    if not 1<=start<=end<=len(reference):raise ValueError('Invalid region')
    length=end-start+1;counts={p:collections.Counter() for p in range(start,end+1)};strand_counts=collections.Counter();events=collections.defaultdict(lambda:collections.Counter());event_reads=collections.defaultdict(set)
    observed_qualities=collections.Counter();observations=[];excluded=collections.Counter();raw_records=0;supplementary=0
    flank=100;left=max(1,start-flank);right=min(len(reference),end+flank);depth_context=collections.Counter();boundary_clips=[];seen=set();alignment_segments={};record_strands={};anchor_records=collections.defaultdict(set)
    for record in sam_records(text):
        raw_records+=1;flag=record['flag']
        if flag&2048:supplementary+=1
        if flag&3844:excluded['flag_filtered']+=1;continue
        if record['seqid']!=seqid:excluded['other_reference']+=1;continue
        if record['mapq']<policy['min_mapq'] or record['mapq']==255:excluded['mapq_filtered_or_unavailable']+=1;continue
        if record['sequence']=='*' or record['cigar']=='*':excluded['missing_sequence_or_cigar']+=1;continue
        key=(record['name'],flag&192) # count a primary segment once; mates remain separate
        if key in seen:excluded['duplicate_primary_segment']+=1;continue
        seen.add(key);ops=cigar_ops(record['cigar']);sequence=record['sequence'];quality=record['quality']
        if sum(n for n,op in ops if op in 'MIS=X')!=len(sequence):raise ValueError('CIGAR/query length mismatch')
        if quality!='*' and len(quality)!=len(sequence):raise ValueError('Quality/query length mismatch')
        strand='-' if flag&16 else '+';r=record['start'];q=0;matches=mismatches=aligned=deletion_bases=insertion_bases=0;covered=set();large_gaps=[];clips=[];segments=[]
        for n,op in ops:
            if op in 'M=X':
                segments.append((r,r+n-1))
                for i in range(n):
                    pos=r+i
                    if not 1<=pos<=len(reference):raise ValueError('Alignment outside reference')
                    if left<=pos<=right:covered.add(pos);depth_context[pos]+=1;anchor_records[pos].add(key)
                    if start<=pos<=end:
                        base=sequence[q+i];ref=reference[pos-1].upper()
                        if quality!='*':observed_qualities[quality[q+i]]+=1
                        if base in 'ACGT' and ref in 'ACGT':
                            counts[pos][base]+=1;strand_counts[(pos,base,strand)]+=1;aligned+=1
                            if base==ref:matches+=1
                            else:mismatches+=1
                r+=n;q+=n
            elif op=='I':
                anchor=r-1
                if start<=anchor<=end:
                    inserted=sequence[q:q+n];event=('I',anchor,inserted);events[event][strand]+=1;event_reads[event].add(key);insertion_bases+=n
                q+=n
            elif op=='D':
                overlap=max(0,min(end,r+n-1)-max(start,r)+1);deletion_bases+=overlap
                if overlap:
                    event=('D',r,reference[r-1:r+n-1]);events[event][strand]+=1;event_reads[event].add(key)
                if n>=policy['large_gap_bp']:large_gaps.append((r,r+n-1))
                r+=n
            elif op=='N':
                large_gaps.append((r,r+n-1));r+=n
            elif op in 'SH':
                if n>=policy['clip_bp']:clips.append(dict(position=r,length=n,type=op))
                if op=='S':q+=n
            elif op=='P':pass
        if not any(start<=p<=end for p in covered) and not deletion_bases:continue
        boundary=[c for c in clips if min(abs(c['position']-start),abs(c['position']-end),abs(c['position']-(start+end)//2))<=policy['clip_radius']]
        boundary_clips.extend(dict(read=record['name'],strand=strand,**c) for c in boundary)
        spans=all(p in covered for p in [left,right]) and not any(a<=right and b>=left for a,b in large_gaps)
        alignment_segments[key]=segments;record_strands[key]=strand
        observations.append(dict(read=record['name'],strand=strand,mapq=record['mapq'],start=record['start'],end=r-1,aligned_callable_bases=aligned,mismatches=mismatches,insertion_bases=insertion_bases,deletion_bases=deletion_bases,spans_region_plus100=spans,boundary_clip=bool(boundary),cigar=record['cigar']))
    base_rows=[];coherent=[];mixed=[]
    for pos,c in counts.items():
        ref=reference[pos-1].upper();depth=sum(c.values());alt,n=max(((b,n) for b,n in c.items() if b!=ref),key=lambda x:x[1],default=(None,0));frac=n/depth if depth else 0
        forward=strand_counts[(pos,alt,'+')] if alt else 0;reverse=strand_counts[(pos,alt,'-')] if alt else 0
        row=dict(position=pos,reference=ref,depth=depth,counts=dict(c),top_alt=alt,alt_records=n,alt_fraction=frac,alt_forward=forward,alt_reverse=reverse)
        base_rows.append(row)
        supported=depth>=policy['min_depth'] and n>=policy['min_alt_records'] and min(forward,reverse)>=policy['min_strand_records']
        if supported and frac>=policy['coherent_fraction']:coherent.append(row)
        elif supported and policy['mixed_fraction']<=frac<policy['coherent_fraction']:mixed.append(row)
    indels=[]
    for event,ss in events.items():
        kind,pos,seq=event;anchor_left=pos-1 if kind=='D' else pos;anchor_right=pos+len(seq) if kind=='D' else pos+1
        callable_keys=anchor_records[anchor_left]&anchor_records[anchor_right]
        support=event_reads[event]&callable_keys;denominator=len(callable_keys);count=len(support);fraction=count/denominator if denominator else 0
        ss=collections.Counter(record_strands[k] for k in support)
        row=dict(kind=kind,position=pos,sequence=seq,length=len(seq),raw_supporting_records=len(event_reads[event]),supporting_records=count,spanning_record_denominator=denominator,support_fraction=fraction,forward=ss['+'],reverse=ss['-'])
        row['coherent']=denominator>=policy['min_depth'] and count>=policy['min_alt_records'] and fraction>=policy['coherent_fraction'] and min(ss['+'],ss['-'])>=policy['min_strand_records'];indels.append(row)
    depths=[r['depth'] for r in base_rows];aligned=sum(r['aligned_callable_bases'] for r in observations);mismatch=sum(r['mismatches'] for r in observations);spanning=[r for r in observations if r['spans_region_plus100']]
    clipping_fraction=sum(r['boundary_clip'] for r in observations)/len(observations) if observations else 0
    left_depth=statistics.mean(depth_context[p] for p in range(left,start)) if start>left else 0;right_depth=statistics.mean(depth_context[p] for p in range(end+1,right+1)) if right>end else 0
    n_fraction=reference[start-1:end].upper().count('N')/length
    signals=[]
    if n_fraction>=.1:signals.append('assembly_gap_or_ambiguity')
    if not observations or statistics.median(depths)<policy['min_depth']:signals.append('insufficient_read_support')
    if coherent:signals.append('coherent_substitution_disagreement')
    if any(r['coherent'] for r in indels):signals.append('coherent_indel_disagreement')
    if mixed:signals.append('mixed_allele_signal')
    if not spanning and min(left_depth,right_depth)>=policy['min_boundary_depth'] and clipping_fraction>=.2:signals.append('junction_discontinuity_signal')
    rate=mismatch/aligned if aligned else None
    if rate is not None and rate>policy['diffuse_mismatch_fraction'] and not coherent:signals.append('diffuse_read_disagreement')
    if not signals:signals=['read_assembly_supported']
    qualities_unresolved=not observed_qualities or set(observed_qualities)=={'!'}
    return dict(schema_version=1,region=dict(seqid=seqid,start=start,end=end),partial_input=partial_input,policy=policy,signals=signals,causal_diagnosis='unresolved',biological_loss_inferred=False,functional_status='not_assessed',variant_calls_made=False,quality_information_unresolved=qualities_unresolved,base_quality_threshold=0,metrics=dict(raw_alignment_records=raw_records,usable_primary_records=len(observations),excluded=dict(excluded),supplementary_records=supplementary,mean_depth=statistics.mean(depths),median_depth=statistics.median(depths),covered_fraction=sum(x>0 for x in depths)/length,mismatch_fraction=rate,aligned_callable_bases=aligned,mismatches=mismatch,insertion_bases=sum(r['insertion_bases'] for r in observations),deletion_bases=sum(r['deletion_bases'] for r in observations),spanning_plus100_records=len(spanning),spanning_forward=sum(r['strand']=='+' for r in spanning),spanning_reverse=sum(r['strand']=='-' for r in spanning),boundary_clipped_fraction=clipping_fraction,left_flank_mean_depth=left_depth,right_flank_mean_depth=right_depth,n_fraction=n_fraction),coherent_substitutions=coherent,mixed_alleles=mixed,indels=indels,per_base=base_rows,reads=observations,boundary_clips=boundary_clips,limitations=['Signals are not confirmed sequencing errors, assembly errors or biological variants.','Thresholds are provisional flags and are not calibrated to every sequencing technology.','Records may not be independent molecules; quality ! does not establish poor sequencing.','Reference skips are not deletion evidence; RNA/splicing is not modeled as a genomic variant.','Indels are alignment representations, not left-normalized variant calls; repeat-shifted events may be split.']+(['Partial input cannot diagnose low sequencing coverage or biological absence.'] if partial_input else []))
