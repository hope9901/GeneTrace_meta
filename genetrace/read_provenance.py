"""Conservative authenticity review: regular tiling is not independent read proof."""
import collections

def assess_read_provenance(sam_text,quality_counts):
    starts=collections.defaultdict(set);templates=collections.Counter();nm=[];records=0
    for line in sam_text.splitlines():
        if not line or line.startswith('@'):continue
        f=line.split('\t');flag=int(f[1])
        if flag&3844 or f[2]=='*':continue
        records+=1;starts[(f[2],bool(flag&16))].add(int(f[3]))
        if int(f[8]):templates[abs(int(f[8]))]+=1
        value=next((int(x[5:]) for x in f[11:] if x.startswith('NM:i:')),None)
        if value is not None:nm.append(value)
    steps=collections.Counter()
    for positions in starts.values():
        ordered=sorted(positions);steps.update(b-a for a,b in zip(ordered,ordered[1:]))
    mode_step,mode_count=steps.most_common(1)[0] if steps else (None,0)
    step_fraction=mode_count/sum(steps.values()) if steps else 0
    template_fraction=max(templates.values())/sum(templates.values()) if templates else 0
    quality_fraction=max(quality_counts.values())/sum(quality_counts.values()) if quality_counts else 0
    zero_nm_fraction=sum(n==0 for n in nm)/len(nm) if nm else 0
    review=records>=100 and sum(steps.values())>=100 and step_fraction>=.95 and template_fraction>=.95 and quality_fraction>=.99 and zero_nm_fraction>=.999
    return dict(sampled_primary_records=records,modal_same_strand_start_step_bp=mode_step,modal_start_step_fraction=step_fraction,modal_absolute_template_length=templates.most_common(1)[0][0] if templates else None,modal_template_length_fraction=template_fraction,modal_quality_character_fraction=quality_fraction,zero_edit_distance_fraction=zero_nm_fraction,provenance_review_required=review,independent_physical_read_validation_established=False,synthetic_data_confirmed=False,note='Conjunction of regular tiling, fixed template length, near-constant quality and exact matches triggers review. It does not prove simulation. Registration and BioSample matching alone do not establish independent physical sequencing.')
