#!/usr/bin/env python3
"""Join DNA/synteny checks with core-tree neighbors and technical strata."""
import collections,csv,json,pathlib,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.cohort_comparison import platform_group
OUT=ROOT/'outputs/ecoli100_candidate_audit_20261102';metadata=json.loads((ROOT/'data/bacteria/cohort.json').read_text());by={r['assembly_accession']:r for r in metadata};counts={r['group']:r for r in csv.DictReader((ROOT/'outputs/ecoli100_gene_groups_20261102/group_counts.csv').open())};neighbors=json.loads((ROOT/'outputs/ecoli100_context_20261102/nearest_neighbors.json').read_text());audits=json.loads((OUT/'results.json').read_text());results=[]
for a in audits:
 acc=a['assembly_accession'];group=counts[a['group']];row=by[acc];nearest=[dict(**n,protein_cluster_copies=int(group[n['assembly']])) for n in neighbors[acc]];strata=[]
 for field in ['bioproject','sequencing_group']:
  focal=row.get(field,'N/A') if field=='bioproject' else platform_group(row)
  same=[m for m in metadata if m['assembly_accession']!=acc and (m.get(field,'N/A') if field=='bioproject' else platform_group(m))==focal] if focal!='N/A' else []
  strata.append(dict(field=field,focal=focal,other_assemblies=len(same),missing_protein_cluster=sum(int(group[m['assembly_accession']])==0 for m in same),interpretation='Descriptive; related strains/project batches are not independent and no causal platform effect is inferred.'))
 results.append(dict(assembly=acc,gene=a['reference_gene'],group=a['group'],dna_classification=a['classification'],query_coverage=a['best_genome_hit']['query_coverage'] if a['best_genome_hit'] else None,identity=a['best_genome_hit']['pident'] if a['best_genome_hit'] else None,both_syntenic_flanks=a.get('both_flanks',False),recorded_disruption=any('frameshift' in str(x['attrs']).lower() or x['attrs'].get('pseudo')=='true' for x in a.get('overlapping_cds',[])),nearest_neighbors=nearest,technical_strata=strata,raw_run_candidate=a['run_candidate']['run_accession'],causal_diagnosis='unresolved'))
(OUT/'context.json').write_text(json.dumps(results,indent=2)+'\n')
lines=['# 대장균 전체 단백질 후보와 계통·주변 유전자 맥락','', '100개 487,647개 단백질 전부를 90% AA identity/80% 양방향 coverage 조건으로 묶었다. 9,233개 상동성 후보 그룹, 3,027개 single-copy 공통 그룹, 높은 빈도 그룹의 누락 조합 915건을 찾았다. 이는 확정 orthogroup 분석이 아니다. 모든 입력 단백질의 1회 배정을 확인했다.','','공통 single-copy 후보 30개 단백질, 정렬 7,664 columns로 잠정 계통수를 만들었다. 재조합 보정·gene-tree reconciliation은 하지 않았으며 동일 단백질에서 생성한 계통수는 독립 검증이 아니다.','','915건 중 작은 WGS run이 연결되는 표본을 우선하여 9건을 DNA 검색·양쪽 주변 유전자 검사로 확인했다. 이 선택은 편의 표본이며 검출 정확도 추정용이 아니다.','','| Assembly | 유전자 | DNA identity | query coverage | 양쪽 주변 유전자 | 주석 disruption | 가까운 5개 중 단백질 그룹 누락 |','|---|---|---|---|---|---|']
for r in results:lines.append(f"| {r['assembly']} | {r['gene']} | {r['identity']} | {r['query_coverage']:.3f} | {r['both_syntenic_flanks']} | {r['recorded_disruption']} | {sum(n['protein_cluster_copies']==0 for n in r['nearest_neighbors'])}/5 |")
lines+=['','tonB는 DNA 상동 서열·단백질 주석이 존재해 단백질 클러스터 미검출을 유전자 소실로 해석할 수 없다. nimT는 frameshifted 주석을 갖고 있어 raw read와 전체 코딩 모델 재검증이 필요하다. 낮은 query coverage의 HSP는 전체 유전자 또는 deletion breakpoint를 뜻하지 않는다.','','기술 조합 및 BioProject별 집계와 가장 가까운 5개 표본의 거리·그룹 상태는 context.json에 있다. 유의성 검정·계통 보정 회귀를 한 결과가 아니며 플랫폼 편향을 입증하지 않는다. raw read 결과는 ../ecoli100_raw_validation_20261102/에서 별도 기록한다. 모든 원인은 미해결로 남긴다.']
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print('Joined contexts',len(results))
