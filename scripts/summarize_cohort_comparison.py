#!/usr/bin/env python3
"""Render the completed population pilot without causal error claims."""
import json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
out=ROOT/'outputs/cohort_comparison_20261102/YKL094W'
s=json.loads((out/'summary.json').read_text())
lines=['# HN1/YJU3 동일 종 1:N 비교','',f"전체 {s['total_assemblies']}개 중 단일 후보 위치와 양쪽 주변 유전자 배치가 지지된 유전체는 {s['unambiguous_synteny_supported']}개입니다. HN1을 제외한 278개를 빈도 계산 대상으로 삼았습니다.",'','| HN1 좌표 | HN1 서열 | 다른 균주 서열별 개수 | 비교 불가 |','|---|---|---|---|']
for m in s['markers']:
 p=m['population'];lines.append(f"| {m['focal_reference_position']} | {p['focal_allele']} | "+', '.join(f'{a}: {n}' for a,n in p['allele_counts'].items())+f" | {p['uncallable_assemblies']} |")
lines+=['','앞서 read 정렬에서 관찰한 네 위치 모두에서 HN1 assembly 서열은 집단에서 흔한 서열입니다. 287486에는 AGG와 ATG라는 다른 서열도 있습니다. 이 변이를 오류로 분류하지 않습니다. 앞선 일부 read의 deletion 관찰은 HN1 assembly 오류를 입증하지 않습니다.','','## 기술별 관찰','', '기술 분류는 공개 등록 및 기존 문헌 보완 메타데이터의 조합입니다. 복합 기술은 하나의 조합으로 유지합니다. HiFi/CLR 같은 하위 유형을 이 요약만으로 추정하지 않습니다.','', '| 좌표 | 기술 조합 | 비교 가능 | 서열별 개수 |','|---|---|---|---|']
for m in s['markers']:
 for g in m['strata']:
  if g['field']=='sequencing_group':lines.append(f"| {m['focal_reference_position']} | {g['value']} | {g['callable_assemblies']} | "+', '.join(f'{a}: {n}' for a,n in g['allele_counts'].items())+' |')
lines+=['','## 해석 범위','','이 결과는 assembly 간 서열의 집단적 배경입니다. 다른 균주의 raw read까지 비교한 결과가 아니므로 특정 플랫폼의 염기 오류율이나 체계적 bias를 측정하지 않습니다. 동일 프로젝트와 근연 균주는 독립 표본이 아니며 계통 보정은 아직 하지 않았습니다. 다수 서열도 정답으로 간주하지 않습니다.','','후보는 단백질 상동성 및 양쪽 주변 유전자 방향·순서로 찾았습니다. 전체 orthogroup 확정과 계통수 검증은 별도 작업입니다. 위치가 모호하거나 정렬에서 염기를 비교할 수 없는 경우는 빈도 분모에서 제외했습니다. 이 파일럿의 네 위치는 기존 HN1 read 관찰에서 선택했으며 전 유전체 변이 분석은 아닙니다.','','조사 순서: 동일 종 집단에서 위치 확인 → 변이·주변 유전자·기술/프로젝트별 분포 확인 → 의심 구간의 raw read 및 독립 기술 검증 → 원인 가설과 미해결 사항 제시. 희귀 변이는 자동 오류가 아니며 원인 확정에는 추가 근거가 필요합니다.','','재현 자료: `assemblies.csv`, `marker_alleles.csv`, `summary.json`, `alignment_manifest.json`; 실행 스크립트 `scripts/compare_locus_cohort.py`.']
(out/'report.md').write_text('\n'.join(lines)+'\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(9,4)); labels=[str(m['focal_reference_position']) for m in s['markers']]
alleles=sorted({a for m in s['markers'] for a in m['population']['allele_counts']});bottom=[0]*len(labels)
for allele in alleles:
 vals=[m['population']['allele_counts'].get(allele,0) for m in s['markers']]
 ax.bar(labels,vals,bottom=bottom,label=allele);bottom=[a+b for a,b in zip(bottom,vals)]
ax.set(ylabel='Other assemblies (callable)',xlabel='HN1 coordinate',title='YJU3 population context: assembly sequences, not error calls')
ax.legend(title='Sequence',bbox_to_anchor=(1.02,1),loc='upper left');fig.tight_layout();fig.savefig(out/'allele_distribution.png',dpi=180);plt.close(fig)
