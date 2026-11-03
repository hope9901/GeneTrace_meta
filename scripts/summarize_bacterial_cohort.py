#!/usr/bin/env python3
"""Report the four-strain E. coli population pilot, without error inference."""
import csv,json,pathlib,os
ROOT=pathlib.Path(__file__).resolve().parents[1];out=ROOT/'outputs/cohort_comparison_20261102/b3987'
s=json.loads((out/'summary.json').read_text());assemblies=json.loads((out/'assemblies.json').read_text());rows=list(csv.DictReader((out/'marker_alleles.csv').open()))
lines=['# 大腸菌 / E. coli rpoB 1:N 비교','', 'MG1655를 대상으로 W3110, REL606, BL21(DE3)를 비교했습니다. 네 유전체 모두 NCBI Complete Genome이며 로컬 genome/GFF/CDS/protein 자료를 사용했습니다. 신규 다운로드는 하지 않았습니다.','','| 균주 | Assembly | rpoB 방향 | 단백질 identity / coverage | 겹치는 CDS |','|---|---|---|---|---|']
for a in assemblies:lines.append(f"| {a['strain']} | {a['assembly_accession']} | {a['best_hit']['strand']} | {a['best_hit']['pident']}% / {a['best_hit']['query_coverage']*100:.0f}% | {a['overlapping_cds']} |")
lines+=['',f"4/4에서 단일 후보와 양쪽 주변 유전자 배치가 지지됐습니다. 비교된 단백질은 동일하지만 DNA 정렬에서 {len(s['markers'])}개 변이 구간이 검출됐습니다. reference 반복 구간별로 묶었으므로 이것을 일반적인 VCF variant 개수로 해석하지 않습니다.",'','| MG1655 좌표 | MG1655 서열 | 다른 세 균주 서열별 개수 | 비교 불가 |','|---|---|---|---|']
for m in s['markers']:
 p=m['population'];lines.append(f"| {m['focal_reference_position']} | {p['focal_allele']} | "+', '.join(f'{a}: {n}' for a,n in p['allele_counts'].items())+f" | {p['uncallable_assemblies']} |")
lines+=['','## 판단과 제한','','이 비교에서는 유전자 소실 또는 assembly 오류를 뒷받침하는 근거를 찾지 못했습니다. 서열 차이는 균주 간 변이 후보이며 오류로 분류하지 않습니다. 단백질 동일성도 기능의 동일성을 입증하지 않습니다.','','표본은 네 균주뿐이며 그중 MG1655와 W3110은 K-12 계열입니다. 독립 표본의 모집단 빈도를 추정할 수 없습니다. 전체 orthogroup 확정 및 계통 보정은 하지 않았습니다. 모든 자료가 Complete Genome이므로 contig assembly의 문제 검출 성능도 평가하지 못했습니다.','','세균의 sequencing 기술은 현재 자료에서 조사되지 않아 N/A로 유지했습니다. 기술별 비교와 raw read 검증은 수행하지 않았으며 플랫폼 오류율도 추정하지 않습니다. SNP, indel의 자연 변이 여부 및 원인 확정에는 근연 균주 확대와 독립 read 근거가 필요합니다.','','변이 검색 범위는 상동성으로 찾은 rpoB 코딩 span(단백질 정렬이므로 종결 코돈 제외)이며 주변 ±500bp는 정렬 context에만 사용했습니다. 모든 비교 가능한 차이를 검색하고 reference homopolymer 단위로 묶었습니다. 정렬 의존적이고 VCF 정규화·haplotype 분석은 아닙니다.','','재현: `python3 scripts/compare_locus_cohort.py --assembly GCF_000005845.2 --gene rpoB --workers 2`','', '입력별 해시·명령은 실행 폴더의 result.json과 alignment_manifest.json, 균주별 상세 결과는 assemblies.csv 및 marker_alleles.csv에 있습니다.']
(out/'report.md').write_text('\n'.join(lines)+'\n')
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runs/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
focal=s['focal_assembly'];accs=[a['assembly_accession'] for a in assemblies];pos=[m['focal_reference_position'] for m in s['markers']];lookup={(r['assembly_accession'],int(r['position'])):r['allele'] for r in rows};matrix=[];texts=[]
for acc in accs:
 vals=[];txt=[]
 for m in s['markers']:
  allele=m['population']['focal_allele'] if acc==focal else lookup.get((acc,m['focal_reference_position']))
  vals.append(-1 if allele is None else int(allele!=m['population']['focal_allele']));txt.append('N/A' if allele is None else allele or '-')
 matrix.append(vals);texts.append(txt)
fig,ax=plt.subplots(figsize=(12,3.5));ax.imshow(np.array(matrix),vmin=-1,vmax=1,cmap='coolwarm',aspect='auto')
for i in range(len(accs)):
 for j in range(len(pos)):ax.text(j,i,texts[i][j],ha='center',va='center',fontsize=9)
ax.set_xticks(range(len(pos)),[str(p) for p in pos],rotation=55,ha='right');ax.set_yticks(range(len(accs)),[a['strain'] for a in assemblies]);ax.set_title('E. coli rpoB: sequence blocks (red = differs from MG1655, not an error)');ax.set_xlabel('MG1655 coordinate');fig.tight_layout();fig.savefig(out/'strain_alleles.png',dpi=170);plt.close(fig)
