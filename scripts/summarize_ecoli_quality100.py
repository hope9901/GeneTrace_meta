#!/usr/bin/env python3
"""Summarize the completed quality-selected bacterial cohort and population pilot."""
import csv,json,pathlib,os,collections
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/ecoli_quality100_20261102';comparison=ROOT/'outputs/cohort_comparison_20261102/b3987'
def main():
 selection=json.loads((OUT/'selection.json').read_text());rows=json.loads((OUT/'downloaded.json').read_text());s=json.loads((comparison/'summary.json').read_text());assemblies=json.loads((comparison/'assemblies.json').read_text())
 if len(rows)!=100 or s['total_assemblies']!=100:raise ValueError('Expected 100 verified assemblies and 100 locus searches')
 columns=['assembly_accession','strain','biosample','bioproject','assembly_level','genome_size_bp','gc_percent','number_of_proteins','sequence_count','sequencing_group','assembly_method','annotation_provider','completeness','contamination']
 with (OUT/'genomes.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
  for row in sorted(rows,key=lambda r:r['assembly_accession']):
   flat={k:row.get(k,'N/A') for k in columns};flat.update(completeness=row['checkm_info']['completeness'],contamination=row['checkm_info']['contamination']);writer.writerow(flat)
 c=[r['checkm_info']['completeness'] for r in rows];cont=[r['checkm_info']['contamination'] for r in rows]
 summary=dict(verified_assemblies=len(rows),source_files_verified=sum(len(r['source_files']) for r in rows),candidate_count=selection['candidate_count'],eligible_count=selection['eligible_count'],completeness_range=[min(c),max(c)],contamination_range=[min(cont),max(cont)],rpoB_locus_status_counts=s['locus_status_counts'],variable_reference_blocks=len(s['markers']),raw_reads_examined=False,causal_error_diagnosis=False)
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 lines=['# 대장균 고품질 100개 패널','',f"NCBI Complete Genome·RefSeq·exclude_atypical 조건의 {selection['candidate_count']:,}개 전체를 조회했습니다. CheckM completeness ≥99%, contamination ≤1%, 구조 주석의 protein-coding count >0 등 조건을 만족한 {selection['eligible_count']:,}개에서 MG1655 기준 및 순위에 따른 99개를 선정했습니다.",'','GenBank/RefSeq 짝은 RefSeq로 제한했고 BioSample 또는 명시된 동일 균주 이름을 중복 제거했습니다. 이름이 없는 표본은 accession 단위로 유지했습니다. 이것은 유전적 clone 중복 제거가 아니며 전체 종을 대표하는 무작위 표본도 아닙니다.','','순위: contamination 오름차순 → completeness 내림차순 → gap 비율 오름차순 → contig N50 내림차순 → accession. MG1655는 고정 비교 기준으로 포함했습니다. 이것은 명시한 proxy에 따른 순위이며 절대적인 품질 정답은 아닙니다.','',f"100개 genome·GFF·CDS·protein·assembly report, 총 {summary['source_files_verified']}개 파일을 NCBI MD5와 비교했고 SHA256을 기록했습니다. completeness 범위 {min(c):.2f}–{max(c):.2f}%, contamination 범위 {min(cont):.2f}–{max(cont):.2f}%입니다.",'','## rpoB 1:N 결과','',f"MG1655를 기준으로 다른 99개를 비교했습니다. 위치 상태: {json.dumps(s['locus_status_counts'],ensure_ascii=False)}. 코딩 span에서 발견된 reference 반복 단위 변이 구간은 {len(s['markers'])}개입니다. 이것은 정규화된 VCF 변이 개수가 아닙니다.",'','| MG1655 좌표 | MG1655 서열 | 다른 99개 중 비교 가능 | 비교 불가 | 다른 균주의 서열별 개수 |','|---|---|---|---|---|']
 for m in s['markers']:
  p=m['population'];lines.append(f"| {m['focal_reference_position']} | {p['focal_allele']} | {p['callable_assemblies']} | {p['uncallable_assemblies']} | "+', '.join(f'{a or "(deletion)"}: {n}' for a,n in p['allele_counts'].items())+' |')
 lines+=['','## 해석 제한','','집단에서 드문 서열도 자연 변이일 수 있으며 다수 서열도 정답은 아닙니다. 위치·주변 유전자·DNA 서열을 비교했으며 원인에 대한 오류 진단을 내리지 않았습니다. raw read 검증·독립 기술 검증·계통 보정·전체 orthogroup 확정은 수행하지 않았습니다. 고품질 Complete Genome 선택으로 contig assembly에서의 검출 성능은 평가하지 못합니다. 기술 정보는 공개 자료에 없으면 N/A입니다. 대장균 전체 유전자/전 유전체 변이를 분석한 결과가 아니라 rpoB의 1:N 파일럿입니다.','','선정표: selection.csv. 전체 후보 메타데이터: candidate_reports.json. 다운로드 목록: downloaded.json. 상세 서열 결과: ../cohort_comparison_20261102/b3987/marker_alleles.csv. 이전 4균주 결과는 pilot4_comparison/에 보존했습니다.']
 (OUT/'report.md').write_text('\n'.join(lines)+'\n');(comparison/'report.md').write_text('\n'.join(lines)+'\n')
 os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'runs/matplotlib-cache'))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 markers=s['markers'];fig,ax=plt.subplots(figsize=(14,5));x=range(len(markers));same=[m['population']['focal_allele_count'] for m in markers];other=[m['population']['callable_assemblies']-n for m,n in zip(markers,same)];missing=[m['population']['uncallable_assemblies'] for m in markers]
 ax.bar(x,same,label='Same sequence as MG1655',color='#7fa89a');ax.bar(x,other,bottom=same,label='Different sequence (not an error call)',color='#dc9a61');ax.bar(x,missing,bottom=[a+b for a,b in zip(same,other)],label='Uncallable',color='#cccccc')
 if len(markers)<=40:ax.set_xticks(list(x),[str(m['focal_reference_position']) for m in markers],rotation=70,ha='right')
 else:ax.set_xlabel('Variable reference repeat blocks, ordered by MG1655 coordinate')
 ax.set(ylabel='Other assemblies (n = 99)',title='Quality-selected E. coli panel: rpoB population context');ax.legend(loc='upper center',bbox_to_anchor=(.5,1.2),ncol=3);fig.tight_layout();fig.savefig(OUT/'rpoB_population.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
