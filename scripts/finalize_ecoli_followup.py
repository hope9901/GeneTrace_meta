#!/usr/bin/env python3
"""Finalize reports after the already-running raw and FCS workflows terminate."""
import datetime,json,pathlib,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
def load(p):
 try:return json.loads(p.read_text())
 except (FileNotFoundError,json.JSONDecodeError):return {}
def main():
 out=ROOT/'outputs/ecoli_followup_20261102';out.mkdir(exist_ok=True)
 for _ in range(720):
  raw=load(ROOT/'outputs/ecoli100_raw_validation_20261102/status.json');fcs=load(ROOT/'outputs/contamination_validation_20261102/status.json')
  (out/'status.json').write_text(json.dumps(dict(stage='waiting_for_existing_workflows',raw=raw,fcs=fcs),indent=2)+'\n')
  if raw.get('stage') in ['complete','failed_resumable'] and fcs.get('stage') in ['complete','failed_resumable']:break
  time.sleep(10)
 else:raise RuntimeError('Workflow finalization timeout; original checkpoints retained')
 if (ROOT/'outputs/ecoli100_raw_validation_20261102/summary.json').exists():
  subprocess.run([sys.executable,str(ROOT/'scripts/expand_ecoli_raw_gene_regions.py')],check=True)
  subprocess.run([sys.executable,str(ROOT/'scripts/audit_existing_ecoli_read_provenance.py')],check=True)
 raw_summary=load(ROOT/'outputs/ecoli100_raw_validation_20261102/summary.json');fcs_summary=load(ROOT/'outputs/contamination_validation_20261102/summary.json');lines=['# 大腸菌 / 대장균 후속 검증','', '## 전체 유전자와 집단 맥락','','전체 487,647 단백질의 빠른 상동성 그룹화는 완료했다(9,233 groups). 이는 확정 orthogroup이 아니다. 정식 OrthoFinder 100개·487,637 코딩 유전자 분석은 별도 백그라운드 작업이며 최종 accounting 검증 전에는 완료로 처리하지 않는다. 차이 10개는 같은 gene ID에서 최장 단백질 모델 하나를 사용하는 데 따른 것이다.','','30개 core 단백질·7,664 정렬 columns의 잠정 계통수와 가까운 표본·기술/프로젝트 맥락을 만들었다. 9개 후보를 DNA와 양쪽 주변 유전자로 감사했다. 재조합 보정과 완전한 orthology reconciliation은 하지 않았다.','','## Read 검사','']
 for a in raw_summary.get('assemblies',[]):
  review=a.get('read_provenance',{});lines.append(f"- {a['assembly']} / {a['run']}: whole-file MD5 검증. provenance review required={review.get('provenance_review_required')}. 정확한 assembly-version 사용 run은 미확정.")
  for r in a['regions']:lines.append(f"  - {r['name']}: mean depth {r['metrics']['mean_depth']:.3f}, mismatch {r['metrics']['mismatch_fraction']:.6f}, signals={r['signals']}.")
 lines+=['','SRR11451201에서 일정한 시작 간격·고정 template 길이·거의 단일 품질 문자·정확한 매치 패턴을 함께 관찰했다. ENA는 Illumina로 등록했지만 파생/정형화 과정이 미확정이므로 독립적인 실제 sequencing 검증 성공으로 집계하지 않는다. 이 패턴만으로 합성 데이터라고 확정하지 않는다. Paired records는 독립 분자 수가 아니다.','','## NCBI FCS 외부 오염 좌표 대조','','대장균 contaminated 경고와 FCS 보고서를 확인한 양성 후보에서 세 사례를 골랐다. 100개 reference panel 서열에 낮게 대응하는 contig를 미리 정한 기준으로 표시한 후 FCS EXCLUDE/TRIM/FIX 대상 contig와 대조했다. FCS 좌표를 검출 입력으로 사용하지 않았다.','', '| Assembly | FCS 조치 contig 전체 | 2kb 이상 검사 대상 | 그중 검출 | FCS 조치 없는 contig | 그중 추가 표시 |','|---|---|---|---|---|---|']
 for e in fcs_summary.get('evaluations',[]):lines.append(f"| {e['assembly']} | {e['FCS_action_contigs']} | {e.get('FCS_action_contigs_eligible','N/A')} | {e['FCS_action_contigs_flagged']} | {e['non_FCS_action_contigs']} | {e['non_FCS_action_contigs_flagged']} |")
 lines+=['','이는 고오염 양성 사례에서의 contig 단위 concordance다. 정밀한 오염 span 검출·독립적인 분류기 정확도·전체 sequencing/assembly 오류 진단 성능이 아니다. FCS가 조치하지 않은 contig도 확정 음성이 아니며 정상 accessory/mobile DNA가 panel에 없을 수 있다. 가까운 세균 간 오염은 서열 대응 검사로 놓칠 수 있다. 유전체 파일은 NCBI MD5로 검증했다. 일부 FCS 보조 보고서는 제공 MD5가 없어 해당 사실과 내려받은 파일 SHA256을 기록했다.','','공식 기준: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/quality/contamination/fcs-contamination/','','## 현재 상태','',f'Raw workflow: {raw.get("stage")}; FCS workflow: {fcs.get("stage")}. 실패 항목은 각 status.json을 확인한다. OrthoFinder는 별도 status에 기록되며 이 보고서가 orthogroup 완료를 뜻하지 않는다.']
 (out/'report.md').write_text('\n'.join(lines)+'\n');(out/'status.json').write_text(json.dumps(dict(stage='reports_written',raw_stage=raw.get('stage'),fcs_stage=fcs.get('stage'),orthogroup_complete=False,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2)+'\n')
if __name__=='__main__':main()
