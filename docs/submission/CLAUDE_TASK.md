# Claude에 그대로 전달할 작업 지시

프로젝트 경로: /Users/hyg/Desktop/project/meta_hackathon
목표: 2026-11-04 오전 내부 제출 준비 완료. 사용자가 제공한 가이드에는 공식 마감이 11/4 21:00 EST라고 적혀 있지만 변경 여부를 확인하지 않았다. 사용자 목표 시각을 우선한다.

당신의 역할은 제출 문구와 영상 구성이다. Codex는 코드·실험·모델 연결·재현 패키지를 담당한다. 충돌을 막기 위해 **docs/submission/claude/** 아래에만 파일을 작성한다. genetrace/, web/, scripts/, tests/, README.md, 데이터, 실행 중인 프로세스, .env는 변경하거나 읽지 않는다. 프로세스를 중지하거나 추가 분석/대규모 다운로드를 시작하지 않는다. API 키를 요구·복사·출력하지 않는다.

먼저 docs/submission/FACTS.md와 docs/submission/RELEASE_PLAN.md를 읽는다. 필요하면 docs/CONTAMINATION_VALIDATION.md, outputs/ecoli_followup_20261102/report.md, outputs/ecoli100_candidate_audit_20261102/report.md를 읽는다. 수치가 충돌하면 추측하지 말고 QA_NOTES.md에 적는다.

산출물:
1. PITCH_EN.md: 한 줄 소개, 문제/사용자, 제품, Meta의 역할, 기술적 차별점, 검증과 한계. 짧은 소개 50단어와 상세 소개 250단어 두 버전.
2. VIDEO_180S.md: 영어 내레이션 320–380단어. 0:00–0:25 문제, 0:25–0:50 에이전트 조사, 0:50–1:35 실제 nimT read 불일치 및 주변 대조, 1:35–2:05 정형적 read 출처 경고, 2:05–2:30 NCBI FCS 외부 검증과 추가 표시, 2:30–3:00 가치·Meta·한계. 화면/클릭/내레이션을 표로 작성.
3. QA_NOTES.md: 과장·수치 혼동·불명확한 Meta 기여·촬영 누락을 우선순위별 지적.

제품 메시지: GeneTrace는 생물학자와 계산 분석 사이를 연결하여, 후속 연구 전에 의심 유전체 구간의 신뢰성을 근거로 조사한다. 전 유전체 통계만으로 오류를 판정하지 않는다. minority allele나 indel을 자동 오류로 분류하지 않는다. evidence ID, 실제 명령/버전/해시, 불확실성과 다음 검사를 보여준다.

금지 표현: 100% 정확도, 모든 오류 탐지, 진짜 소실 확정, FCS로 자체 정확도 입증, 단백질 상동 그룹을 확정 orthogroup으로 부르기, Meta 연결 코드만 있는 상태를 실제 호출 성공이라고 쓰기, 저장된 분석 결과 조회를 방금 새 분석을 실행했다고 표현하기. Meta 실제 호출 완료 여부는 사용자가 검증 결과를 전달하기 전까지 [PENDING LIVE META VERIFICATION]으로 둔다.

제공 데이터는 공개 NCBI/ENA 자료다. demo에서는 저장된 실제 분석과 새 모델 호출을 명확히 구분해야 한다. 결과 문서 작성 후 사용자에게 세 파일 경로와 수정이 필요한 주장만 간단히 알려라. 외부 게시/업로드/저장소 push는 하지 않는다.
