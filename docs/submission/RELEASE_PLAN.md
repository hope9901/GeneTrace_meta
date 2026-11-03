> Update 2026-11-03: fresh nimT execution is implemented and verified. See [FRESH_EXECUTION.md](FRESH_EXECUTION.md) and [EXECUTION_COMPARISON.md](EXECUTION_COMPARISON.md). DIAMOND/OrthoFinder jobs are paused. The final video is generated from actual execution captures plus explanatory cards; older snapshot-only notes below describe the earlier stage.

# 11/4 오전 제출 준비 계획

사용자 개인 목표인 2026-11-04 오전을 내부 마감으로 둔다. 제공 가이드에는 공식 마감 Nov4 21:00 EST(계산상 한국 Nov5 11:00)라고 적혀 있으나 플랫폼의 최신 요건은 사용자가 직접 확인한다. 내부 목표를 늦추지 않는다.

## 사용자 — 지금
1. Meta Model API 키와 정확한 모델 ID를 로컬 .env의 MODEL_API_KEY / META_MODEL에 설정. 키는 채팅·Claude·화면녹화·공개 저장소에 노출하지 않는다. 신용카드/크레딧 설정은 직접 처리.
2. 공개 GitHub 저장소를 준비하고 주소 전달. 현재 프로젝트에는 .git이 없다. 검토한 코드 패키지만 올리며 데이터/도구/venv/키를 통째로 올리지 않는다.
3. 제출 플랫폼에 로그인해 필수 입력·영상 길이·URL 접근성·최종 마감 확인. 3분 데모 영상 필요.
4. 생물학적 주장과 사례 선택 검토. nimT 불일치를 오류 확정으로 말하지 않는다.
5. 오늘 저녁 실제 화면을 녹화하고 접근 가능한 영상 링크 준비. 내일 오전은 새 기능 개발 대신 재생·권한·제출 점검.

## Codex — 코드와 재현 패키지
- 실제 근거를 재배포 가능한 작은 데모 snapshot으로 구성. 저장된 분석과 새 모델 호출을 명확히 구분.
- 다운로드/대규모 orthogroup 완료가 없어도 켤 수 있는 제출용 evidence review demo 완성.
- Meta 키 설정 후 실제 호출·모델명·호출 trace·인용 근거를 확인. 실패는 숨기거나 모델 미사용으로 대체하지 않음.
- 기존 research UI의 세균 read/provenance 근거 연결과 위험한 과장 문구 수정.
- README를 현재 상태로 정리, 실행/검증/한계/데이터 및 외부 도구 attribution, 깨끗한 배포 패키지/secret 제외 검사.
- 장기 분석은 계속 진행하되 최종 orthogroup 완료를 제출 시연의 전제조건으로 만들지 않음.

## Claude — docs/submission/claude/만 편집
CLAUDE_TASK.md를 그대로 전달한다. 영어 소개 50/250단어, 180초 대본/shot list, 과장/누락 QA를 맡긴다. 코드·분석·README는 Codex가 담당하므로 수정 금지.

## 출시 기준
- 실제 Meta 모델 호출과 유효한 evidence citation을 포함한 최소 한 번의 성공 run. 없으면 필수 요건 충족이라고 주장할 수 없음.
- 실제 사례 3개(지역 불일치·read 출처 검토·오염 대조), 재현 가능한 명령/해시와 한계 공개.
- 공개 repo 및 접근 가능한 3분 영상. 심사자가 소규모 demo를 원본 대규모 데이터 없이 실행 가능.
- 제출 설명에서 이번 행사 중 만든 기능, 사용한 외부 도구/데이터와 신규 기여를 구분.
