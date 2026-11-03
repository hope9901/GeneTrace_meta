> Update 2026-11-03: fresh nimT execution is implemented and verified. See [FRESH_EXECUTION.md](FRESH_EXECUTION.md) and [EXECUTION_COMPARISON.md](EXECUTION_COMPARISON.md). DIAMOND/OrthoFinder jobs are paused. The final video is generated from actual execution captures plus explanatory cards; older snapshot-only notes below describe the earlier stage.

# 提出前 / 제출 전 점검

## 제출을 막는 필수 항목
- [ ] Meta 키/정확한 모델 ID 설정, 실제 Meta run 성공, evidence citation과 호출 횟수·모델명 확인.
- [ ] 공개 저장소 생성·업로드 후 로그아웃 상태에서 접근 확인. dist ZIP 자체를 raw data처럼 올리지 말고 소스 파일로 게시.
- [ ] 소유자가 배포 라이선스/공개 범위 결정. 원본 API 키·raw read·설치된 바이너리는 업로드하지 않기.
- [ ] 3분 영상 제작·접근 가능한 링크 확인. 이번 행사에서 만든 부분과 기존 외부 도구 구분.
- [ ] 제출 플랫폼의 최신 필수 필드/마감 직접 확인. 사용자 내부 목표는 11/4 오전.

## 촬영 순서
1. `python3 -m genetrace.demo --port 8766` 후 http://127.0.0.1:8766/.
2. Meta 설정 상태 확인, 실제 Meta 모드로 nimT 조사. 모델이 검사 도구를 선택하는 trace와 EID 인용 표시.
3. nimT 불일치뿐 아니라 주변 대조 불일치도 보여주고 원인은 미해결이라고 설명.
4. read provenance 사례의 정형적 패턴과 독립 검증 보류를 보여주기.
5. FCS 사례의 eligible denominator와 추가 표시를 함께 보여주기.
6. source SHA256/원시 근거 JSON 내보내기.
7. 필요시 full workspace 8765에서 백그라운드 OrthoFinder 상태를 보여주되 완료 전 결과를 주장하지 않기.

녹화 중 키를 표시하지 않는다. 저장된 분석 조회를 새 alignment 실행으로 소개하지 않는다. API 호출이 실패하면 실패 상태를 숨기지 말고 원인을 해결한 후 다시 녹화한다. 편집 영상에서 기다리는 시간을 줄이더라도 실제 지연시간을 오해하게 하지 않는다.

## 재현 확인
- 깨끗한 압축 해제 디렉터리에서 Python 표준 라이브러리만으로 demo 시작.
- `python3 -m unittest discover -s tests -p test_demo.py -v` 통과.
- `python3 scripts/submission_preflight.py`에서 실제 Meta 호출과 URL 상태 점검.
- docs/submission/claude/의 영어 문구는 FACTS.md와 대조 후 채택.
