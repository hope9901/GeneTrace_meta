# GeneTrace local investigation workspace

## 실행

```sh
python3 -m genetrace.web --port 8765
```

브라우저에서 `http://127.0.0.1:8765`를 엽니다. 실행 중인 orthogroup·read 분석과 별도 프로세스입니다. Mac을 종료하면 서버와 분석 모두 중단되므로 재개 안내를 참고하세요. 외부 서버에 배포하지 않은 로컬 인터페이스입니다.

## 현재 지원 범위

- 279개 대표 assembly와 S288C의 6,002개 단백질 코딩 유전자에서 선택합니다. 기준 유전자명, locus tag, protein accession을 정확히 지정할 수 있습니다.
- 실제 BLASTp 및 TBLASTN 검색을 단일 스레드로 실행합니다. 기준은 identity >=90%, single-HSP query coverage >=90%, E<=1e-20입니다. 최대 100개 subject hit를 요청하므로 hit 수는 유전자 copy number가 아닙니다. intron-rich 또는 divergent homolog는 놓칠 수 있습니다.
- 후보 위치의 GFF feature·CDS overlap, frameshift/pseudogene 기록, 기존 annotation scope audit, N 비율, GC와 contig 경계까지의 거리를 점검합니다.
- raw read 연결 등록과 이용 가능한 지역 분석 결과를 조회합니다. 현재 지역 read 결과 자동 연결은 HN1/YJU3에만 구현했습니다. 다른 대상은 metadata만 반환하며 read 검사를 수행한 것으로 표시하지 않습니다.
- 전체 orthogroup 계산의 상태를 조회합니다. 계산 완료 전 소속·확장·소실을 추정하지 않습니다.
- 검사마다 근거 ID, 명령·도구 버전, 입력 SHA256 및 원시 alignment를 `runs/agent/<id>/`에 저장하고 전체 근거 JSON을 내려받을 수 있습니다.

구간 도식은 BLAST 상동 후보와 겹치거나 가까운 GFF gene feature를 좌표 비례로 표시합니다. 이 도식만으로 synteny 또는 orthology를 입증하지 않습니다. 사용자가 임의의 raw FASTQ를 업로드하거나 지역 read 분석을 새로 시작하는 기능은 아직 없습니다.

## 두 가지 실행 방식

**근거 검사 / 모델 미사용**: 고정된 일곱 검사 순서를 실행합니다. 실제 생물정보 분석을 수행하지만 AI가 계획을 세웠다고 표시하지 않습니다.

**Meta 에이전트**: 실제 서열 검색 후 Meta 모델이 허용된 도구 중 다음 검사를 선택하고, 제공된 evidence ID를 인용한 가설을 제안합니다. 도구는 한 번씩만 실행하고 생물학적 소실·기능 판정은 확정하지 않습니다. 모델이 임의 경로·명령·검사 인자를 만들 수 없습니다. 알려지지 않은 도구, 재실행 요청과 잘못된 근거 인용은 거부합니다. Meta 요청 실패를 성공이나 모델 미사용 모드로 조용히 대체하지 않습니다.

Meta 설정은 아직 없으며 실제 API 호출을 검증하지 않았습니다. 공식 API에 맞춘 연결 코드와 응답 검증을 구현했습니다. 키가 없는 상태에서는 서버가 Meta 조사를 시작하지 않습니다.

```sh
cp .env.example .env
```

`.env`의 `MODEL_API_KEY`와 `META_MODEL`을 직접 채웁니다. 정확한 모델 ID는 Meta 대시보드 또는 `/v1/models`에서 확인하세요. 키는 채팅이나 저장소에 넣지 않습니다. 공식 주소는 `https://api.meta.ai/v1`이며 키는 서버에서만 읽고 화면·근거 JSON에 저장하지 않습니다. 실제 모델 호출 시에는 공개 assembly·유전자 정보와 수집한 요약 근거가 Meta API로 전송됩니다. 전체 raw FASTQ는 전송하지 않습니다.

공식 문서: [API reference](https://dev.meta.ai/docs/api-reference), [Chat completions](https://dev.meta.ai/docs/api-reference/chat-completions).

## 확인한 실제 사례

- HN1 / YJU3: DNA identity 99.042%, query coverage 100%, 단백질 strong hit 0, CDS overlap 없음 → annotation gap candidate. 기능 미평가.
- JXXY16.1 / ACS1: 핵 서열에 비해 제공 주석이 소기관에 제한됨 → annotation scope mismatch. 개별 유전자 소실로 해석하지 않음.
- FostersB / SWC3: frameshift 기록 → disrupted model unresolved. 서열 오류·assembly·실제 변이의 구분에 추가 근거 필요.

단위 테스트 39개와 실제 BLAST 사례를 확인했습니다. Meta live API, 전체 orthogroup 후속 배정, 전체 read에 대한 variant/consensus 평가는 완료로 간주하지 않습니다.

## 안전한 로컬 실행 범위

서버는 127.0.0.1에만 바인딩합니다. API는 동일 origin 요청을 받고 한 번에 조사 한 건만 실행합니다. 외부 공개, 임의 파일 경로 조회, shell 명령 실행 API는 없습니다. 프런트엔드는 모델·데이터 내용을 HTML이 아니라 텍스트로 표시합니다. 메타데이터와 모델 출력은 지시가 아닌 데이터로 다룹니다.


## 2026-11-01 도메인·세균 확장

E. coli 네 모델 균주와 MG1655 기준 코딩 단백질 4,300개의 선택을 추가했습니다. Pfam 6-profile 패널의 전체 proteome 검색 및 세균 6-frame ORF 검색 결과를 연결합니다. `lookup_domain_evidence`는 기준·대상 단백질의 도메인 근거와 검색 범위를 조회하며, 미검출을 유전자 소실로 해석하지 않습니다. 현재 테스트 45개가 통과했습니다. 자세한 방법은 [도메인·세균 보고서](DOMAIN_AND_BACTERIA.md)를 참고하세요.


## 2026-11-02 구간 QC 확장

`lookup_regional_qc`를 추가했습니다. 기존 full BAM의 HN1 5개 구간에 대한 mismatch·indel·clipping·read spanning 근거를 조회합니다. 결과가 없는 대상은 정상으로 표시하지 않습니다. 시뮬레이션 검증은 `/api/validation`과 로컬 화면에서 범위·한계를 함께 표시합니다. [검증 보고서](REGIONAL_VALIDATION.md)를 참고하세요.
