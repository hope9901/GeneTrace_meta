> Update 2026-11-03: fresh nimT execution is implemented and verified. See [FRESH_EXECUTION.md](FRESH_EXECUTION.md) and [EXECUTION_COMPARISON.md](EXECUTION_COMPARISON.md). DIAMOND/OrthoFinder jobs are paused. The final video is generated from actual execution captures plus explanatory cards; older snapshot-only notes below describe the earlier stage.

# 제출용 사실표 — 2026-11-03

## 제품/대상
GeneTrace: evidence-first research agent for investigating genome reliability before biological interpretation. Model organism cohorts: S. cerevisiae and E. coli. Original contribution is the hackathon-built investigation pipeline, data audit, evidence-linked agent and UI. Existing third-party tools and public data must be credited separately.

## 데이터와 완료된 분석
- 효모: GFF/CDS/protein 제공 435개에서 운영상 strain/BioSample 중복 정리 후 279개. sequence-verified clone deduplication은 아니다.
- 대장균: Complete RefSeq·exclude-atypical 후보 5,512개 전체 평가. CheckM completeness ≥99%, contamination ≤1% 등 필터로 2,161개 적격. MG1655 기준 포함, 나머지 99개를 명시한 품질 순위와 중복 정리로 선정. 종 전체 대표성/전 세계 절대 최상 품질을 주장하지 않는다.
- 100개 genome/GFF/CDS/protein/assembly report 500개 파일 publisher MD5 검증. 단백질 487,647 records. coding gene ID별 최장 모델 487,637개가 정식 OrthoFinder 입력.
- 효모 YJU3 279개 1:N: 단일 위치·양쪽 주변 유전자 지지 277개. 네 위치에서 비교 가능한 다른 균주 각 268개. 나머지 비교 불가는 빈도 분모에서 제외.
- 대장균 rpoB 100개: 위치·양쪽 주변 유전자 지지 100개, 변이 reference repeat block 50개. 정규화 VCF 변이 개수가 아니다.
- 대장균 100개 전체 단백질 Pfam 6-profile 검색 완료. 전 Pfam도, orthogroup도 아니다.
- MMseqs 90% AA identity/80% 양방향 coverage 빠른 상동성 그룹: 9,233개, single-copy 공통 후보 3,027개, 높은 빈도 그룹 누락 조합 915건. 확정 orthogroup이 아니다.
- 30개 공통 single-copy 후보, 7,664 amino-acid alignment columns로 잠정 계통수. 재조합 보정/완전한 orthology reconciliation 미실행.
- 편의 표본 9건 DNA·synteny·기술/BioProject·가까운 표본 감사 완료.

## 실제 read 검증
- 100개 중 63개 exact BioSample에 WGS run 220개. assembly-version에 실제 사용된 run인지 미확정.
- GCF_001558995.2 / JEONG-1266 / SRR2762258: nimT 평균 depth 약31.37, local substitution mismatch 약3.09%, 1bp deletion disagreement 21/23 alignment records(15 forward/6 reverse). 주변 nimR/yeaO도 약2% 불일치. 국소 assembly 오류 확정이 아니라 sample/version/mixture 등 검증이 필요한 사례. Paired records는 독립 분자 수가 아니다.
- GCF_013112335.1 / GF3-3 / SRR11451201: 네 구간 depth25, mismatch0, same-strand start step10bp가 정형적이며 template length450 고정·거의 단일 quality 문자. provenance review required=true. ENA Illumina 등록. 합성이라고 확정하지 않으며 독립 physical sequencing 검증 성공으로 집계하지 않는다.
- HN1 full PacBio run ERR4194695 MD5 검증·정렬 완료. 품질 문자열 !만으로 저품질 sequencing이라고 결론 내리지 않는다. ZMW 독립성·정확한 run-version 사용 미확정.

## 외부 오염 양성 검증
NCBI contaminated 경고 + FCS 보고서 확인 후 3개 사례. detector는 100개 패널 대응률(2kb 이상 contig, 1kb 이상/90% identity 정렬, query coverage<50%)을 먼저 계산하고 FCS EXCLUDE/TRIM/FIX contig와 대조. FCS 좌표를 detector 입력으로 사용하지 않았다. 고오염 편의 표본이며 contig 단위 concordance다.
| Assembly | FCS 조치 전체 | 2kb 이상 검사 대상 | 그중 검출 | FCS 미조치 추가 표시 |
|---|---:|---:|---:|---:|
| GCA_002933585.1 | 4 | 4 | 4 | 0 |
| GCA_027889185.1 | 2922 | 1 | 1 | 31 |
| GCA_009888865.1 | 169 | 126 | 126 | 63 |
짧은 contig를 제외했다. 추가 표시가 많아 specificity가 확립되지 않았다. FCS 미조치는 확정 음성이 아니다. 근연 세균 오염/정상 accessory DNA/모바일 요소 감별은 추가 작업이다. 일부 FCS auxiliary report는 publisher MD5 미제공으로 SHA256 snapshot만 기록했다.

## 진행 중 / 미검증
- 효모279개 정식 OrthoFinder: core inference 진행. core64+additional215; 전체 유전자 포함. 최종 결과 미완료.
- 대장균100개 정식 OrthoFinder: core16 완료, additional84 assignment 진행. 최종 accounting 미완료.
- Meta: 2026-11-03 muse-spark-1.3 실제 nimT run 성공. API 응답 4회, 근거 4개, 인용 해석 4개. Provider 반환 model ID와 인용 ID 유효성 확인. demo/live_meta_verification.json 및 demo/verified_meta_report.json이 검증 기록이다. 기본 모드는 여전히 모델 미사용이며 저장된 보고서 열기는 새 호출이 아니다.
- Hosted live product/public repository/video URL은 아직 확보하지 않았다.
- 기존 80개 테스트 통과는 당시 개발 검증이며 real-world diagnosis accuracy가 아니다. 제출 직전 최신 테스트/skip 수를 별도 기록한다.

증거 원본은 outputs/ecoli_followup_20261102/report.md 및 각 analysis manifest에 있다. 이 표의 사실이 새 결과로 바뀌면 제출 직전에 갱신한다.
