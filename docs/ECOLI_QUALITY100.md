# 대장균 고품질 100개 패널

기존 4균주는 기능 확인용 파일럿이었다. 2026-11-02 사용자 요청에 따라 NCBI의 전체 Complete RefSeq 후보를 조회하고 고품질 100개 패널로 확장한다. 모든 대장균 assembly를 다운로드하는 계획은 저장 공간 제약으로 사용자가 100개 패널로 변경했다.

## 선정

- Taxon 562 하위 분류, current, Complete Genome, RefSeq, exclude_atypical=true 조건의 전체 후보를 조회한다. 첫 페이지만으로 선정하지 않는다.
- NCBI CheckM completeness ≥99%, contamination ≤1%, protein-coding annotation >0, 명시적인 taxonomy check 문제 없음.
- MG1655는 비교 기준으로 고정한다. 나머지는 contamination 오름차순, completeness 내림차순, gap fraction 오름차순, contig N50 내림차순, accession 순이다.
- BioSample 또는 비어 있지 않은 균주명을 기준으로 중복을 제거한다. 균주명은 대소문자와 양쪽 공백만 정리한다. 이름의 동의어·계통적 clone 중복 제거는 아니다.
- 지표 기반 순위는 절대적인 품질 정답이 아니며, 표본은 종 전체를 대표하는 다양성 패널이 아니다. 기술·project·계통 confounding이 남는다.

## 수집과 재현

` .venv/bin/python -u scripts/prepare_ecoli_quality100.py `

전체 후보 및 선정 근거는 `outputs/ecoli_quality100_20261102/candidate_reports.json`, `selection.json`, `selection.csv`에 기록한다. genome, GFF, CDS, protein, assembly report를 각 assembly 폴더에 다운로드하고 NCBI MD5를 검증한 뒤 SHA256을 기록한다. 체크섬 실패는 성공으로 집계하지 않는다. 6GiB 여유 공간을 확보할 수 없으면 중단하고 파일은 재개에 사용한다.

100개가 모두 검증된 후 활성 세균 cohort를 교체하고 gene map 및 rpoB 1:N 검사를 실행한다. 원래 4균주 cohort와 결과는 `previous_cohort.json`, `pilot4_comparison/`에 보존한다. global yeast orthology 계산은 별도이다.

이 분석은 DNA 상동성·양쪽 주변 유전자·rpoB 코딩 span의 변이 구간 비교다. 원인 진단, 전체 세균 orthogroup, 전체 유전체 변이, 타 균주 raw read 검증은 아직 수행하지 않는다. sequencing technology는 NCBI 메타데이터 그대로 기록하고 미제공은 N/A로 둔다.
