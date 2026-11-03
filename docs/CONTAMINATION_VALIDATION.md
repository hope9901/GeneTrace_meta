# NCBI contaminated 데이터 기반 검증

NCBI FCS summary와 assembly metadata의 `atypical.warnings: contaminated`를 확인한 대장균 양성 후보 100개를 확보했다. 이 중 3개 genome 및 FCS 좌표 보고서를 가져와 비교했다. 의도적으로 고오염 양성 자료를 택했으므로 모집단 표본이나 전체 진단 정확도 benchmark가 아니다.

검출 입력은 고품질 대장균 100개 패널과 테스트 genome 서열이다. 완전히 같은 reference contig는 SHA256으로 묶되 100개 모두의 서열을 포함했다. minimap2 asm20에서 길이 1kb 이상·identity 90% 이상 정렬의 query coverage를 합집합으로 계산했다. 2kb 이상 contig 중 대응 비율 50% 미만인 contig를 후보로 표시했다. 이 기준은 FCS 좌표와 대조하기 전에 고정했으며 FCS 좌표는 detector 입력으로 쓰지 않았다.

그 후 FCS EXCLUDE/TRIM/FIX 표시 contig와 대조했다. contig 전체의 대응성 검사이며 FCS span 경계를 복원한 것은 아니다. 짧은 contig는 현재 검사 범위 밖이다.

| Assembly | FCS 조치 contig 전체 | 2kb 이상 검사 대상 | 그중 표시 | FCS 조치 없는 contig의 추가 표시 |
|---|---|---|---|---|
| GCA_002933585.1 | 4 | 4 | 4 | 0 |
| GCA_027889185.1 | 2922 | 1 | 1 | 31 |
| GCA_009888865.1 | 169 | 126 | 126 | 63 |

검사 가능한 FCS 조치 contig에서는 일치했지만, 두 번째·세 번째 genome에서 추가 표시가 많았다. 이것을 모두 오염이나 확정 false positive로 부르지 않는다. 정상 accessory DNA·plasmid·divergence와 FCS 미표시 오염 등이 남을 수 있다. FCS 조치가 없다는 사실도 확정 음성을 뜻하지 않는다. 이 결과로 독립적인 contamination classifier의 정확도, 자연 변이 감별력, sequencing/assembly 오류 진단 성능을 주장하지 않는다.

다음 개선은 짧은 contig에 대한 별도 분석, mobile element 및 정상 accessory DNA 구분, 세균 간 근연 오염을 위한 독립적인 분류 근거 추가, 정확한 span 단위 평가 및 별도 음성/보류 검증군이다. thresholds를 바꿀 경우 현재 결과와 분리한 탐색적 분석으로 기록해야 한다.

원본 genome은 NCBI 제공 MD5와 대조했다. 일부 FCS 보조 보고서는 publisher MD5 목록에 없어 그 사실과 내려받은 snapshot의 SHA256을 기록했다. 공식 label은 오염 외부 근거이며 단일 염기·indel 오류의 정답으로 사용하지 않는다.

공식 출처:
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/quality/contamination/fcs-contamination/
- https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_fcs_summary_and_details.txt

실행 결과: outputs/contamination_validation_20261102/
