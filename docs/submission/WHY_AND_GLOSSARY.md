# 비전공자를 위한 프로젝트 설명과 영상 도입부

## 한 문장

유전자가 없어 보인다는 결과를 믿고 실험을 시작하기 전에, 진짜 생물학적 차이인지 데이터가 만들어지는 과정의 문제인지 근거를 확인하는 연구 보조 에이전트입니다.

English: **Before researchers act on an apparent gene difference, GeneTrace helps them inspect whether the evidence supports biology, a data problem, or an unresolved question.**

## 왜 필요한가

DNA sequencing 결과는 바로 완성된 유전자 지도가 아닙니다. DNA 조각을 읽고(reads), 이어 붙여 유전체를 구성하고(assembly), 유전자 위치와 구조를 표시하는(annotation) 단계를 거칩니다. 실제 생물학적 차이뿐 아니라 이 과정의 문제도 유전자가 없어지거나 달라진 것처럼 보이게 할 수 있습니다.

이 차이를 확인하려면 생물학적 맥락과 계산 분석을 함께 봐야 합니다. GeneTrace는 여러 균주, 주변 유전자, read와 assembly의 불일치, 데이터 출처를 연결해 연구자가 다음에 무엇을 검증해야 하는지 돕습니다. 시간 절감이나 진단 정확도를 정량적으로 입증한 제품이라는 주장은 하지 않습니다.

## 25–35초 영어 도입부 초안

> A gene appears to be missing. Is that a real biological difference, or a problem in the data? Sequencing gives us DNA fragments. Software reconstructs a genome, then marks where genes are expected. A missing label does not prove the DNA is gone. Before researchers plan their next experiment, GeneTrace brings population comparisons, gene annotations and sequencing evidence together—and shows what still needs checking.

말하는 속도에 따라 시간을 직접 재고 조정합니다. 기존 180초 영상의 첫 도입부를 대체하는 용도이며 총 길이에 추가하지 않습니다.

## 비유

책의 조각을 모아 원문을 복원한 다음 장과 절을 표시한다고 생각하면 됩니다. 장 제목 표시가 없다고 해서 그 내용까지 없어진 것은 아닙니다. Reads=문장 조각, Assembly=복원한 텍스트, Annotation=장·절 표시라는 개념적 비유입니다. 유전자와 책의 장이 생물학적으로 일대일 대응한다는 뜻은 아닙니다.

## 화면 용어 빠른 안내

| 용어 | 쉬운 설명 | 해석 주의 |
|---|---|---|
| Genome / gene / protein | 유전체는 DNA 전체, 유전자는 기능성 RNA·단백질 생성에 기여하는 DNA 구간, 단백질은 세포의 여러 작업을 수행하는 분자 | 모든 유전자가 단백질을 만드는 것은 아님 |
| Read | sequencing 데이터에 담긴 DNA 서열 조각 | 등록된 read도 생성·가공 경로 확인 필요 |
| Assembly / contig | 조각으로 복원한 유전체 / 연속된 복원 서열 한 덩어리 | 완성된 assembly도 무오류 보장은 아님 |
| Annotation / CDS | 유전자 위치·구조 표시 / 단백질 코딩 구간 표시 | 표시가 없다고 DNA가 없는 것은 아님 |
| Strain / reference / 1:N | 같은 종 내 계통 / 비교 기준 / 하나를 여러 관련 유전체와 비교 | 기준·다수 서열도 자동 정답 아님 |
| Allele | 같은 위치의 다른 DNA 서열 버전 | 소수·희귀 allele도 자연 변이일 수 있음 |
| Depth / mismatch | 위치를 덮는 aligned read 수 / read와 assembly의 염기 불일치 | 높은 depth나 불일치만으로 오류 확정 불가 |
| Indel / frameshift | 상대적 삽입·결실 / 코딩 서열을 읽는 묶음의 변화 | 코딩 구간에서 길이가 3의 배수가 아닌 indel은 reading frame을 바꿀 수 있음. 모든 indel이 오류는 아님 |
| Identity / query coverage | 정렬된 부분의 일치 비율 / 검색한 서열 중 정렬된 비율 | 짧은 구간의 높은 일치율과 전체 유전자 일치는 다름 |
| Synteny / flanking genes | 유전자 배치의 보존 / 양옆 유전자 | 같은 위치를 비교하는지 확인하는 보조 근거 |
| Homolog / orthogroup / domain | 공통 조상 유래 관계 / 한 조상 유전자에서 유래한 집합 / 단백질의 인식 가능한 영역 | 유사성 그룹·공유 도메인이 orthology나 기능을 입증하지 않음 |
| Provenance / contamination / FCS | 생성·가공 출처 / 의도하지 않은 서열 혼입 / NCBI 외부 오염 검사 | FCS는 모든 sequencing·assembly 오류의 정답이 아님 |

## 3분 영상에 넣을 최소 개념

1. Reads → Assembly → Annotation을 한 번만 설명합니다.
2. nimT 사례에서는 “이 위치의 read와 완성 서열이 다르다. 주변도 달라서 원인 확인이 더 필요하다”라고 말합니다.
3. Provenance 사례에서는 “완벽한 일치도 독립적인 실제 측정이라는 보장은 아니다”라고 말합니다.
4. FCS 사례에서는 “알려진 오염 자료와 비교했지만 짧은 조각과 추가 표시라는 한계가 있다”라고 말합니다.
5. Orthogroup, synteny, HMM, 모델 token 가격 등은 질문을 받으면 설명하고 도입부를 채우지 않습니다.

일반 정의 참고: [NHGRI glossary](https://www.genome.gov/genetics-glossary), [NCBI annotation guide](https://support.nlm.nih.gov/knowledgebase/article/KA-03574/en-us). 제품 기능·결과는 프로젝트의 실제 분석 기록을 기준으로 합니다.
