# 2026-10-30 작업 중단 및 재개

사용자의 Mac 종료 요청으로 로컬 orthogroup 계산, FASTQ 다운로드, 완료 후 분석 대기 프로세스를 중단했습니다. 원격 서버에서 실행 중인 작업은 없습니다.

- Core 검색 완료 체크포인트: 986개. OrthoFinder 3.2는 체크포인트와 파일 크기를 확인하여 완료된 검색을 재사용하고 미완료 검색을 다시 실행합니다.
- HN1 FASTQ 완료 조각: 327,471,122 / 927,324,001 bytes. `.chunks` 폴더를 보존하세요.
- 중단 기록: `outputs/shutdown_checkpoint_20261030.json`

Mac을 다시 켠 뒤 프로젝트 폴더에서 명시적으로 재개합니다. 다음 명령은 아직 실행하지 않았습니다.

```sh
python3 scripts/run_scalable_orthology.py --resume --threads 6 --analysis-threads 1
```

FASTQ 다운로드는 다른 터미널에서 재개할 수 있습니다.

```sh
python3 scripts/download_verified_read.py --url 'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/005/ERR4194695/ERR4194695_subreads.fastq.gz' --bytes 927324001 --md5 1d3414c2defd37fd93d3f90342ebb569 --output data/read_dataset_20261030/fastq/ERR4194695_subreads.fastq.gz --workers 4
```

전체 FASTQ 검증 후의 자동 분석 대기도 필요하면 별도 터미널에서 재개합니다.

```sh
.venv/bin/python scripts/finish_hn1_read_qc.py
```

원격 이전은 아직 하지 않았습니다. 원격 장비와 저장 공간·접속 권한을 정한 뒤 입력·체크포인트를 복사하고, Mac 전용 실행 환경은 원격 운영체제에 맞게 다시 설치해야 합니다.
