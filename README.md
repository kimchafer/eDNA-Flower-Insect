# eDNA
eDNA analyses for insect pollinators

## Reference DB Setup

분석에 앞서 편의를 위해 작업 디렉토리 구조를 설정하자.

DB와 분석용 디렉토리를 우선 별개로 생성한다. ('.'으로 표현된 베이스 디렉토리 이름은 'eDNA'라고 가정하자)
```
.
├── analyses
└── DB
```
esearch 기능의 사용 불가로 인해 gb파일로부터 각 서열에 상응하는 accesstion number와 taxonomy ID를 추출해야만 했다.

fasta파일의 헤더에는 taxonomy ID 정보가 없기 때문에 gb형식을 사용해야 한다.

Cytochrome c oxidase subunit I의 축약 형태가 다양하므로 모든 키워드 검색 결과를 개별 gb파일로 다운로드하자.
```
.
├── analysis
└── DB
    ├── Co1_Insecta.gb
    ├── CoI_Insecta.gb
    ├── Cox1_Insecta.gb
    └── CoxI_Insecta.gb
```
DB 디렉토리에 다운로드된 파일들을 모으고, 다음의 파이썬 스크립트로 분석을 수행하면 하나의 파일로 합칠 수 있다: `1_combine_gb_files.py`

모든 작업 스크립트는 작업 디렉토리에 저장된 상태임을 전제로 한다.
```
import os

# 입력 폴더와 출력 파일 경로 설정
input_folder = "/eDNA/DB/"  # 여러 .gb 파일이 저장된 폴더 경로
output_file = "COI_all_Insecta.gb"         # 합쳐진 파일의 출력 경로

# 파일 합치기
with open(output_file, 'w') as outfile:
    for filename in os.listdir(input_folder):
        if filename.endswith(".gb"):  # .gb 확장자 확인
            file_path = os.path.join(input_folder, filename)
            with open(file_path, 'r') as infile:
                content = infile.read()
                # 각 파일 내용을 추가
                outfile.write(content)
                # 구분자 추가 (중복 방지용으로 "//" 추가 확인)
                if not content.strip().endswith("//"):
                    outfile.write("\n//\n")

print(f"모든 .gb 파일이 {output_file}로 합쳐졌습니다.")
```
작업이 완료된 작업 디렉토리 구성은 다음과 같다:
```
.
├── analysis
└── DB
    ├── 1_combine_gb_files.py
    ├── Co1_Insecta.gb
    ├── CoI_Insecta.gb
    ├── Cox1_Insecta.gb
    ├── CoxI_Insecta.gb
    └── COI_all_Insecta.gb
```

'COI_all_Insecta.gb'가 만들어지면 성공이다. 하나로 합친 파일 외에 나머지 개별 파일들은 삭제하여도 좋다.

이제 두 가지 정보를 추출해야 한다.

gb파일로부터 추출한 accesstion number와 taxonomy ID는 다음과 같이 정렬되어야 한다:
```
#taxidmap.txt 예시
KX423731.2	1987132
LC797984.1	2726158
LC797983.1	2726158
LC797982.1	2726158
LC797981.1	2726158
LC797980.1	2726158
LC797979.1	2726158
LC797978.1	2726158
LC797977.1	2726158
LC797976.1	2726158
```

이를 위해 `2_taxidmap.py` 파이썬 스크립트를 실행시켜보자.

Output file로 'taxidmap.txt'가 생성되면 추출이 잘 완료된 것이다.
```
import re

# 입력 및 출력 파일 경로 명시
input_file = "COI_all_Insecta.gb"  # 입력 파일 (GenBank 파일)
output_file = "taxidmap.txt"       # 출력 파일 (Taxonomy ID 매핑 결과 저장)

# 추출 결과를 저장할 리스트
results = []

# 파일 읽기
with open(input_file, 'r') as f:
    lines = f.readlines()

# VERSION과 Taxon ID 추출
version = None
taxon_found = False  # 현재 `VERSION`에 대해 taxon ID가 발견되었는지 추적
for line in lines:
    # VERSION 추출
    if line.startswith("VERSION"):
        # 이전 `VERSION` 처리
        if version is not None and not taxon_found:
            results.append(f"{version}\t0")  # taxon ID가 없는 경우 0으로 처리

        # 현재 `VERSION` 업데이트
        version = line.strip().split()[1]  # VERSION의 두 번째 요소
        taxon_found = False  # 새로운 VERSION의 taxon ID 여부 초기화

    # taxon ID 추출
    elif "/db_xref=\"taxon:" in line and version is not None:
        if not taxon_found:  # taxon ID가 발견되지 않은 경우만 추가
            taxon_id = re.search(r'taxon:(\d+)', line).group(1)
            results.append(f"{version}\t{taxon_id}")
            taxon_found = True  # taxon ID가 발견됨

# 마지막 VERSION 처리
if version is not None and not taxon_found:
    results.append(f"{version}\t0")

# 파일 쓰기
with open(output_file, 'w') as f:
    for result in results:
        f.write(result + '\n')

print(f"VERSION과 Taxon ID가 {output_file}에 저장되었습니다.")
```
작업이 정상적으로 완료되었다면 작업 디렉토리 내 파일 구성은 다음과 같다:
```
.
├── analysis
└── DB
    ├── 1_combine_gb_files.py
    ├── 2_taxidmap.py
    ├── COI_all_Insecta.gb
    └── taxidmap.txt
```

이제 레퍼런스 DB 제작을 위해 gb 파일을 fasta 형식으로 전환한다: `3_convert_gb_to_fasta.py`

파일 용량에 따라 시간이 꽤 소요될 수 있다.
```
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

# 입력 및 출력 파일 경로 설정
genbank_file = "/eDNA/DB/COI_all_Insecta.gb"
output_fasta_file = "/eDNA/DB/COI_all_Insecta.fasta"
skipped_ids_file = "eDNA/DB/convert_skipped_ids.txt"

# 통계 변수 초기화
total_records = 0
converted_records = 0
skipped_records = 0

# 스킵된 ID 기록용
skipped_ids = []

# GenBank에서 FASTA로 변환
with open(genbank_file, "r") as gb_file, \
     open(output_fasta_file, "w") as fasta_out, \
     open(skipped_ids_file, "w") as skipped_out:

    for record in SeqIO.parse(gb_file, "genbank"):
        total_records += 1
        try:
            # 서열 데이터가 정의되지 않은 경우 예외 처리
            sequence = str(record.seq)
            
            # FASTA 헤더: >ACCESSION VERSION DESCRIPTION
            header = f">{record.id} {record.description}"
            fasta_out.write(f"{header}\n{sequence}\n")
            converted_records += 1

        except (UndefinedSequenceError, ValueError) as e:
            # 에러가 발생한 경우 스킵
            skipped_ids.append(record.id)
            skipped_out.write(f"{record.id}\n")
            skipped_records += 1

# 통계 출력
print(f"Total records processed: {total_records}")
print(f"Successfully converted records: {converted_records}")
print(f"Skipped records: {skipped_records}")
print(f"Skipped IDs saved to: {skipped_ids_file}")
```

작업이 정상적으로 완료되었다면 작업 디렉토리 내 파일 구성은 다음과 같다:
```
.
├── analysis
└── DB
    ├── 1_combine_gb_files.py
    ├── 2_taxidmap.py
    ├── 3_convert_gb_to_fasta.py
    ├── COI_all_Insecta.gb
    ├── COI_all_Insecta.fasta
    ├── convert_skipped_ids.txt
    └── taxidmap.txt
```

fasta 파일 전환이 완료되면 해당 파일과 taxonomy ID 정보를 결합한 레퍼런스 DB 제작을 진행한다.
```
mmseqs createdb COI_all_Insecta.fasta cox1.refDB

mmseqs createtaxdb cox1.refDB tmp —ncbi-tax-dump taxonomy/ —tax-mapping-file taxidmap.txt
```
작업이 완료된 작업 디렉토리 구성은 다음과 같다:
```
.
├── analysis
└── DB
    ├── 1_combine_gb_files.py
    ├── 2_taxidmap.py
    ├── 3_convert_gb_to_fasta.py
    ├── COI_all_Insecta.gb
    ├── COI_all_Insecta.fasta
    ├── convert_skipped_ids.txt
    ├── cox1.refDB
    ├── cox1.refDB.dbtype
    ├── cox1.refDB_h
    ├── cox1.refDB_h.dbtype
    ├── cox1.refDB_h.index
    ├── cox1.refDB.index
    ├── cox1.refDB.lookup
    ├── cox1.refDB.source
    └── taxidmap.txt
```

## Taxonomy Analysis
