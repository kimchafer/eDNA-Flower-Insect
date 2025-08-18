# Flower-Visiting Insect Diversity in Korea Using eDNA Metabarcoding
eDNA analyses for insect pollinators

분석에 앞서 편의를 위해 작업 디렉토리 구조를 설정하자.

DB와 분석용 디렉토리를 우선 별개로 생성한다. ('.'으로 표현된 베이스 디렉토리 이름은 'eDNA'라고 가정하자)
```
.
├── analysis
└── DB
```
위 구조를 기본으로하여 아래의 분석 튜토리얼들이 진행된다.

__MMseqs2 작성자: Jihoon KIM 작성중__   
__Qimme 2 작성자:__ 

<details>
  <summary><h3>MMseqs2</h3></summary>
  <div>
    <ul>
  <details>
    <summary><h3>Reference DB Setup<h3></summary>


본 연구에 사용된 MMseqs2 전용 DB는 NCBI 기반의 곤충 COI 서열 모음이다.  
해당 데이터는 2025년 1월 16일 다운로드받아졌다.  
아래는 NCBI gb파일로부터 MMseqs2 전용 DB까지의 제작 과정 튜토리얼이며, 이미 만들어진 전용 DB는 리포지토리에 파일이 업로드되어 있으니 바로 사용하면 된다.  

셋업 당시 esearch 기능의 사용 불가로 인해 gb파일로부터 각 서열에 상응하는 accesstion number와 taxonomy ID를 추출해야만 했다.  
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
skipped_ids_file = "/eDNA/DB/convert_skipped_ids.txt"

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
Taxonomy 계층 정보를 통해 불필요한 데이터를 필터링하고, 해당 정보를 실제 mmseqs2 전용 DB에 입히기 위해서, taxdump = NCBI Taxonomy 원천 데이터 묶음이 필요하다.
```
mkdir taxdump
cd taxdump
wget -N https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz names.dmp nodes.dmp merged.dmp delnodes.dmp
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
    ├── taxidmap.txt
    └── taxdump
        ├── delnodes.dmp
        ├── merged.dmp
        ├── names.dmp
        ├── nodes.dmp
        └── taxdump.tar.gz
```
NCBI와 같은 공공데이터베이스에는 taxonomy info가 불완전한 데이터도 많이 존재한다.

이번 연구에선 비교적 구체적인 수준의 동정이 요구되므로, genus 수준까지의 데이터들만 레퍼런스로 인정하고, 그보다 정보력이 떨어지는, 즉 family 수준 이상의 정보 밖에 없는 서열들은 DB에서 제거하고자 한다.

위에서 만들어진 COI_all_Insecta.fasta, taxidmap.txt, taxdump를 베이스로 필터링을 수행한다: `4_DB_tax_filtering.sh`

이는 taxidmap.txt의 taxid를 taxdump의 계층 정보에 검색하여, 속 수준 이하의 taxonomy info를 가진 accession만 필터링하고, 해당 accession만으로 이뤄진 새로운 fasta 파일을 생성하는 과정이다.
```
#!/bin/bash
#SBATCH --job-name=DB_tax_filtering          # Job name
#SBATCH --time=24:00:00                      # Time limit
#SBATCH --mem=16G                            # Memory limit
#SBATCH -p Node7                             # Partition
#SBATCH -n 12                                # Number of CPUs
#SBATCH --output=DB_tax_filtering_%j.out     # Standard log
#SBATCH --error=DB_tax_filtering_%j.err      # Error log

###############################################################################
# 안전 설정
# - 에러 시 즉시 중단(-e), 설정 안 된 변수 사용 금지(-u), 파이프 에러 전파(-o pipefail)
# - 에러가 발생한 라인 번호를 출력하도록 trap 설정
###############################################################################
set -euo pipefail
trap 'echo "[ERROR] at line $LINENO"; exit 1' ERR

###############################################################################
# 경로/파일 설정
# - DB_DIR: 작업 디렉토리 (FASTA와 taxidmap.txt가 있는 곳)
# - FASTA: 원본 전체 COI FASTA
# - MAP  : Accession \t TaxID 매핑 파일 (mmseqs 전용 DB 만들기 전 단계에서 추출한 것)
# - TAXDUMP: NCBI taxdump 디렉토리 (nodes.dmp, names.dmp 필수 / merged.dmp, delnodes.dmp 선택)
# - OUT_BASENAME: 필터링 후 생성할 FASTA의 기본 파일명(확장자 .fasta 자동 부여)
###############################################################################
DB_DIR="/storage2/flower_eDNA/DB"
FASTA="${DB_DIR}/COI_all_Insecta.fasta"
MAP="${DB_DIR}/taxidmap.txt"
TAXDUMP="/storage2/flower_eDNA/DB/taxdump"
OUT_BASENAME="COI_good_tax"   # 결과: ${DB_DIR}/${OUT_BASENAME}.fasta

###############################################################################
# 산출물(중간/최종)
# - good_tax.ids.txt: 남길(KEEP) Accession 목록 (Genus/Species/종하위)
# - bad_tax.ids.txt : 제거할(DROP) Accession 목록 (Family 이상 / no rank / invalid 등)
# - tax_filter.report.tsv: 각 Accession에 대한 판정 로그(Accession, TaxID, rank, name, decision)
# - good_tax.regex.txt: 버전 유무 허용 정규식 목록 (^ACC(\.\d+)?$)
# - COI_good_tax.fasta: 최종 KEEP 서열만 추출된 FASTA
###############################################################################
GOOD_IDS="${DB_DIR}/good_tax.ids.txt"
BAD_IDS="${DB_DIR}/bad_tax.ids.txt"
REPORT="${DB_DIR}/tax_filter.report.tsv"
REGEX_FILE="${DB_DIR}/good_tax.regex.txt"
OUT_FASTA="${DB_DIR}/${OUT_BASENAME}.fasta"

###############################################################################
# 스레드 설정
# - seqkit 등 멀티스레드 도구에서 사용할 스레드 수를 SLURM 환경 변수에서 우선 획득
# - 없으면 기본 12로 설정
###############################################################################
THREADS="${SLURM_CPUS_ON_NODE:-${SLURM_CPUS_PER_TASK:-12}}"

###############################################################################
# 사전 점검(필수)
# - seqkit 설치 여부
# - FASTA / MAP / taxdump 기본 파일 존재 확인
# - 참고: taxdump는 최소 nodes.dmp, names.dmp 필요. 있으면 merged.dmp, delnodes.dmp도 활용
###############################################################################
echo "== 설정 확인 =="
echo "DB_DIR   : $DB_DIR"
echo "FASTA    : $FASTA"
echo "MAP      : $MAP"
echo "TAXDUMP  : $TAXDUMP"
echo "THREADS  : $THREADS"
echo "OUT_FASTA: $OUT_FASTA"
echo

command -v seqkit >/dev/null 2>&1 || { echo "[에러] seqkit이 설치되어 있어야 합니다."; exit 1; }
[[ -s "$FASTA" ]] || { echo "[에러] FASTA 파일이 보이지 않습니다: $FASTA"; exit 1; }
[[ -s "$MAP"   ]] || { echo "[에러] taxidmap.txt 파일이 보이지 않습니다: $MAP"; exit 1; }
[[ -s "$TAXDUMP/nodes.dmp" && -s "$TAXDUMP/names.dmp" ]] || {
  echo "[에러] $TAXDUMP에 nodes.dmp 또는 names.dmp가 없습니다."; exit 1; }

###############################################################################
# 1) TaxID → rank 판정 및 KEEP/DROP 분류 (Python)
#    - nodes.dmp: TaxID의 parent_tax_id, rank 정보를 제공
#    - names.dmp: TaxID의 scientific name(보고용)
#    - merged.dmp(선택): 옛 TaxID → 현재 유효 TaxID로 매핑(존재 시 적용)
#    - delnodes.dmp(선택): 삭제된 TaxID 목록(존재 시 DROP 처리)
#    - 규칙:
#       * KEEP: rank ∈ {species, genus}
#       * KEEP: subspecies/varietas/subgenus 등 '종하위'라도 species 조상이 존재하면 KEEP
#       * DROP: family 이상, no rank, 삭제된 TaxID, invalid TaxID
#    - 메모리 사용을 줄이기 위해 라인 단위 처리 + 이후 sort -u로 ID 중복 제거
###############################################################################
python3 - "$MAP" "$TAXDUMP" "$GOOD_IDS" "$BAD_IDS" "$REPORT" <<'PY'
import sys, os, csv

MAP, TAXDUMP, OUT_KEEP, OUT_DROP, OUT_REPORT = sys.argv[1:6]
NODES = os.path.join(TAXDUMP, "nodes.dmp")
NAMES = os.path.join(TAXDUMP, "names.dmp")
MERGED = os.path.join(TAXDUMP, "merged.dmp")
DELNODES = os.path.join(TAXDUMP, "delnodes.dmp")

ALLOW = {"genus","species"}
SUBLIKE = {"subspecies","forma","varietas","subvariety","infraspecies",
           "biotype","morph","pathovar","serovar","subgenus"}

# 1) nodes.dmp 로드: parent, rank
parent, rank = {}, {}
with open(NODES, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        p=[x.strip() for x in line.split("|")]
        if len(p)>=3:
            parent[p[0]] = p[1]
            rank[p[0]]   = p[2]

# 2) names.dmp 로드: scientific name (리포트용)
name = {}
with open(NAMES, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        p=[x.strip() for x in line.split("|")]
        if len(p)>=4 and p[3]=="scientific name":
            name[p[0]] = p[1]

# 3) merged.dmp 로드(선택): 옛 TaxID → 신규 TaxID 매핑
merged = {}
if os.path.exists(MERGED):
    with open(MERGED, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p=[x.strip() for x in line.split("|")]
            if len(p)>=2:
                old_id, new_id = p[0], p[1]
                merged[old_id] = new_id

# 4) delnodes.dmp 로드(선택): 삭제된 TaxID
deleted = set()
if os.path.exists(DELNODES):
    with open(DELNODES, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p=[x.strip() for x in line.split("|")]
            if p and p[0].isdigit():
                deleted.add(p[0])

def canonical_tid(tid: str) -> str:
    """merged/delnodes를 반영하여 최종 유효 TaxID 반환(없으면 원본 유지)"""
    if tid in merged:
        return merged[tid]
    return tid

def has_species_ancestor(tid: str) -> bool:
    """주어진 TaxID에서 위로 타고 올라가며 species 조상이 있는지 확인"""
    seen=set(); cur=tid
    while cur in parent and cur not in seen:
        seen.add(cur)
        if rank.get(cur,"") == "species":
            return True
        np = parent.get(cur)
        if not np or np == cur:
            break
        cur = np
    return False

def decide_keep(tid: str):
    """KEEP 여부와 실제 랭크 반환"""
    r = rank.get(tid, "")
    if r in ALLOW:
        return True, r
    if r in SUBLIKE:
        return has_species_ancestor(tid), r
    return False, r

# 라인 단위로 판정 및 파일에 스트리밍 기록
with open(OUT_KEEP,"w",encoding="utf-8") as f_keep, \
     open(OUT_DROP,"w",encoding="utf-8") as f_drop, \
     open(OUT_REPORT,"w",encoding="utf-8") as f_rep:

    w = csv.writer(f_rep, delimiter="\t", lineterminator="\n")
    w.writerow(["accession","taxid","rank","name","decision"])

    keep_cnt = drop_cnt = 0
    with open(MAP,"r",encoding="utf-8") as f:
        rdr = csv.reader(f, delimiter="\t")
        for rec in rdr:
            if not rec:
                continue
            acc = rec[0].strip()
            tid = rec[1].strip() if len(rec)>1 else ""

            # 숫자 여부 확인
            if not tid.isdigit():
                w.writerow([acc, tid, "NA", "", "drop_invalid_taxid"])
                f_drop.write(acc + "\n"); drop_cnt += 1
                continue

            # 삭제된 TaxID는 DROP
            if tid in deleted:
                w.writerow([acc, tid, "deleted", "", "drop_deleted_taxid"])
                f_drop.write(acc + "\n"); drop_cnt += 1
                continue

            # merged 매핑(옛→신규) 반영
            tid2 = canonical_tid(tid)

            # 랭크 판정
            ok, r = decide_keep(tid2)
            nm = name.get(tid2, "")
            w.writerow([acc, tid2, r if r else "no_rank", nm, "keep" if ok else "drop"])
            if ok:
                f_keep.write(acc + "\n"); keep_cnt += 1
            else:
                f_drop.write(acc + "\n"); drop_cnt += 1

print("[INFO] streaming classification done.")
PY

###############################################################################
# 1-보강) ID 중복 제거
# - 동일 Accession이 여러 번 등장하는 경우를 대비하여 고유화
# - LC_ALL=C: sort 성능 최적화
###############################################################################
export LC_ALL=C
sort -u -o "$GOOD_IDS" "$GOOD_IDS" || true
sort -u -o "$BAD_IDS"  "$BAD_IDS"  || true

echo "[INFO] KEEP IDs: $(wc -l < "$GOOD_IDS")"
echo "[INFO] DROP IDs: $(wc -l < "$BAD_IDS")"

###############################################################################
# 2) Accession.버전 혼재 대응 정규식 생성
# - FASTA 헤더의 첫 토큰이 Accession 또는 Accession.버전(MN123456 혹은 MN123456.1)
# - 정규식 형태: ^ACC(\.\d+)?$
#   -> 버전 유무 상관없이 매칭되도록 패턴 생성
# - core 기준(버전 제거)으로 중복 없이 생성
###############################################################################
awk '{
  split($1,a,"."); core=a[1];
  if(!(core in seen)){ print "^" core "(\\.\\d+)?$"; seen[core]=1 }
}' "$GOOD_IDS" > "$REGEX_FILE"

echo "[INFO] regex patterns: $(wc -l < "$REGEX_FILE")"

###############################################################################
# 3) FASTA에서 KEEP만 추출 (seqkit grep)
# - -n : FASTA의 "ID(name)"만 검색 (헤더 전체가 아니라 첫 토큰 기반)
# - -r : 정규식 사용
# - -f : 패턴 파일 입력
# - -j : 멀티스레드
# - 주의: 헤더 첫 토큰이 Accession(.버전)이어야 정확히 매칭됨 (일반적인 NCBI/BOLD FASTA는 보통 OK)
###############################################################################
echo "[INFO] original FASTA sequences: $(grep -c '^>' "$FASTA")"
seqkit grep -n -r -f "$REGEX_FILE" -j "$THREADS" "$FASTA" > "$OUT_FASTA"
echo "[INFO] filtered FASTA sequences: $(grep -c '^>' "$OUT_FASTA")"

###############################################################################
# 완료 안내 + 후속 팁
# - mmseqs DB 생성: mmseqs createdb ${OUT_FASTA} ${DB_DIR}/COI_good_tax.DB
# - 더 엄격하게 하고 싶다면:
#   * 종하위(SUBLIKE)도 제외하려면 Python의 SUBLIKE 처리 부분을 DROP으로 바꾸면 됨
# - 메모리 최적화:
#   * 16G → 8G로 줄여도 동작 가능할 수 있음(환경/파일 구조에 따라 다름)
###############################################################################
echo
echo "[DONE] 생성 파일"
echo " - KEEP ID 목록      : $GOOD_IDS"
echo " - DROP ID 목록      : $BAD_IDS"
echo " - 판정 리포트       : $REPORT"
echo " - 정규식 패턴       : $REGEX_FILE"
echo " - 최종 FASTA (KEEP) : $OUT_FASTA"


```
작업이 정상적으로 완료되었다면 작업 디렉토리 내 파일 구성은 다음과 같다:
```
.
├── analysis
└── DB
    ├── 1_combine_gb_files.py
    ├── 2_taxidmap.py
    ├── 3_convert_gb_to_fasta.py
    ├── 4_DB_tax_filtering.sh
    ├── COI_all_Insecta.gb
    ├── COI_all_Insecta.fasta
    ├── convert_skipped_ids.txt
    ├── taxidmap.txt
    └── taxdump
```
fasta 파일 전환이 완료되면 해당 파일과 taxonomy ID 정보를 결합한 레퍼런스 DB 제작을 진행한다.
```
mmseqs createdb COI_good_tax.fasta cox1.refDB

mmseqs createtaxdb cox1.refDB tmp --ncbi-tax-dump taxdump/ --tax-mapping-file good_tax.ids.txt
```
작업이 완료된 작업 디렉토리 구성은 다음과 같다:
```
.
├── analysis
└── DB
    ├── 1_combine_gb_files.py
    ├── 2_taxidmap.py
    ├── 3_convert_gb_to_fasta.py
    ├── 4_DB_tax_filtering.sh
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
    ├── taxidmap.txt
    ├── taxdump
    └── tmp
```
  </details>

  <details>
    <summary><h3>Taxonomy Analysis</h3></summary>

데이터베이스 구축이 완료되었으면 우리의 샘플 데이터를 동정하고 분류학적 결과를 만들어보자.

우선 편의를 위해 분석 디렉토리에 아래와 같이 하위 디렉토리를 생성한다:

```
.
├── analysis
│   ├── processed_data
│   ├── raw_data
│   ├── results
│   └── tmp
└── DB

```
raw_data 디렉토리에는 예시로 본 연구에서 사용된 제주도 벚나무의 eDNA 샘플을 사용하였다.
제주도 벚나무 총 6개 반복 가운데 2개 반복만이 라이브러리 제작에 이용되었다.  
* PS32  
* PS36
<br><br>  
```
.
├── analysis
│   ├── processed_data
│   ├── raw_data
│   │   ├── S2_PS32_Lib_S66_L001_R1_001.fastq.gz
│   │   ├── S2_PS32_Lib_S66_L001_R2_001.fastq.gz
│   │   ├── S2_PS36_Lib_S66_L001_R1_001.fastq.gz
│   │   └── S2_PS36_Lib_S66_L001_R2_001.fastq.gz
│   ├── results
│   └── tmp
└── DB

```
Paired-end 데이터이기 때문에 각 샘플별로 두 개의 fastq 파일이 준비되었다.  
다음 과정은 QC, 어댑터 서열 제거 및 저퀄리티 리드 필터링, Paired-End 데이터 병합, .fastq -> .fasta 변환 과정을 포함하고 있다: `A1_filtering.sh`

```
#!/bin/bash
#SBATCH --job-name=filter_merge_fasta        # Job name
#SBATCH --output=filter_merge_fasta_%j.out   # Standard output
#SBATCH --error=filter_merge_fasta_%j.err    # Error output
#SBATCH --time=UNLIMITED                     # No time limit
#SBATCH -p Node7                             # Partition
#SBATCH -n 12                                # Number of CPUs

# 샘플 이름 설정
sample=$1
if [ -z "$sample" ]; then
  echo "Error: No sample name provided. Run as: sbatch script_name.sh <sample>"
  exit 1
fi

# 작업 시작 시간 기록
start_time=$(date +%s)

# 디렉토리 설정
work_dir=/eDNA/analysis
raw_data_dir=${work_dir}/raw_data
processed_data_dir=${work_dir}/processed_data

# 1. Raw paired-end 데이터 필터링 (fastp 사용)
# 아래 코드에서 이름 형식에 관한 부분은 시퀀싱 이후 전달받은 샘플명에 따라 유동적으로 수정
echo "==> Running fastp for adapter trimming and quality filtering: $sample"
R1=$(find $raw_data_dir -name "S2_${sample}_Lib_*_R1_001.fastq.gz")
R2=$(find $raw_data_dir -name "S2_${sample}_Lib_*_R2_001.fastq.gz")

filtered_R1=${processed_data_dir}/${sample}_filtered_R1.fastq.gz
filtered_R2=${processed_data_dir}/${sample}_filtered_R2.fastq.gz
html_report=${processed_data_dir}/${sample}_fastp.html
json_report=${processed_data_dir}/${sample}_fastp.json

fastp -i $R1 -I $R2 \
      -o $filtered_R1 -O $filtered_R2 \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 100 \
      -w 16 -h $html_report -j $json_report

if [ $? -ne 0 ]; then
  echo "Error: fastp failed for sample $sample"
  exit 1
fi

# 2. Paired-End 데이터 병합 (PEAR 사용)
echo "==> Merging filtered paired-end reads for sample: $sample"
merged_prefix=${processed_data_dir}/${sample}_merged

pear -f $filtered_R1 -r $filtered_R2 -o $merged_prefix -j 8

merged_file=${merged_prefix}.assembled.fastq
if [ ! -f "$merged_file" ]; then
  echo "Error: Merged file not found for sample $sample."
  exit 1
fi

# 3. FASTQ -> FASTA 변환
echo "==> Converting FASTQ to FASTA for sample: $sample"
fasta_file=${processed_data_dir}/${sample}.fasta
seqtk seq -A $merged_file > $fasta_file
```
작업이 정상적으로 완료되었다면 작업 디렉토리 내 파일 구성은 다음과 같다:
```
.
├── analysis
│   ├── A1_filtering.sh
│   ├── processed_data
│   │   ├── PS32_fastp.json
│   │   ├── PS36_fastp.json
│   │   ├── PS32_filtered_R1.fastq.gz
│   │   ├── PS32_filtered_R2.fastq.gz
│   │   ├── PS36_filtered_R1.fastq.gz
│   │   ├── PS36_filtered_R2.fastq.gz
│   │   ├── PS32_merged.assembled.fastq
│   │   ├── PS36_merged.assembled.fastq
│   │   ├── PS32_merged.discarded.fastq
│   │   ├── PS36_merged.discarded.fastq
│   │   ├── PS32_merged.unassembled.forward.fastq
│   │   ├── PS36_merged.unassembled.forward.fastq
│   │   ├── PS32_merged.unassembled.reverse.fastq
│   │   ├── PS36_merged.unassembled.reverse.fastq
│   │   ├── PS32.fasta
│   │   └── PS36.fasta
│   ├── raw_data
│   ├── results
│   └── tmp
└── DB
```
두 샘플의 최종 fasta 파일이 PS32.fasta와 PS36.fasta로 생성되었다.  
둘 다 같은 제주도 벚나무 시료이므로 두 파일을 합쳐준다:
```
cat /eDNA/analysis/processed_data/PS32.fasta \
    /eDNA/analysis/processed_data/PS36.fasta \
    > /eDNA/analysis/processed_data/PS_JJ.fasta
```
작업이 완료되면 작업 디렉토리 내 PS_JJ.fasta가 새로 생성되었다:
```
.
├── analysis
│   ├── A1_filtering.sh
│   ├── processed_data
│   │   ├── PS32_fastp.json
│   │   ...
│   │   └── PS_JJ.fasta
│   ├── raw_data
│   ├── results
│   └── tmp
└── DB
```
Assembly없이 raw read를 그대로 동정에 사용하는 메타바코딩 특성 상 리드들의 길이와 그에 따른 정보력이 결과에 큰 영향을 미친다.  
불필요한 짧은 서열을 제거하면서 데이터의 용량을 줄이기 위해 리드들 사이에서 포함관계가 성립하는 경우 짧은 리드를 제거하는 과정을 거친다.  
길이와 서열 모두 100% 일치하거나 길이가 그보다 짧은 경우에만 제거되며, 해당 코드는 `A2_dereplication.sh`에서 찾아볼 수 있다:  



위 과정을 통해 PS_JJ_derep.fasta가 생성된 것을 확인할 수 있다:
```
.
├── analysis
│   ├── A1_filtering.sh
│   ├── A2_dereplication.sh
│   ├── processed_data
│   │   ├── PS32_fastp.json
│   │   ...
│   │   ├── PS_JJ_derep.fasta
│   │   └── PS_JJ.fasta
│   ├── raw_data
│   ├── results
│   └── tmp
└── DB
```
실제 PS_JJ.fasta에서 PS_JJ_derep.fasta으로 용량 변화가 108.9 MB -> 73.0 MB로 이뤄진 것을 확인할 수 있다.
<br><br>
<img width="3000" height="1800" alt="read_length_distribution" src="https://github.com/user-attachments/assets/e56d3872-a44e-4cce-af47-257579d47798" />
<br><br>
두 파일의 리드 길이 분포 그래프에서도 PS_JJ_derep.fasta의 불필요한 리드들이 제거된 것을 확인할 수 있다.

마지막 생성 파일을 인풋으로 하여 서열 동정 및 분류 분석을 수행해보자.



```
#!/bin/bash
#SBATCH --job-name=tax_analysis          # Job name
#SBATCH --output=full_analysis_%j.out    # Standard output
#SBATCH --error=full_analysis_%j.err     # Error output
#SBATCH --time=48:00:00                  # Time limit
#SBATCH -p Node7                         # Partition
#SBATCH -n 32                            # Number of CPUs

# 샘플 이름 설정
sample=$1
if [ -z "$sample" ]; then
  echo "Error: No sample name provided. Run as: sbatch script_name.sh <sample>"
  exit 1
fi

# 작업 시작 시간 기록
start_time=$(date +%s)

# 디렉토리 설정
work_dir=/eDNA/analysis
processed_data_dir=${work_dir}/processed_data
results_dir=${work_dir}/results
tmp_dir=${work_dir}/tmp
ref_db=/eDNA/DB/cox1.refDB  # Reference DB

# 1. MMseqs2 쿼리 데이터베이스 생성
echo "==> Creating MMseqs2 database for sample: $sample"
query_db=${processed_data_dir}/${sample}_queryDB
mmseqs createdb ${sample}_derep.fasta $query_db

# 2. Taxonomy 결과 경로 설정
taxonomy_result=${results_dir}/${sample}_taxonomy_result
tsv_result=${results_dir}/${sample}_taxonomy_result.tsv
krona_report=${results_dir}/${sample}_krona_report.html
taxonomy_tmp_dir=${tmp_dir}/${sample}_taxonomy_tmp
mkdir -p $taxonomy_tmp_dir

# 기존 Taxonomy 결과 파일 삭제
if [ -f "${taxonomy_result}.dbtype" ]; then
  echo "==> Previous taxonomy result exists. Deleting..."
  rm -rf "${taxonomy_result}"*
fi

# 3. MMseqs2 Taxonomy 분석 수행
echo "==> Running taxonomy analysis for sample: $sample"
mmseqs taxonomy $query_db $ref_db $taxonomy_result $taxonomy_tmp_dir \
  --tax-lineage 1 --threads 32 --search-type 3 --orf-filter 0 --min-seq-id 0.98 \
  --alignment-mode 4 --min-aln-len 300
if [ $? -ne 0 ]; then
  echo "Error: Taxonomy analysis failed for sample $sample."
  exit 1
fi

# 4. Taxonomy 결과를 TSV로 변환
echo "==> Converting taxonomy results to TSV for sample: $sample"
mmseqs createtsv $query_db $taxonomy_result $tsv_result --threads 16
if [ $? -ne 0 ]; then
  echo "Error: Failed to convert taxonomy results to TSV for sample $sample."
  exit 1
fi

# 5. Krona HTML 보고서 생성
echo "==> Generating Krona HTML report for sample: $sample"
mmseqs taxonomyreport $ref_db $taxonomy_result $krona_report --threads 16 --report-mode 1
if [ $? -ne 0 ]; then
  echo "Error: Failed to generate Krona report for sample $sample."
  exit 1
fi

# 작업 종료 시간 기록 및 소요 시간 계산
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "==> Taxonomy analysis completed for sample: $sample"
echo "Results summary:"
echo "  - Taxonomy DB: $taxonomy_result"
echo "  - TSV Result: $tsv_result"
echo "  - Krona HTML Report: $krona_report"
echo "Total elapsed time: $((elapsed_time / 60)) minutes and $((elapsed_time % 60)) seconds"
```





```
.
├── analysis
│   ├── A1_filtering.sh
│   ├── A2_dereplication.sh
│   ├── A3_tax_analysis.sh
│   ├── processed_data
│   ├── raw_data
│   ├── results
│   └── tmp
└── DB
```




  </details>

  </div>
</details>








<details>
  <summary><h3>QIMME 2</h3></summary>
  <!-- 내용 -->
</details>
