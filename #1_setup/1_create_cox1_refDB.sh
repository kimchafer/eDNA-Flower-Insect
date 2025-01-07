#!/bin/sh
#SBATCH -J cox1_db
#SBATCH -p Node7
#SBATCH -n 4
#SBATCH --time=UNLIMITED

# taxonomy ID를 포함한 NCBI cox1 nt DB 생성 및 MMseqs2 전용 변환

# ================================
# 기본 경로 설정
# ================================
BASE_DIR="/storage2/jihoonkim/eDNA"              # 기본 작업 디렉토리
NCBI_DIR="$BASE_DIR/ncbi"                        # NCBI 관련 데이터 디렉토리
COX1_DIR="$NCBI_DIR/cox1"                        # Cox1 DB 작업 디렉토리
MMSEQS_DIR="$BASE_DIR/mmseqsDB"                  # MMseqs2 전용 데이터 디렉토리
MMSEQS_COX1_DIR="$MMSEQS_DIR/cox1"               # MMseqs2 Cox1 DB 디렉토리

# ================================
# 파일 경로 변수 지정
# ================================
GENBANK_FILE="$COX1_DIR/cox1_sequences.gb"       # GenBank 형식 파일
FASTA_FILE="$COX1_DIR/cox1_sequences_with_taxid.fasta" # FASTA 파일
MMSEQS_DB="$MMSEQS_COX1_DIR/cox1_refDB"         # MMseqs2 전용 DB

# ================================
# 작업 시작 시간 기록
# ================================
start_time=$(date +%s)

# ================================
# GenBank 파일 다운로드
# ================================
echo "Starting cox1 gene sequence download in GenBank format..."
mkdir -p "$COX1_DIR"  # 작업 디렉토리 생성
esearch -db nucleotide -query "cox1[Gene] NOT mitogenome[Title]" | \
efetch -format gb > "$GENBANK_FILE"

# 다운로드 확인
if [ -s "$GENBANK_FILE" ]; then
    echo "Download complete! GenBank results saved in $GENBANK_FILE"
else
    echo "Download failed or no sequences found. Please check your query."
    exit 1
fi

# ================================
# GenBank -> FASTA 변환
# ================================
echo "==> Converting GenBank to FASTA with Taxonomy ID"
python3 <<EOF
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

# 입력 및 출력 파일 경로 설정
genbank_file = "${GENBANK_FILE}"
output_fasta = "${FASTA_FILE}"

# GenBank에서 FASTA로 변환
with open(genbank_file, "r") as gb_file, open(output_fasta, "w") as fasta_out:
    for record in SeqIO.parse(gb_file, "genbank"):
        try:
            # 서열 데이터가 정의되지 않은 경우 예외 처리
            sequence = str(record.seq)  # 서열 데이터 추출
            
            # Extract metadata from GenBank
            seq_id = record.id  # 서열 ID
            organism = record.annotations.get("organism", "unknown")  # 생물 학명
            
            # Taxonomy ID 추출
            tax_id = "unknown"  # 기본값
            for feature in record.features:
                if feature.type == "source":  # "source" Feature에서 Taxonomy ID 찾기
                    taxon_info = feature.qualifiers.get("db_xref", [])
                    for entry in taxon_info:
                        if entry.startswith("taxon:"):
                            tax_id = entry.split(":")[1]  # TaxID 추출
                            break

            # Write to FASTA
            header = f">{seq_id} {organism} [taxid={tax_id}]"
            fasta_out.write(f"{header}\n{sequence}\n")

        except (UndefinedSequenceError, ValueError) as e:
            print(f"Warning: Skipping record {record.id} due to error: {e}")

print(f"FASTA with Taxonomy ID saved to: {output_fasta}")
EOF

# 변환 성공 여부 확인
if [ $? -ne 0 ]; then
    echo "Error: GenBank to FASTA conversion failed."
    exit 1
fi

# ================================
# MMseqs2 전용 DB로 변환
# ================================
echo "==> Creating MMseqs2 database"
mkdir -p "$MMSEQS_COX1_DIR"  # MMseqs2 Cox1 디렉토리 생성
mmseqs createdb "$FASTA_FILE" "$MMSEQS_DB"

# 변환 성공 여부 확인
if [ $? -ne 0 ]; then
    echo "Error: MMseqs2 database creation failed."
    exit 1
fi

# ================================
# 작업 종료 시간 기록
# ================================
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

# ================================
# 작업 시간 출력
# ================================
echo "Total elapsed time: $((elapsed_time / 60)) minutes and $((elapsed_time % 60)) seconds"

