#!/bin/sh
#SBATCH -J cox1_db
#SBATCH -p Node7
#SBATCH -n 4
#SBATCH --time=UNLIMITED

# taxonomy ID를 포함한 NCBI cox1 nt DB 생성 및 MMseqs2 전용 변환

# 디렉토리 설정
cd storage2/jihoonkim/eDNA/ncbi/cox1

# 출력 파일 경로
GENBANK_OUTPUT="cox1_sequences.gb"
FASTA_OUTPUT="cox1_sequences_with_taxid.fasta"

# 작업 시작 시간 기록
start_time=$(date +%s)

# 조건에 맞는 서열 검색 및 다운로드 (GenBank 형식)
echo "Starting cox1 gene sequence download in GenBank format..."
esearch -db nucleotide -query "cox1[Gene] NOT mitogenome[Title]" | \
efetch -format gb > $GENBANK_OUTPUT

# 다운로드 확인
if [ -s $GENBANK_OUTPUT ]; then
    echo "Download complete! GenBank results saved in $GENBANK_OUTPUT"
else
    echo "Download failed or no sequences found. Please check your query."
    exit 1
fi

# Python 스크립트를 사용하여 Taxonomy ID 포함 FASTA로 변환
echo "Converting GenBank to FASTA with Taxonomy ID..."
python << EOF
from Bio import SeqIO
from Bio.Entrez import efetch, read

def get_taxid(organism):
    """Fetch Taxonomy ID for an organism using NCBI Entrez API."""
    handle = efetch(db="taxonomy", term=organism, retmode="xml")
    records = read(handle)
    if records:
        return records[0]["TaxId"]
    return "unknown"

with open("$FASTA_OUTPUT", "w") as fasta_out:
    for record in SeqIO.parse("$GENBANK_OUTPUT", "genbank"):
        organism = record.annotations["organism"]
        taxid = get_taxid(organism)
        header = f">{record.id} {organism} [taxid={taxid}]"
        fasta_out.write(f"{header}\n{record.seq}\n")

print("Conversion complete! FASTA with Taxonomy ID saved in", "$FASTA_OUTPUT")
EOF

# mmseqs 전용 reference DB로 변환
mmseqs createdb /storage2/jihoonkim/eDNA/ncbi/cox1/cox1_sequences_with_taxid.fasta \
/storage2/jihoonkim/eDNA/ncbi/cox1/cox1_refDB

# 작업 종료 시간 기록
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

# 작업 시간 출력
echo "Total elapsed time: $((elapsed_time / 60)) minutes and $((elapsed_time % 60)) seconds"
