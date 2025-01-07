#!/bin/bash
#SBATCH -J tax_db                    # 작업 이름
#SBATCH -p Node7                      # 실행 노드
#SBATCH -n 4                          # 사용 CPU 수
#SBATCH --time=UNLIMITED              # 실행 시간 제한 없음

# 2_create_taxDB.sh
# 설명:
# NCBI에서 Taxonomy 데이터를 다운로드한 뒤, MMseqs2에서 사용할 수 있는 Taxonomy DB를 생성하는 스크립트입니다.

# 주요 작업:
# 1. NCBI Taxonomy 데이터를 다운로드(taxdump.tar.gz) 및 압축 해제.
# 2. MMseqs2의 createtaxdb 명령어를 이용해 Taxonomy DB 생성.
# 3. 생성된 Taxonomy DB를 지정된 경로에 저장.

# 2_create_taxDB.sh
# Description: 
# This script downloads the NCBI Taxonomy database and creates a taxonomy database for MMseqs2.
# The output taxonomy database is used for annotating sequences with taxonomic information.

# Steps:
# 1. Download and extract NCBI Taxonomy data (taxdump.tar.gz).
# 2. Use the NCBI taxonomy data to create an MMseqs2-compatible taxonomy database.
# 3. Outputs the taxonomy database at the specified location.

# ================================
# 경로 설정
# ================================
BASE_DIR="/storage2/jihoonkim/eDNA"              # 기본 디렉토리
NCBI_DIR="$BASE_DIR/ncbi"                        # NCBI 데이터 디렉토리
NCBI_TAX_DUMP="$NCBI_DIR/taxonomy"               # NCBI Taxonomy 디렉토리
MMSEQS_DIR="$BASE_DIR/mmseqsDB"                  # MMseqs2 디렉토리
COX1_DB="$MMSEQS_DIR/cox1/cox1_refDB"            # MMseqs2 Cox1 DB
TAXONOMY_DB="$MMSEQS_DIR/taxonomy/mmseqs_taxDB"  # MMseqs2 Taxonomy DB

# ================================
# NCBI Taxonomy 데이터 다운로드 및 압축 해제
# ================================
mkdir -p "$NCBI_TAX_DUMP"
wget -q -O "$NCBI_TAX_DUMP/taxdump.tar.gz" https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && \
tar -xzf "$NCBI_TAX_DUMP/taxdump.tar.gz" -C "$NCBI_TAX_DUMP" && \
echo "NCBI Taxonomy data downloaded and extracted to $NCBI_TAX_DUMP"

# ================================
# MMseqs2 Taxonomy DB 생성
# ================================
mkdir -p "$(dirname "$TAXONOMY_DB")"
mmseqs createtaxdb "$COX1_DB" "$TAXONOMY_DB" --ncbi-tax-dump "$NCBI_TAX_DUMP"

# ================================
# 완료 메시지
# ================================
echo "MMseqs2 Taxonomy database created at $TAXONOMY_DB"

