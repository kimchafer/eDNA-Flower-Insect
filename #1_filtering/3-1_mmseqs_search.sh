#!/bin/bash
#SBATCH --job-name=3-1_mmseqs
#SBATCH --output=3-1_mmseqs_%j.out
#SBATCH --error=3-1_mmseqs_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --partition=Node6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=400G

# ================================================================================================
# PURPOSE
#   Screen DADA2-derived ASVs against a nucleotide reference
#   database using MMseqs2 and retrieve taxonomic information for downstream biological filtering.
#
# ANALYSIS
#   1. Convert ASV FASTA sequences into an MMseqs2 query database.
#   2. Assign taxonomy using the top nucleotide hit.
#   3. Perform nucleotide similarity searches against the reference DB.
#   4. Export alignment statistics and taxonomic information for all hits.
#   5. Select the best-scoring hit for each ASV.
#
# OUTPUT
#   3_mmseqs2_filtering/results/ASVs_taxonomy_result.tsv
#       : MMseqs2 taxonomy assignments
#   3_mmseqs2_filtering/results/ASVs_hits.tsv
#       : all nucleotide-search hits with alignment and taxonomy information
#   3_mmseqs2_filtering/results/ASVs_best_hits.tsv
#       : highest-scoring reference hit for each ASV
#   3_mmseqs2_filtering/db/
#       : MMseqs2 query database
# ================================================================================================

set -euo pipefail
IFS=$'\n\t'

# ==== Environment and Paths ====
THREADS=${SLURM_CPUS_PER_TASK:-30}
export OMP_NUM_THREADS="$THREADS"

source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate mm

WORK_DIR="/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
INPUT_FASTA="${WORK_DIR}/2_DADA2/ASVs.nochim.fasta"
REF_DB="/storage/woojunbang/mmseq2_for_taxonomy/nt.fnaDB"

OUT_DIR="${WORK_DIR}/3_mmseqs2_filtering"
DB_DIR="${OUT_DIR}/db"
RESULTS_DIR="${OUT_DIR}/results"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$DB_DIR" "$RESULTS_DIR" "$TMP_DIR"

QUERY_DB="${DB_DIR}/ASVs_queryDB"
TAX_RESULT="${RESULTS_DIR}/ASVs_taxonomy_result"
TAX_TSV="${RESULTS_DIR}/ASVs_taxonomy_result.tsv"
HIT_RESULT="${RESULTS_DIR}/ASVs_search_result"
HIT_TSV="${RESULTS_DIR}/ASVs_hits.tsv"
BEST_HIT_TSV="${RESULTS_DIR}/ASVs_best_hits.tsv"

[[ -s "$INPUT_FASTA" ]] || { echo "[ERROR] Missing input: $INPUT_FASTA"; exit 1; }
[[ -e "${REF_DB}.dbtype" || -e "$REF_DB" ]] || {
    echo "[ERROR] Missing reference DB: $REF_DB"; exit 1;
}

echo "[INFO] Input : $INPUT_FASTA"
echo "[INFO] DB    : $REF_DB"
echo "[INFO] Output: $OUT_DIR"

# ==== 1. Create Query Database ====
echo "[STEP 1] Creating query DB..."
mmseqs createdb "$INPUT_FASTA" "$QUERY_DB"

# ==== 2. Taxonomic Assignment ====
echo "[STEP 2] Running taxonomy..."
mmseqs taxonomy "$QUERY_DB" "$REF_DB" "$TAX_RESULT" "$TMP_DIR/taxonomy" \
    --tax-lineage 1 \
    --threads "$THREADS" \
    --search-type 3 \
    --orf-filter 0 \
    --lca-mode 4 \
    --alignment-mode 3 \
    --strand 2 \
    -s 5 \
    -e 1e-5 \
    --min-seq-id 0.80 \
    -c 0.80 \
    --cov-mode 2 \
    --max-accept 20 \
    --max-rejected 50

mmseqs createtsv \
    "$QUERY_DB" "$TAX_RESULT" "$TAX_TSV" \
    --threads "$THREADS"

# ==== 3. Nucleotide Search ====
echo "[STEP 3] Running nucleotide search..."
mmseqs search "$QUERY_DB" "$REF_DB" "$HIT_RESULT" "$TMP_DIR/search" \
    --threads "$THREADS" \
    --search-type 3 \
    --alignment-mode 3 \
    --strand 2 \
    -s 5 \
    -e 1e-5 \
    --min-seq-id 0.80 \
    -c 0.80 \
    --cov-mode 2 \
    --max-accept 20 \
    --max-rejected 50

# ==== 4. Export All Hits ====
echo "[STEP 4] Exporting hit information..."
mmseqs convertalis \
    "$QUERY_DB" "$REF_DB" "$HIT_RESULT" "$HIT_TSV" \
    --format-output \
"query,target,pident,alnlen,qcov,tcov,evalue,bits,taxid,taxname,taxlineage,theader"

sed -i \
'1iquery\ttarget\tpident\talnlen\tqcov\ttcov\tevalue\tbits\ttaxid\ttaxname\ttaxlineage\ttarget_header' \
"$HIT_TSV"

# ==== 5. Extract Best Hit per ASV ====
echo "[STEP 5] Extracting best hits..."
python3 <<PY
import pandas as pd

hit_file = "$HIT_TSV"
best_file = "$BEST_HIT_TSV"

hits = pd.read_csv(hit_file, sep="\t", dtype=str)

for c in ["pident", "alnlen", "qcov", "tcov", "evalue", "bits"]:
    hits[c] = pd.to_numeric(hits[c], errors="coerce")

hits = hits.sort_values(
    ["query", "bits", "evalue"],
    ascending=[True, False, True]
)

best = hits.drop_duplicates("query", keep="first").copy()
best.to_csv(best_file, sep="\t", index=False)

print(f"[INFO] Total alignment hits : {len(hits)}")
print(f"[INFO] ASVs with best hit   : {len(best)}")
PY

N_INPUT=$(grep -c '^>' "$INPUT_FASTA")

echo "========================================================"
echo "            MMseqs2 ANALYSIS COMPLETE"
echo "========================================================"
echo "Input ASVs     : $N_INPUT"
echo "Taxonomy TSV   : $TAX_TSV"
echo "All hits TSV   : $HIT_TSV"
echo "Best hits TSV  : $BEST_HIT_TSV"
echo "========================================================"

rm -rf "$TMP_DIR"
