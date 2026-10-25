#!/usr/bin/env bash
#SBATCH --job-name=6_swarm_MOTU
#SBATCH --output=6_swarm_%j.out
#SBATCH --error=6_swarm_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=Node4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# ============================================================
# PURPOSE
#   Cluster subsequence-resolved ASVs into MOTUs using Swarm
#   at a specified d-value while preserving sample-wise
#   abundance information.
#
# ANALYSIS
#   1. Prepare Swarm input FASTA with global ASV abundances.
#   2. Run Swarm clustering at the specified d-value.
#   3. Assign ASVs to MOTUs and rename MOTUs systematically.
#   4. Export representative sequences for each MOTU.
#   5. Aggregate ASV abundances into sample-by-MOTU tables.
#   6. Verify ASV membership, MOTU counts, and read preservation.
#
# OUTPUT
#   6_swarm_MOTU/ALL_dX/ALL_dX.swarm
#       : raw Swarm clustering result
#   6_swarm_MOTU/ALL_dX/ALL_dX_members.tsv
#       : ASV-to-MOTU membership table
#   6_swarm_MOTU/ALL_dX/ALL_dX_MOTU.fasta
#       : representative sequence for each MOTU
#   6_swarm_MOTU/ALL_dX/ALL_dX_MOTU_table.tsv
#       : sample-by-MOTU abundance table
#   6_swarm_MOTU/ALL_dX/ALL_dX_MOTU_abundance_long.tsv
#       : long-format MOTU abundance table
#   6_swarm_MOTU/ALL_dX/ALL_dX_summary.txt
#       : clustering and read-count summary
# ============================================================

set -euo pipefail
start_time=$(date +%s)

# ==== Environment and Paths ====
source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate dada2_env

WORK_DIR="/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
D_VAL="${1:-1}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

RESOLVED_FASTA="${WORK_DIR}/5_subsequence/ASVs.subsequence_resolved.fasta"
ABUND_RDS="${WORK_DIR}/5_subsequence/ASV_abundance_resolved.rds"

OUT_DIR="${WORK_DIR}/6_swarm_MOTU/ALL_d${D_VAL}"
POOL_NAME="ALL_d${D_VAL}"
mkdir -p "$OUT_DIR"

SWARM_IN_FA="${OUT_DIR}/${POOL_NAME}_swarm_input.fasta"
ASV_SAMPLE_TSV="${OUT_DIR}/${POOL_NAME}_ASV_sample_counts.tsv"
SAMPLE_ORDER="${OUT_DIR}/${POOL_NAME}_sample_order.txt"

SWARM_FILE="${OUT_DIR}/${POOL_NAME}.swarm"
STATS_FILE="${OUT_DIR}/${POOL_NAME}.stats"
SEEDS_RAW="${OUT_DIR}/${POOL_NAME}_seeds_raw.fasta"

RENAME_TSV="${OUT_DIR}/${POOL_NAME}_rename.tsv"
MEMBERS_TSV="${OUT_DIR}/${POOL_NAME}_members.tsv"
MOTU_FASTA="${OUT_DIR}/${POOL_NAME}_MOTU.fasta"

MOTU_LONG="${OUT_DIR}/${POOL_NAME}_MOTU_abundance_long.tsv"
MOTU_TABLE="${OUT_DIR}/${POOL_NAME}_MOTU_table.tsv"
SUMMARY="${OUT_DIR}/${POOL_NAME}_summary.txt"

[[ -s "$RESOLVED_FASTA" ]] || { echo "[ERROR] Missing $RESOLVED_FASTA"; exit 1; }
[[ -s "$ABUND_RDS" ]] || { echo "[ERROR] Missing $ABUND_RDS"; exit 1; }

command -v swarm >/dev/null 2>&1 || { echo "[ERROR] swarm not found"; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Rscript not found"; exit 1; }

echo "========================================================"
echo "[INFO] Swarm MOTU clustering"
echo "[INFO] d value : $D_VAL"
echo "[INFO] threads : $THREADS"
echo "[INFO] output  : $OUT_DIR"
echo "========================================================"

# ==== 1. Prepare Swarm Input ====
echo "[STEP 1] Preparing Swarm input..."

Rscript - "$RESOLVED_FASTA" "$ABUND_RDS" \
    "$SWARM_IN_FA" "$ASV_SAMPLE_TSV" "$SAMPLE_ORDER" <<'R_EOF'

args <- commandArgs(trailingOnly=TRUE)

fasta_file <- args[1]
abund_rds <- args[2]
swarm_fa <- args[3]
asv_tsv <- args[4]
sample_order <- args[5]

read_fasta <- function(path) {
  x <- readLines(path)
  h <- which(substr(x, 1, 1) == ">")
  ids <- sub("^>(\\S+).*", "\\1", x[h])
  ends <- c(h[-1] - 1, length(x))
  seqs <- mapply(function(a, b)
    paste0(x[(a + 1):b], collapse=""),
    h, ends, USE.NAMES=FALSE)
  setNames(toupper(seqs), ids)
}

seqs <- read_fasta(fasta_file)
abund <- as.matrix(readRDS(abund_rds))

if (anyDuplicated(names(seqs)))
    stop("Duplicate ASV IDs in resolved FASTA")

if (!setequal(names(seqs), colnames(abund)))
    stop("Resolved FASTA and abundance table contain different ASVs")

abund <- abund[, names(seqs), drop=FALSE]

if (sum(abund) <= 0)
    stop("Resolved abundance table contains no reads")

sizes <- colSums(abund)

con <- file(swarm_fa, "w")
for (id in names(seqs))
    writeLines(c(paste0(">", id, ";size=", sizes[id], ";"), seqs[id]), con)
close(con)

idx <- which(abund > 0, arr.ind=TRUE)
asv_long <- data.frame(
    ASV=colnames(abund)[idx[,2]],
    sample=rownames(abund)[idx[,1]],
    count=as.numeric(abund[idx]),
    stringsAsFactors=FALSE
)

write.table(
    asv_long, asv_tsv,
    sep="\t", quote=FALSE, row.names=FALSE
)

writeLines(rownames(abund), sample_order)

cat("[R INFO] Samples entering Swarm :", nrow(abund), "\n")
cat("[R INFO] ASVs entering Swarm    :", ncol(abund), "\n")
cat("[R INFO] Reads entering Swarm   :", sum(abund), "\n")
R_EOF

# ==== 2. Run Swarm Clustering ====
echo "[STEP 2] Running Swarm at d=${D_VAL}..."

if [[ "$D_VAL" -eq 1 ]]; then
    swarm \
        -d "$D_VAL" \
        -f \
        -t "$THREADS" \
        -z \
        -s "$STATS_FILE" \
        -o "$SWARM_FILE" \
        -w "$SEEDS_RAW" \
        "$SWARM_IN_FA"
else
    swarm \
        -d "$D_VAL" \
        -t "$THREADS" \
        -z \
        -s "$STATS_FILE" \
        -o "$SWARM_FILE" \
        -w "$SEEDS_RAW" \
        "$SWARM_IN_FA"
fi

[[ -s "$SWARM_FILE" ]] || { echo "[ERROR] Swarm output missing"; exit 1; }
[[ -s "$SEEDS_RAW" ]] || { echo "[ERROR] Swarm seed FASTA missing"; exit 1; }

# ==== 3. Build MOTU Membership and Representative FASTA ====
echo "[STEP 3] Building MOTU membership..."

awk -v d="$D_VAL" '
{
    seed=$1
    sub(/;size=[0-9]+;$/, "", seed)
    motu=sprintf("d%d_MOTU%06d", d, NR)
    print seed "\t" motu
}' "$SWARM_FILE" > "$RENAME_TSV"

awk -v d="$D_VAL" '
{
    motu=sprintf("d%d_MOTU%06d", d, NR)
    for (i=1; i<=NF; i++) {
        id=$i
        sub(/;size=[0-9]+;$/, "", id)
        print id "\t" motu
    }
}' "$SWARM_FILE" > "$MEMBERS_TSV"

awk '
BEGIN { FS="\t" }
NR==FNR {
    map[$1]=$2
    next
}
/^>/ {
    id=substr($0,2)
    sub(/;size=[0-9]+;.*/, "", id)

    if (!(id in map)) {
        print "[ERROR] Seed not found in rename table: " id > "/dev/stderr"
        exit 1
    }

    print ">" map[id]
    next
}
{ print }
' "$RENAME_TSV" "$SEEDS_RAW" > "$MOTU_FASTA"

# ==== 4. Aggregate ASV Abundance to MOTUs ====
echo "[STEP 4] Building MOTU abundance tables..."

Rscript - "$ASV_SAMPLE_TSV" "$MEMBERS_TSV" "$SAMPLE_ORDER" \
    "$MOTU_LONG" "$MOTU_TABLE" <<'R_EOF'

args <- commandArgs(trailingOnly=TRUE)

asv <- read.delim(args[1], stringsAsFactors=FALSE)
members <- read.delim(
    args[2],
    header=FALSE,
    col.names=c("ASV", "MOTU"),
    stringsAsFactors=FALSE
)
samples <- readLines(args[3])
long_out <- args[4]
table_out <- args[5]

if (anyDuplicated(members$ASV))
    stop("An ASV occurs in multiple Swarm MOTUs")

missing <- setdiff(asv$ASV, members$ASV)
if (length(missing))
    stop("ASVs missing from membership: ",
         paste(head(missing,10), collapse=", "))

x <- merge(asv, members, by="ASV", all.x=TRUE)

motu_long <- aggregate(
    count ~ MOTU + sample,
    data=x,
    FUN=sum
)
motu_long <- motu_long[
    order(motu_long$MOTU, motu_long$sample),
]

write.table(
    motu_long, long_out,
    sep="\t", quote=FALSE, row.names=FALSE
)

motus <- sort(unique(members$MOTU))
mat <- matrix(
    0,
    nrow=length(motus),
    ncol=length(samples),
    dimnames=list(motus, samples)
)

for (i in seq_len(nrow(motu_long)))
    mat[motu_long$MOTU[i], motu_long$sample[i]] <- motu_long$count[i]

out <- data.frame(
    "#OTU ID"=rownames(mat),
    mat,
    check.names=FALSE
)

write.table(
    out, table_out,
    sep="\t", quote=FALSE, row.names=FALSE
)
R_EOF

# ==== 5. Sanity Checks ====
echo "[STEP 5] Running sanity checks..."

INPUT_ASVS=$(grep -c '^>' "$SWARM_IN_FA" || true)
MEMBER_ASVS=$(wc -l < "$MEMBERS_TSV")
N_MOTUS=$(wc -l < "$SWARM_FILE")
N_SEEDS=$(grep -c '^>' "$MOTU_FASTA" || true)
N_SAMPLES=$(wc -l < "$SAMPLE_ORDER")

INPUT_READS=$(awk -F';size=' '/^>/ {
    split($2,a,";")
    sum+=a[1]
} END {printf "%.0f\n", sum}' "$SWARM_IN_FA")

MOTU_READS=$(awk 'BEGIN{FS="\t"} NR>1 {
    for(i=2;i<=NF;i++) sum+=$i
} END {printf "%.0f\n", sum}' "$MOTU_TABLE")

[[ "$INPUT_ASVS" -eq "$MEMBER_ASVS" ]] || {
    echo "[ERROR] Swarm input ASVs != membership ASVs"
    exit 1
}

[[ "$N_MOTUS" -eq "$N_SEEDS" ]] || {
    echo "[ERROR] MOTU count != representative sequence count"
    exit 1
}

[[ "$INPUT_READS" -eq "$MOTU_READS" ]] || {
    echo "[ERROR] Input reads != MOTU table reads"
    exit 1
}

{
    echo "============================================================"
    echo "Swarm MOTU summary"
    echo "d value                  : $D_VAL"
    echo "Samples                  : $N_SAMPLES"
    echo "ASVs entering Swarm      : $INPUT_ASVS"
    echo "MOTUs                    : $N_MOTUS"
    echo "Input reads              : $INPUT_READS"
    echo "MOTU table reads         : $MOTU_READS"
    echo "Membership ASVs          : $MEMBER_ASVS"
    echo "Representative sequences : $N_SEEDS"
    echo "============================================================"
} | tee "$SUMMARY"

elapsed=$(($(date +%s)-start_time))
printf "Total elapsed time: %dm %ds\n" "$((elapsed/60))" "$((elapsed%60))"
echo "Results directory: $OUT_DIR"
