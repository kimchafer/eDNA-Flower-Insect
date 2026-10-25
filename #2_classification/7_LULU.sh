#!/usr/bin/env bash
#SBATCH --job-name=7_LULU
#SBATCH --output=7_LULU_%j.out
#SBATCH --error=7_LULU_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=Node6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# ================================================================================================
# PURPOSE
#   Curate Swarm d=3 MOTUs using LULU to merge likely erroneous or redundant daughter MOTUs
#   based on sequence similarity, abundance, and sample co-occurrence.
#
# ANALYSIS
#   1. Perform an all-against-all VSEARCH comparison of Swarm d=3 representative MOTU sequences.
#   2. Remove samples containing zero total reads before LULU.
#   3. Run LULU using a minimum sequence match of 97%.
#   4. Use abundance and co-occurrence information to identify and merge likely daughter MOTUs.
#   5. Export the curated abundance table, MOTU map, representative FASTA, and LULU result object.
#   6. Verify MOTU counts and preservation of total reads.
#
# OUTPUT
#   7_LULU/ALL_d3_vsearch_matchlist.tsv
#       : pairwise MOTU matches with >=97% sequence identity
#   7_LULU/LULU_curated_MOTU_table.tsv
#       : final sample-by-MOTU abundance table after LULU
#   7_LULU/LULU_MOTU_map.tsv
#       : LULU parent-child curation map
#   7_LULU/LULU_curated_MOTU.fasta
#       : representative sequences of retained MOTUs
#   7_LULU/LULU_result.rds
#       : complete LULU result object
#   7_LULU/LULU_summary.tsv
#       : final LULU curation summary
#   7_LULU/all_zero_samples_removed.txt
#       : samples removed because of zero total reads
# ================================================================================================

set -euo pipefail
start_time=$(date +%s)

# ==== Environment and Paths ====
source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate dada2_env

WORK_DIR="/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
INPUT_DIR="${WORK_DIR}/6_swarm_MOTU/ALL_d3"
OUT_DIR="${WORK_DIR}/7_LULU"

FASTA="${INPUT_DIR}/ALL_d3_MOTU.fasta"
TABLE="${INPUT_DIR}/ALL_d3_MOTU_table.tsv"

MATCHLIST="${OUT_DIR}/ALL_d3_vsearch_matchlist.tsv"
ZERO_SAMPLE_LOG="${OUT_DIR}/all_zero_samples_removed.txt"

CURATED_TABLE="${OUT_DIR}/LULU_curated_MOTU_table.tsv"
MOTU_MAP="${OUT_DIR}/LULU_MOTU_map.tsv"
CURATED_FASTA="${OUT_DIR}/LULU_curated_MOTU.fasta"
LULU_RDS="${OUT_DIR}/LULU_result.rds"
SUMMARY="${OUT_DIR}/LULU_summary.tsv"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MINIMUM_MATCH=97

mkdir -p "$OUT_DIR"

[[ -s "$FASTA" ]] || {
    echo "[ERROR] Missing $FASTA"
    exit 1
}

[[ -s "$TABLE" ]] || {
    echo "[ERROR] Missing $TABLE"
    exit 1
}

command -v vsearch >/dev/null 2>&1 || {
    echo "[ERROR] vsearch not found"
    exit 1
}

command -v Rscript >/dev/null 2>&1 || {
    echo "[ERROR] Rscript not found"
    exit 1
}

Rscript -e '
if (!requireNamespace("lulu", quietly=TRUE))
    stop("R package lulu is not installed")
if (!requireNamespace("dplyr", quietly=TRUE))
    stop("R package dplyr is not installed")
'

echo "============================================================"
echo "[INFO] LULU curation"
echo "[INFO] Swarm input          : d=3"
echo "[INFO] LULU minimum_match   : ${MINIMUM_MATCH}%"
echo "[INFO] FASTA                : $FASTA"
echo "[INFO] MOTU table           : $TABLE"
echo "[INFO] Output               : $OUT_DIR"
echo "============================================================"

# ==== 1. VSEARCH Self-Comparison ====
echo "[STEP 1] Building VSEARCH match list..."

vsearch \
    --usearch_global "$FASTA" \
    --db "$FASTA" \
    --strand plus \
    --self \
    --id 0.97 \
    --iddef 1 \
    --query_cov 0.90 \
    --maxaccepts 0 \
    --maxrejects 0 \
    --threads "$THREADS" \
    --userout "$MATCHLIST" \
    --userfields query+target+id

[[ -s "$MATCHLIST" ]] || {
    echo "[ERROR] VSEARCH match list is empty"
    exit 1
}

echo "[INFO] VSEARCH match pairs: $(wc -l < "$MATCHLIST")"
echo "[INFO] First five matches:"
head -n 5 "$MATCHLIST"

# ==== 2. LULU Curation ====
echo "[STEP 2] Running LULU with minimum_match=${MINIMUM_MATCH}..."

Rscript - \
    "$TABLE" \
    "$FASTA" \
    "$MATCHLIST" \
    "$ZERO_SAMPLE_LOG" \
    "$CURATED_TABLE" \
    "$MOTU_MAP" \
    "$CURATED_FASTA" \
    "$LULU_RDS" \
    "$SUMMARY" \
    "$MINIMUM_MATCH" <<'R_EOF'

suppressPackageStartupMessages(library(lulu))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly=TRUE)

table_file <- args[1]
fasta_file <- args[2]
match_file <- args[3]
zero_sample_log <- args[4]
curated_table_file <- args[5]
motu_map_file <- args[6]
curated_fasta_file <- args[7]
lulu_rds_file <- args[8]
summary_file <- args[9]
minimum_match <- as.numeric(args[10])

# ==== Read FASTA ====
read_fasta <- function(path) {
    x <- readLines(path)
    h <- which(substr(x, 1, 1) == ">")

    if (!length(h))
        stop("No FASTA records found")

    ids <- sub("^>(\\S+).*", "\\1", x[h])
    ends <- c(h[-1] - 1, length(x))

    seqs <- mapply(
        function(a, b) paste0(x[(a + 1):b], collapse=""),
        h,
        ends,
        USE.NAMES=FALSE
    )

    setNames(toupper(seqs), ids)
}

seqs <- read_fasta(fasta_file)

if (anyDuplicated(names(seqs)))
    stop("Duplicate MOTU IDs in FASTA")

# ==== Read MOTU Table ====
x <- read.delim(
    table_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

if (ncol(x) < 2)
    stop("Invalid MOTU abundance table")

motu_ids <- x[[1]]

if (anyDuplicated(motu_ids))
    stop("Duplicate MOTU IDs in abundance table")

otu <- as.data.frame(
    x[, -1, drop=FALSE],
    check.names=FALSE,
    stringsAsFactors=FALSE
)

rownames(otu) <- motu_ids
otu[] <- lapply(otu, as.numeric)

if (anyNA(otu))
    stop("NA introduced while converting abundance values to numeric")

if (any(otu < 0))
    stop("Negative abundance detected")

# ==== FASTA and Table Consistency ====
if (!setequal(rownames(otu), names(seqs))) {
    missing_fasta <- setdiff(rownames(otu), names(seqs))
    missing_table <- setdiff(names(seqs), rownames(otu))

    if (length(missing_fasta))
        cat(
            "[ERROR] Missing from FASTA:",
            paste(missing_fasta, collapse=", "),
            "\n"
        )

    if (length(missing_table))
        cat(
            "[ERROR] Missing from table:",
            paste(missing_table, collapse=", "),
            "\n"
        )

    stop("MOTU IDs differ between FASTA and abundance table")
}

seqs <- seqs[rownames(otu)]

input_motus <- nrow(otu)
input_reads <- sum(as.matrix(otu))

# ==== Remove All-Zero Samples ====
sample_totals <- colSums(otu)
zero_samples <- names(sample_totals)[sample_totals == 0]

if (length(zero_samples)) {
    otu <- otu[, sample_totals > 0, drop=FALSE]
    writeLines(zero_samples, zero_sample_log)
} else {
    writeLines("NONE", zero_sample_log)
}

if (ncol(otu) == 0)
    stop("No non-zero samples remain")

if (sum(as.matrix(otu)) != input_reads)
    stop("Read count changed after zero-sample removal")

# ==== Read VSEARCH Match List ====
matchlist <- read.delim(
    match_file,
    header=FALSE,
    col.names=c("OTUid", "hit", "match"),
    stringsAsFactors=FALSE
)

if (!nrow(matchlist))
    stop("VSEARCH match list is empty")

matchlist$match <- as.numeric(matchlist$match)

if (anyNA(matchlist$match))
    stop("Invalid sequence identity values in match list")

unknown_query <- setdiff(matchlist$OTUid, rownames(otu))
unknown_hit <- setdiff(matchlist$hit, rownames(otu))

if (length(unknown_query) || length(unknown_hit))
    stop("Match list contains MOTU IDs absent from abundance table")

# ==== Run LULU ====
res <- lulu::lulu(
    otu,
    matchlist,
    minimum_ratio_type="min",
    minimum_ratio=1,
    minimum_match=minimum_match,
    minimum_relative_cooccurence=0.95
)

curated <- res$curated_table

if (nrow(curated) != res$curated_count)
    stop("Curated MOTU count mismatch")

if (res$curated_count + res$discarded_count != input_motus)
    stop("LULU MOTU count sanity check failed")

reads_after <- sum(as.matrix(curated))

if (reads_after != input_reads)
    stop("Read count changed during LULU")

retained_ids <- rownames(curated)

if (!all(retained_ids %in% names(seqs)))
    stop("Curated MOTU absent from input FASTA")

# ==== Export Curated MOTU Table ====
curated_out <- data.frame(
    "#OTU ID"=retained_ids,
    curated,
    check.names=FALSE
)

write.table(
    curated_out,
    curated_table_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ==== Export MOTU Map ====
otu_map_out <- data.frame(
    MOTU=rownames(res$otu_map),
    res$otu_map,
    check.names=FALSE
)

write.table(
    otu_map_out,
    motu_map_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ==== Export Curated Representative FASTA ====
con <- file(curated_fasta_file, "w")

for (id in retained_ids)
    writeLines(
        c(
            paste0(">", id),
            seqs[[id]]
        ),
        con
    )

close(con)

# ==== Save Complete LULU Result ====
saveRDS(
    res,
    lulu_rds_file
)

# ==== Summary ====
summary_df <- data.frame(
    swarm_d=3,
    minimum_match=minimum_match,
    minimum_ratio_type="min",
    minimum_ratio=1,
    minimum_relative_cooccurence=0.95,
    input_MOTUs=input_motus,
    retained_MOTUs=res$curated_count,
    merged_MOTUs=res$discarded_count,
    percent_merged=round(
        100 * res$discarded_count / input_motus,
        3
    ),
    samples=ncol(otu),
    reads_before=input_reads,
    reads_after=reads_after
)

write.table(
    summary_df,
    summary_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n============================================================\n")
cat("LULU curation summary\n")
cat("Swarm d               : 3\n")
cat("Minimum match          :", minimum_match, "%\n")
cat("Input MOTUs            :", input_motus, "\n")
cat("Retained MOTUs         :", res$curated_count, "\n")
cat("Merged MOTUs           :", res$discarded_count, "\n")
cat(
    "Percent merged        :",
    round(100 * res$discarded_count / input_motus, 3),
    "%\n"
)
cat("Samples retained       :", ncol(otu), "\n")
cat("Input reads            :", input_reads, "\n")
cat("Reads after LULU       :", reads_after, "\n")
cat(
    "Zero samples removed  :",
    if (length(zero_samples))
        paste(zero_samples, collapse=", ")
    else
        "NONE",
    "\n"
)
cat("============================================================\n")

R_EOF

# ==== 3. Final Sanity Checks ====
echo "[STEP 3] Running final sanity checks..."

[[ -s "$CURATED_TABLE" ]] || {
    echo "[ERROR] Curated MOTU table missing"
    exit 1
}

[[ -s "$MOTU_MAP" ]] || {
    echo "[ERROR] LULU MOTU map missing"
    exit 1
}

[[ -s "$CURATED_FASTA" ]] || {
    echo "[ERROR] Curated MOTU FASTA missing"
    exit 1
}

[[ -s "$LULU_RDS" ]] || {
    echo "[ERROR] LULU result RDS missing"
    exit 1
}

[[ -s "$SUMMARY" ]] || {
    echo "[ERROR] LULU summary missing"
    exit 1
}

N_CURATED_TABLE=$(
    awk 'NR > 1 {n++} END {print n+0}' "$CURATED_TABLE"
)

N_CURATED_FASTA=$(
    grep -c '^>' "$CURATED_FASTA" || true
)

if [[ "$N_CURATED_TABLE" -ne "$N_CURATED_FASTA" ]]; then
    echo "[ERROR] Curated table MOTUs != curated FASTA records"
    exit 1
fi

echo "============================================================"
echo "               LULU CURATION COMPLETE"
echo "============================================================"
echo "Swarm d             : 3"
echo "LULU minimum_match  : ${MINIMUM_MATCH}%"
echo "Final MOTUs         : $N_CURATED_TABLE"
echo "Curated table       : $CURATED_TABLE"
echo "Curated FASTA       : $CURATED_FASTA"
echo "MOTU map            : $MOTU_MAP"
echo "Summary             : $SUMMARY"
echo "============================================================"

elapsed=$(( $(date +%s) - start_time ))

printf "Total elapsed time: %dm %ds\n" \
    "$((elapsed/60))" "$((elapsed%60))"

echo "Results directory: $OUT_DIR"
