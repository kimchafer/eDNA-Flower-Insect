#!/bin/bash
#SBATCH --job-name=1_trimming
#SBATCH --output=1_trimming_%j.out
#SBATCH --error=1_trimming_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=Node6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G

# ==========================================================================================
# PURPOSE
#   Prepare paired-end eDNA reads for downstream DADA2 analysis.
#
# ANALYSIS
#   1. Decompress raw FASTQ files.
#   2. Remove Illumina adapters with cutadapt.
#   3. Remove COI primers (mlCOIintF / jgHCO2198).
#   4. Discard read pairs without the expected primers.
#   5. Standardize sample IDs (e.g., CC03_Lib_S73_L001 -> CC03).
#   6. Summarize read counts retained at each trimming step.
#
# OUTPUT
#   1_trimming/1_unzipped/       : decompressed raw FASTQ files
#   1_trimming/2_adapter_trim/   : adapter-trimmed FASTQ files
#   1_trimming/3_primer_trim/    : primer-trimmed FASTQ files for DADA2
#   1_trimming/logs/             : cutadapt reports, sample-name map, and read-count summary
# ==========================================================================================

set -euo pipefail
shopt -s nullglob

# ==== Activate Conda Environment ====
source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate /home/JihoonKim/miniconda3/envs/cutadapt

# ==== 1. Tools and Paths ====
THREADS=12

WORK="/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
INPUT="/storage2/flower_eDNA/raw_data_merged"

TR="${WORK}/1_trimming"
RAW="${TR}/1_unzipped"
RAW_TRIM="${TR}/2_adapter_trim"
RAW_PRIM="${TR}/3_primer_trim"
LOGS="${TR}/logs"

mkdir -p "$RAW" "$RAW_TRIM" "$RAW_PRIM" "$LOGS"

# ==== 2. Adapter and Primer Configuration ====
ADAPT_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPT_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

FWD="GGWACWGGWTGAACWGTWTAYCCYCC"
REV_RC="TGRTTYTTYGGDCAYCCDGARGTDTA"

revcomp_iupac() {
    local s="$1"
    local r
    r="$(echo -n "$s" | rev)"
    echo -n "$r" | tr \
        "ACGTRYSWKMBDHVNacgtryswkmbdhvn" \
        "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"
}

REV="$(revcomp_iupac "$REV_RC")"

# ------------------------------------------------------------
# [STEP 1] Raw Data Preparation
# ------------------------------------------------------------
echo "[STEP 1] Unzipping raw data..."

for gz in "${INPUT}"/*_R[12]_001.fastq.gz; do
    out="${RAW}/$(basename "${gz%.gz}")"
    [[ -s "$out" ]] || gzip -dc "$gz" > "$out"
done

# ------------------------------------------------------------
# [STEP 2] Adapter Trimming
# Quality-score filtering is performed later in DADA2.
# ------------------------------------------------------------
echo "[STEP 2] Adapter trimming..."

: > "${LOGS}/cutadapt_adapter_report.txt"

for r1 in "${RAW}"/*_R1_001.fastq; do
    s="$(basename "${r1%_R1_001.fastq}")"
    r2="${RAW}/${s}_R2_001.fastq"

    [[ -f "$r2" ]] || {
        echo "[WARNING] R2 not found for ${s}; skipping."
        continue
    }

    cutadapt -j "$THREADS" \
        -a "${ADAPT_R1}" \
        -A "${ADAPT_R2}" \
        -O 15 \
        -e 0.20 \
        --no-indels \
        --trim-n \
        --minimum-length 50 \
        -o "${RAW_TRIM}/${s}_R1_001.fastq" \
        -p "${RAW_TRIM}/${s}_R2_001.fastq" \
        "$r1" "$r2" \
        >> "${LOGS}/cutadapt_adapter_report.txt"
done

# ------------------------------------------------------------
# [STEP 3] Primer Trimming + Sample ID Standardization
# Example: CC03_Lib_S73_L001 -> CC03
# ------------------------------------------------------------
echo "[STEP 3] Primer trimming and sample ID standardization..."

: > "${LOGS}/cutadapt_primer_report.txt"

SAMPLE_MAP="${LOGS}/sample_name_map.tsv"
echo -e "original_sample_ID\tsample_ID" > "$SAMPLE_MAP"

declare -A USED_SAMPLE_IDS

for r1 in "${RAW_TRIM}"/*_R1_001.fastq; do
    s="$(basename "${r1%_R1_001.fastq}")"
    r2="${RAW_TRIM}/${s}_R2_001.fastq"

    [[ -f "$r2" ]] || {
        echo "[WARNING] R2 not found for ${s}; skipping."
        continue
    }

    # Extract two-letter plant code + two-digit sample number
    if [[ "$s" =~ ^([A-Za-z]{2}[0-9]{2}) ]]; then
        short_id="${BASH_REMATCH[1]^^}"
    else
        echo "[ERROR] Unexpected sample name: ${s}" >&2
        exit 1
    fi

    # Prevent duplicate shortened sample IDs
    if [[ -n "${USED_SAMPLE_IDS[$short_id]:-}" ]]; then
        echo "[ERROR] Duplicate shortened sample ID: ${short_id}" >&2
        exit 1
    fi

    USED_SAMPLE_IDS["$short_id"]=1

    echo -e "${s}\t${short_id}" >> "$SAMPLE_MAP"

    cutadapt -j "$THREADS" \
        -g "^${FWD}" \
        -G "^${REV}" \
        -O 15 \
        -e 0.20 \
        --no-indels \
        --match-read-wildcards \
        --discard-untrimmed \
        --trim-n \
        --max-n 0 \
        --minimum-length 50 \
        --pair-filter=both \
        -o "${RAW_PRIM}/${short_id}_R1_001.fastq" \
        -p "${RAW_PRIM}/${short_id}_R2_001.fastq" \
        "$r1" "$r2" \
        >> "${LOGS}/cutadapt_primer_report.txt"
done

# ------------------------------------------------------------
# [SUMMARY]
# ------------------------------------------------------------
echo "========================================================"
echo "                 eDNA TRIMMING COMPLETE                 "
echo "========================================================"
echo " 1. RAW UNZIPPED    : $(ls -1 "${RAW}"/*_R1_001.fastq 2>/dev/null | wc -l) pairs"
echo " 2. ADAPTER TRIMMED : $(ls -1 "${RAW_TRIM}"/*_R1_001.fastq 2>/dev/null | wc -l) pairs"
echo " 3. PRIMER TRIMMED  : $(ls -1 "${RAW_PRIM}"/*_R1_001.fastq 2>/dev/null | wc -l) pairs"
echo "========================================================"
echo "[DONE] DADA2 input dir : ${RAW_PRIM}"
echo "[INFO] Sample name map : ${SAMPLE_MAP}"
echo "========================================================"

# ------------------------------------------------------------
# [TSV] Read counts summary
# ------------------------------------------------------------
echo "[INFO] Generating read-count summary..."

SUMMARY_FILE="${LOGS}/read_counts_summary.tsv"

echo -e \
    "Sample_ID\tStep1_Raw\tStep2_Adapter\tStep3_Primer" \
    > "$SUMMARY_FILE"

get_fq_count() {
    if [[ -f "$1" ]]; then
        local lines
        lines=$(wc -l < "$1")
        echo $((lines / 4))
    else
        echo 0
    fi
}

tail -n +2 "$SAMPLE_MAP" |
while IFS=$'\t' read -r original_id short_id; do

    c1=$(get_fq_count \
        "${RAW}/${original_id}_R1_001.fastq")

    c2=$(get_fq_count \
        "${RAW_TRIM}/${original_id}_R1_001.fastq")

    c3=$(get_fq_count \
        "${RAW_PRIM}/${short_id}_R1_001.fastq")

    echo -e \
        "${short_id}\t${c1}\t${c2}\t${c3}" \
        >> "$SUMMARY_FILE"
done

echo "[INFO] Read-count summary: ${SUMMARY_FILE}"
