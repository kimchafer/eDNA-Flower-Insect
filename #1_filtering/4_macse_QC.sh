#!/usr/bin/env bash
#SBATCH --job-name=4_MACSE_QC
#SBATCH --output=4_macse_%j.out
#SBATCH --error=4_macse_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=Node4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G

# ============================================================
# PURPOSE
#   Perform codon-aware quality control of arthropod COI ASVs
#   using MACSE and remove sequences containing likely
#   frameshift or premature-stop errors.
#
# ANALYSIS
#   1. Align COI ASVs with MACSE using the invertebrate
#      mitochondrial genetic code (genetic code 5).
#   2. Generate per-sequence MACSE statistics.
#   3. Detect internal frameshifts, stop codons, and deletions.
#   4. Remove ASVs containing internal frameshifts or stops.
#   5. Retain ASVs with internal deletions alone.
#   6. Export PASS/FAIL FASTA files and a per-ASV QC table.
#
# OUTPUT
#   4_MACSE/ASVs.macse_NT.fasta
#       : codon-aware nucleotide alignment
#   4_MACSE/ASVs.macse_AA.fasta
#       : translated amino-acid alignment
#   4_MACSE/ASVs.macse_stats.csv
#       : MACSE per-sequence statistics
#   4_MACSE/ASVs.macse_QC.tsv
#       : per-ASV QC classification and failure reason
#   4_MACSE/ASVs.macse_pass.fasta
#       : ASVs passing codon-aware QC
#   4_MACSE/ASVs.macse_fail.fasta
#       : ASVs failing codon-aware QC
# ============================================================

set -euo pipefail
start_time=$(date +%s)

# ==== Environment and Paths ====
source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate macse_env

WORK_DIR="/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
INPUT_FASTA="${WORK_DIR}/3_mmseqs2_filtering/ASVs.filtered_COI_Arthropoda.fasta"
OUT_DIR="${WORK_DIR}/4_MACSE"
mkdir -p "$OUT_DIR"

GENETIC_CODE=5
MACSE_NT="${OUT_DIR}/ASVs.macse_NT.fasta"
MACSE_AA="${OUT_DIR}/ASVs.macse_AA.fasta"
MACSE_STATS="${OUT_DIR}/ASVs.macse_stats.csv"
PASS_FASTA="${OUT_DIR}/ASVs.macse_pass.fasta"
FAIL_FASTA="${OUT_DIR}/ASVs.macse_fail.fasta"
QC_TSV="${OUT_DIR}/ASVs.macse_QC.tsv"
SUMMARY="${OUT_DIR}/ASVs.macse_summary.txt"

[[ -s "$INPUT_FASTA" ]] || {
    echo "[ERROR] Missing input FASTA: $INPUT_FASTA"
    exit 1
}

# ==== MACSE Command Configuration ====
if command -v macse >/dev/null 2>&1; then
    MACSE_CMD=(macse)
else
    MACSE_JAR=$(find "$CONDA_PREFIX" -type f \
        \( -name 'macse*.jar' -o -name 'MACSE*.jar' \) \
        2>/dev/null | head -n 1)

    [[ -n "$MACSE_JAR" ]] || {
        echo "[ERROR] MACSE executable or jar not found"
        exit 1
    }

    command -v java >/dev/null 2>&1 || {
        echo "[ERROR] Java not found"
        exit 1
    }

    MACSE_CMD=(java -jar "$MACSE_JAR")
fi

echo "[INFO] MACSE QC"
echo "[INFO] Input: $INPUT_FASTA"
echo "[INFO] Genetic code: $GENETIC_CODE"

# ==== 1. Codon-Aware Alignment ====
echo "[STEP 1] MACSE alignSequences..."

"${MACSE_CMD[@]}" \
    -prog alignSequences \
    -seq "$INPUT_FASTA" \
    -gc_def "$GENETIC_CODE" \
    -out_NT "$MACSE_NT" \
    -out_AA "$MACSE_AA"

[[ -s "$MACSE_NT" && -s "$MACSE_AA" ]] || {
    echo "[ERROR] MACSE alignment output missing"
    exit 1
}

# ==== 2. Per-Sequence Statistics ====
echo "[STEP 2] MACSE statistics..."

"${MACSE_CMD[@]}" \
    -prog exportAlignment \
    -align "$MACSE_NT" \
    -gc_def "$GENETIC_CODE" \
    -out_stat_per_seq "$MACSE_STATS"

[[ -s "$MACSE_STATS" ]] || {
    echo "[ERROR] MACSE statistics missing"
    exit 1
}

# ==== 3. QC Classification ====
# FAIL: internal frameshift or internal stop
# internal deletion alone is retained
echo "[STEP 3] Classifying ASVs..."

python3 - "$INPUT_FASTA" "$MACSE_STATS" \
    "$PASS_FASTA" "$FAIL_FASTA" "$QC_TSV" "$SUMMARY" <<'PY'

import sys, csv, re

input_fasta, stats_file, pass_fasta, fail_fasta, qc_tsv, summary_file = sys.argv[1:]

def read_fasta(path):
    records, order = {}, []
    name, seq = None, []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq)
                    order.append(name)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)

        if name is not None:
            records[name] = "".join(seq)
            order.append(name)

    return records, order

records, order = read_fasta(input_fasta)

with open(stats_file, newline="") as f:
    sample = f.read(4096)
    f.seek(0)

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        reader = csv.DictReader(f, dialect=dialect)
    except csv.Error:
        reader = csv.DictReader(f)

    rows = list(reader)
    fieldnames = reader.fieldnames or []

if not fieldnames:
    raise RuntimeError("Could not read MACSE statistics header")

def norm(x):
    return re.sub(r"[^a-z0-9]+", "", x.lower())

norm_fields = {norm(x): x for x in fieldnames}

def find_field(*terms):
    for raw in fieldnames:
        n = norm(raw)
        if all(term in n for term in terms):
            return raw
    return None

id_col = norm_fields.get("seqname") or find_field("seq") or find_field("name") or find_field("id")
fs_col = norm_fields.get("internalfs") or find_field("internal", "frameshift") or find_field("frameshift")
stop_col = norm_fields.get("internalstop") or find_field("internal", "stop") or find_field("stop")
del_col = norm_fields.get("internaldel") or find_field("internal", "del") or find_field("deletion")

if id_col is None or fs_col is None or stop_col is None:
    raise RuntimeError(f"Required MACSE columns not found: {fieldnames}")

print(f"[INFO] ID={id_col}, FS={fs_col}, STOP={stop_col}, DEL={del_col or 'NA'}")

def number(value):
    if value is None:
        return 0
    m = re.search(r"-?\d+(?:\.\d+)?", str(value).strip())
    return float(m.group()) if m else 0

qc = {}

for row in rows:
    seq_id = row[id_col].strip().split()[0]
    n_fs = number(row[fs_col])
    n_stop = number(row[stop_col])
    n_del = number(row[del_col]) if del_col else 0

    reasons = []
    if n_fs > 0:
        reasons.append("internal_frameshift")
    if n_stop > 0:
        reasons.append("internal_stop")

    qc[seq_id] = {
        "internal_frameshift": n_fs,
        "internal_stop": n_stop,
        "internal_del": n_del,
        "status": "FAIL" if reasons else "PASS",
        "reason": ";".join(reasons) if reasons else "none"
    }

missing = [x for x in order if x not in qc]
if missing:
    raise RuntimeError(
        f"{len(missing)} input ASVs absent from MACSE statistics. "
        f"Examples: {', '.join(missing[:10])}"
    )

with open(qc_tsv, "w") as out:
    out.write(
        "ASV_ID\tlength\tinternal_frameshift\tinternal_stop\t"
        "internal_del\tstatus\treason\n"
    )

    for seq_id in order:
        q = qc[seq_id]
        out.write(
            f"{seq_id}\t{len(records[seq_id])}\t"
            f"{q['internal_frameshift']:g}\t"
            f"{q['internal_stop']:g}\t"
            f"{q['internal_del']:g}\t"
            f"{q['status']}\t{q['reason']}\n"
        )

with open(pass_fasta, "w") as pass_out, open(fail_fasta, "w") as fail_out:
    for seq_id in order:
        out = pass_out if qc[seq_id]["status"] == "PASS" else fail_out
        out.write(f">{seq_id}\n{records[seq_id]}\n")

n_total = len(order)
n_pass = sum(qc[x]["status"] == "PASS" for x in order)
n_fail = n_total - n_pass
n_fs = sum(qc[x]["internal_frameshift"] > 0 for x in order)
n_stop = sum(qc[x]["internal_stop"] > 0 for x in order)
n_del = sum(qc[x]["internal_del"] > 0 for x in order)

summary = (
    "=" * 60 + "\n"
    f"Total input ASVs              : {n_total}\n"
    f"MACSE QC PASS                 : {n_pass}\n"
    f"MACSE QC FAIL                 : {n_fail}\n"
    f"ASVs with internal frameshift : {n_fs}\n"
    f"ASVs with internal stop       : {n_stop}\n"
    f"ASVs with internal deletion   : {n_del}\n"
    f"PASS FASTA                    : {pass_fasta}\n"
    f"FAIL FASTA                    : {fail_fasta}\n"
    f"QC table                      : {qc_tsv}\n"
    + "=" * 60 + "\n"
)

print(summary)

with open(summary_file, "w") as out:
    out.write(summary)
PY

# ==== 4. Sanity Check ====
INPUT_N=$(grep -c '^>' "$INPUT_FASTA" || true)
PASS_N=$(grep -c '^>' "$PASS_FASTA" || true)
FAIL_N=$(grep -c '^>' "$FAIL_FASTA" || true)

if [[ $((PASS_N + FAIL_N)) -ne "$INPUT_N" ]]; then
    echo "[ERROR] PASS + FAIL != input"
    exit 1
fi

elapsed=$(($(date +%s)-start_time))

echo "================ MACSE QC SUMMARY ================"
echo "Input ASVs : $INPUT_N"
echo "PASS ASVs  : $PASS_N"
echo "FAIL ASVs  : $FAIL_N"
echo "QC table   : $QC_TSV"
echo "PASS FASTA : $PASS_FASTA"
echo "FAIL FASTA : $FAIL_FASTA"
printf "Elapsed     : %dm %ds\n" "$((elapsed/60))" "$((elapsed%60))"
echo "=================================================="
