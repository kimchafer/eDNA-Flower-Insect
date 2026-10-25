#!/usr/bin/env python3

# ==========================================================================================
# PURPOSE
#   Filter DADA2-derived ASVs to retain arthropod COI sequences
#   using MMseqs2 taxonomic assignments and nucleotide best hits.
#
# ANALYSIS
#   1. Load MMseqs2 taxonomy assignments and best-hit results.
#   2. Retain ASVs assigned to Arthropoda.
#   3. Identify COI sequences from explicit COI annotations or mitochondrial genome records.
#   4. Record pass/fail decisions and reasons for every ASV.
#   5. Export the filtered Arthropoda COI ASV dataset.
#
# OUTPUT
#   3_mmseqs2_filtering/results/ASVs_filter_decisions.tsv
#       : filtering decision and supporting evidence for every ASV
#   3_mmseqs2_filtering/ASVs.filtered_COI_Arthropoda.fasta
#       : ASVs passing both Arthropoda and COI filters
# ==========================================================================================

import os
import re
import pandas as pd
from Bio import SeqIO

# ==== Paths and Output Directories ====
WORK_DIR = "/storage2/flower_eDNA/analysis/5_20260819_new_pipe"

INPUT_FASTA = f"{WORK_DIR}/2_DADA2/ASVs.nochim.fasta"

OUT_DIR = f"{WORK_DIR}/3_mmseqs2_filtering"
RESULTS_DIR = f"{OUT_DIR}/results"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

TAX_TSV = f"{RESULTS_DIR}/ASVs_taxonomy_result.tsv"
BEST_HIT_TSV = f"{RESULTS_DIR}/ASVs_best_hits.tsv"

DECISION_TSV = f"{RESULTS_DIR}/ASVs_filter_decisions.tsv"
FILTERED_FASTA = f"{OUT_DIR}/ASVs.filtered_COI_Arthropoda.fasta"

# ==== Input Checks ====
for path in [INPUT_FASTA, TAX_TSV, BEST_HIT_TSV]:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Missing or empty input file: {path}")

# ==== Load Input ASVs ====
records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
all_asvs = pd.DataFrame({"query": [r.id for r in records]})

# ==== Load MMseqs2 Taxonomy Results ====
tax_cols = ["query", "taxid", "rank", "taxname", "lineage"]

tax = pd.read_csv(
    TAX_TSV,
    sep="\t",
    header=None,
    names=tax_cols,
    usecols=[0, 1, 2, 3, 4],
    dtype=str
)

tax = tax.rename(columns={
    "taxid": "taxonomy_taxid",
    "rank": "taxonomy_rank",
    "taxname": "taxonomy_taxname",
    "lineage": "taxonomy_lineage"
})

# ==== Load MMseqs2 Best Hits ====
best = pd.read_csv(BEST_HIT_TSV, sep="\t", dtype=str)

for col in ["pident", "alnlen", "qcov", "tcov", "evalue", "bits"]:
    if col in best.columns:
        best[col] = pd.to_numeric(best[col], errors="coerce")

# ==== Define Locus Annotation Patterns ====
coi_re = re.compile(
    r"(cytochrome(?:\s+c)?\s+oxidase\s+subunit\s+(?:i|1)\b"
    r"|cytochrome(?:\s+c)?\s+oxidase\s+(?:i|1)\b"
    r"|\bcox[\s_-]?(?:i|1)\b"
    r"|\bco[\s_-]?1\b"
    r"|\bcoi\b)",
    re.IGNORECASE
)

mitogenome_re = re.compile(
    r"(mitochondrial\s+DNA,\s*(?:complete|partial)\s+(?:genome|sequence)"
    r"|(?:complete|partial)?\s*mitochondrial\s+genome"
    r"|mitochondrial\s+genome,\s*(?:complete|partial)"
    r"|genome\s+assembly[^;]*organelle:\s*mitochondrion"
    r"|organelle:\s*mitochondrion"
    r"|mitochondrion,?\s*(?:complete|partial)\s+(?:genome|sequence)"
    r"|mitochondrial\s+chromosome"
    r"|mitochondrion\s*$)",
    re.IGNORECASE
)

def locus_status(header):
    if pd.isna(header) or not str(header).strip():
        return "NO_HIT"

    header = str(header)

    # Explicit COI annotation takes precedence
    if coi_re.search(header):
        return "COI"

    if mitogenome_re.search(header):
        return "MITOGENOME"

    return "OTHER"

# ==== Merge Taxonomic and Alignment Information ====
df = all_asvs.merge(tax, on="query", how="left")
df = df.merge(best, on="query", how="left")

# ==== Apply Taxonomic and Locus Filters ====
# Taxonomic filter
df["taxon_pass"] = df["taxonomy_lineage"].str.contains(
    "Arthropoda",
    case=False,
    na=False
)

# Locus filter
df["locus_status"] = df["target_header"].apply(locus_status)
df["locus_pass"] = df["locus_status"].isin(["COI", "MITOGENOME"])

# Final decision
df["final_pass"] = df["taxon_pass"] & df["locus_pass"]

def filter_reason(row):
    if row["final_pass"]:
        return "PASS"
    if not row["taxon_pass"]:
        return "NON_ARTHROPODA"
    if row["locus_status"] == "NO_HIT":
        return "NO_HIT"
    return "NON_COI_LOCUS"

df["filter_reason"] = df.apply(filter_reason, axis=1)

# ==== Write Filtering Decision Table ====
cols = [
    "query",
    "taxonomy_taxid",
    "taxonomy_rank",
    "taxonomy_taxname",
    "taxonomy_lineage",
    "target",
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "evalue",
    "bits",
    "taxid",
    "taxname",
    "taxlineage",
    "target_header",
    "taxon_pass",
    "locus_status",
    "locus_pass",
    "final_pass",
    "filter_reason"
]

df[cols].to_csv(
    DECISION_TSV,
    sep="\t",
    index=False
)

# ==== Export Filtered FASTA ====
passed_ids = set(df.loc[df["final_pass"], "query"])
filtered_records = [r for r in records if r.id in passed_ids]

SeqIO.write(
    filtered_records,
    FILTERED_FASTA,
    "fasta"
)

# ==== Summary ====
print("=" * 60)
print("        COI + ARTHROPODA FILTERING COMPLETE")
print("=" * 60)
print(f"Total input ASVs       : {len(records)}")
print(f"ASVs with best hit     : {df['target'].notna().sum()}")
print(f"Arthropoda PASS        : {df['taxon_pass'].sum()}")

print("\nLocus status (all taxa)")
print(df["locus_status"].value_counts().to_string())

print("\nLocus status (Arthropoda only)")
print(
    df.loc[df["taxon_pass"], "locus_status"]
    .value_counts()
    .to_string()
)

print(f"\nFinal PASS             : {df['final_pass'].sum()}")
print(f"Final FAIL             : {(~df['final_pass']).sum()}")
print(f"Filtered FASTA records : {len(filtered_records)}")
print(f"\nDecision table         : {DECISION_TSV}")
print(f"Filtered FASTA         : {FILTERED_FASTA}")

print("\n[IMPORTANT QC]")
print("Mitochondrial records may use annotation expressions not covered by")
print("the current regular-expression rules.")
print("Always manually inspect locus_status = OTHER in ASVs_filter_decisions.tsv")
print("before accepting the final filtered dataset.")

print("=" * 60)

# ==== Sanity Check ====
if len(filtered_records) != int(df["final_pass"].sum()):
    raise RuntimeError(
        "Final PASS count and filtered FASTA record count differ."
    )
