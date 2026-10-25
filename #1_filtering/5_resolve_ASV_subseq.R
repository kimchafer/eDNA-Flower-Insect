#!/usr/bin/env Rscript

# ==== Paths ====
WORK_DIR <- "/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
INPUT_FASTA <- file.path(WORK_DIR, "4_MACSE", "ASVs.macse_pass.fasta")
SEQTAB_RDS <- file.path(WORK_DIR, "2_DADA2", "seqtab.final.rds")
OUT_DIR <- file.path(WORK_DIR, "5_subsequence")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MATCH_OUT <- file.path(OUT_DIR, "ASV_subsequence_matches.tsv")
RESOLUTION_OUT <- file.path(OUT_DIR, "ASV_subsequence_resolution.tsv")
RESOLVED_FASTA <- file.path(OUT_DIR, "ASVs.subsequence_resolved.fasta")
ABUND_TSV <- file.path(OUT_DIR, "ASV_abundance_resolved.tsv")
ABUND_RDS <- file.path(OUT_DIR, "ASV_abundance_resolved.rds")
SUMMARY_OUT <- file.path(OUT_DIR, "ASV_subsequence_summary.txt")

if (!file.exists(INPUT_FASTA)) stop("Input FASTA not found: ", INPUT_FASTA)
if (!file.exists(SEQTAB_RDS)) stop("DADA2 seqtab not found: ", SEQTAB_RDS)

# ==== FASTA reader ====
read_fasta <- function(path) {
  x <- readLines(path)
  h <- which(substr(x, 1, 1) == ">")
  ids <- sub("^>(\\S+).*", "\\1", x[h])
  ends <- c(h[-1] - 1, length(x))
  seqs <- mapply(function(a, b) paste0(x[(a + 1):b], collapse = ""),
                 h, ends, USE.NAMES = FALSE)
  setNames(toupper(seqs), ids)
}

revcomp <- function(x) {
  paste0(rev(strsplit(chartr("ACGTRYMKWSBDHVN",
                             "TGCAYRKMWSVHDBN", x), "")[[1]]),
         collapse = "")
}

seqs <- read_fasta(INPUT_FASTA)
ids <- names(seqs)
lens <- nchar(seqs)

# ==== Recover sample-wise abundance from DADA2 ====
seqtab <- readRDS(SEQTAB_RDS)
dada_ids <- paste0("ASV_", seq_len(ncol(seqtab)))
names(dada_ids) <- colnames(seqtab)

missing_ids <- setdiff(ids, dada_ids)
if (length(missing_ids))
  stop("MACSE PASS ASVs absent from DADA2 table: ",
       paste(head(missing_ids, 10), collapse = ", "))

keep_idx <- match(ids, dada_ids)
abund <- seqtab[, keep_idx, drop = FALSE]
colnames(abund) <- ids

# Verify FASTA sequences correspond to DADA2 ASV IDs
dada_seqs <- colnames(seqtab)[keep_idx]
if (!all(unname(seqs[ids]) == dada_seqs))
  stop("FASTA sequences and DADA2 ASV IDs are inconsistent")

# ==== Find exact subsequences ====
lengths <- sort(unique(lens))
matches <- list()
k <- 0

for (child_len in lengths) {
  for (parent_len in lengths[lengths > child_len]) {
    child_ids <- ids[lens == child_len]
    parent_ids <- ids[lens == parent_len]

    for (child in child_ids) {
      child_seq <- seqs[[child]]
      child_rc <- revcomp(child_seq)

      for (parent in parent_ids) {
        parent_seq <- seqs[[parent]]
        start0 <- regexpr(child_seq, parent_seq, fixed = TRUE)[1]
        orientation <- "forward"

        if (start0 < 0) {
          start0 <- regexpr(child_rc, parent_seq, fixed = TRUE)[1]
          orientation <- "reverse_complement"
        }
        if (start0 < 0) next

        start <- start0
        end <- start + child_len - 1

        relation <- if (start == 1 && end < parent_len) {
          "3prime_trunc"
        } else if (start > 1 && end == parent_len) {
          "5prime_trunc"
        } else if (start > 1 && end < parent_len) {
          "internal"
        } else {
          "exact"
        }

        k <- k + 1
        matches[[k]] <- data.frame(
          child_ASV = child,
          child_length = child_len,
          parent_ASV = parent,
          parent_length = parent_len,
          length_difference = parent_len - child_len,
          orientation = orientation,
          start = start,
          end = end,
          relation = relation,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

match_df <- if (length(matches)) do.call(rbind, matches) else
  data.frame(
    child_ASV=character(), child_length=integer(),
    parent_ASV=character(), parent_length=integer(),
    length_difference=integer(), orientation=character(),
    start=integer(), end=integer(), relation=character()
  )

write.table(match_df, MATCH_OUT, sep="\t", quote=FALSE, row.names=FALSE)

# ==== Resolve candidate parents ====
resolution <- list()
k <- 0

for (child in unique(match_df$child_ASV)) {
  cand <- match_df[match_df$child_ASV == child, , drop=FALSE]
  terminal <- cand[cand$relation %in% c("5prime_trunc", "3prime_trunc"), , drop=FALSE]
  preferred <- if (nrow(terminal)) terminal else cand

  nearest <- preferred[
    preferred$parent_length == min(preferred$parent_length), , drop=FALSE
  ]
  action <- if (nrow(nearest) == 1) "merge" else "ambiguous"

  for (i in seq_len(nrow(nearest))) {
    k <- k + 1
    resolution[[k]] <- data.frame(
      child_ASV = nearest$child_ASV[i],
      direct_parent_ASV = nearest$parent_ASV[i],
      parent_ASV = nearest$parent_ASV[i],
      child_length = nearest$child_length[i],
      parent_length = nearest$parent_length[i],
      length_difference = nearest$length_difference[i],
      orientation = nearest$orientation[i],
      relation = nearest$relation[i],
      action = action,
      stringsAsFactors = FALSE
    )
  }
}

res_df <- if (length(resolution)) do.call(rbind, resolution) else
  data.frame(
    child_ASV=character(), direct_parent_ASV=character(),
    parent_ASV=character(), child_length=integer(),
    parent_length=integer(), length_difference=integer(),
    orientation=character(), relation=character(), action=character()
  )

# Resolve chains to final parent
merge_rows <- res_df$action == "merge"
merge_map <- setNames(res_df$parent_ASV[merge_rows],
                      res_df$child_ASV[merge_rows])

final_parent <- function(x) {
  seen <- character()
  while (x %in% names(merge_map)) {
    if (x %in% seen) stop("Circular subsequence relationship: ", x)
    seen <- c(seen, x)
    x <- merge_map[[x]]
  }
  x
}

if (any(merge_rows))
  res_df$parent_ASV[merge_rows] <-
    vapply(res_df$parent_ASV[merge_rows], final_parent, character(1))

write.table(res_df, RESOLUTION_OUT, sep="\t", quote=FALSE, row.names=FALSE)

# ==== Merge sample-wise abundance ====
reads_before <- sum(abund)
merge_children <- unique(res_df$child_ASV[res_df$action == "merge"])

if (length(merge_children)) {
  for (child in merge_children) {
    parent <- unique(res_df$parent_ASV[
      res_df$child_ASV == child & res_df$action == "merge"
    ])
    if (length(parent) != 1) stop("Invalid resolved parent for ", child)

    abund[, parent] <- abund[, parent] + abund[, child]
  }
  abund <- abund[, !colnames(abund) %in% merge_children, drop=FALSE]
}

reads_after <- sum(abund)
if (reads_before != reads_after)
  stop("Read-count mismatch after subsequence merge")

# ==== Export resolved FASTA and abundance ====
resolved_ids <- colnames(abund)

con <- file(RESOLVED_FASTA, "w")
for (id in resolved_ids)
  writeLines(c(paste0(">", id), seqs[[id]]), con)
close(con)

saveRDS(abund, ABUND_RDS)

abund_out <- data.frame(
  Sample_ID = rownames(abund),
  abund,
  check.names = FALSE
)
write.table(abund_out, ABUND_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# ==== Summary ====
n_merge <- length(merge_children)
n_ambiguous <- length(unique(res_df$child_ASV[res_df$action == "ambiguous"]))

summary <- c(
  paste(rep("=", 65), collapse=""),
  sprintf("Input ASVs                   : %d", length(ids)),
  "ASV length distribution",
  sprintf("  %d bp : %d", lengths, sapply(lengths, function(x) sum(lens == x))),
  paste(rep("-", 65), collapse=""),
  sprintf("Exact-subsequence child ASVs : %d", length(unique(match_df$child_ASV))),
  sprintf("Automatically resolved       : %d", n_merge),
  sprintf("Ambiguous                    : %d", n_ambiguous),
  sprintf("Total matching pairs         : %d", nrow(match_df)),
  sprintf("Resolved ASVs                : %d", ncol(abund)),
  sprintf("Reads before merge           : %d", reads_before),
  sprintf("Reads after merge            : %d", reads_after),
  paste(rep("-", 65), collapse=""),
  sprintf("Resolved FASTA               : %s", RESOLVED_FASTA),
  sprintf("Resolved abundance TSV       : %s", ABUND_TSV),
  sprintf("Resolved abundance RDS       : %s", ABUND_RDS),
  sprintf("Resolution table             : %s", RESOLUTION_OUT),
  paste(rep("=", 65), collapse="")
)

cat(paste(summary, collapse="\n"), "\n")
writeLines(summary, SUMMARY_OUT)
