#!/usr/bin/env Rscript

# ========================================================================================
# PURPOSE
#   Infer high-quality COI amplicon sequence variants (ASVs)
#   from primer-trimmed paired-end eDNA reads using DADA2.
#
# ANALYSIS
#   1. Inspect read-quality profiles.
#   2. Filter and truncate reads based on quality.
#   3. Learn forward and reverse error models.
#   4. Infer ASVs using pseudo-pooling.
#   5. Merge paired reads with strict overlap criteria.
#   6. Remove chimeric sequences.
#   7. Retain expected COI amplicon lengths.
#   8. Track read retention and ASV counts across steps.
#
# OUTPUT
#   2_DADA2/seqtab.final.rds     : final sample-by-ASV count table
#   2_DADA2/ASVs.nochim.fasta    : final ASV sequences
#   2_DADA2/filtered/             : quality-filtered FASTQ files
#   2_DADA2/plots/                : quality and error-model plots
#   2_DADA2/logs/                 : read tracking, ASV statistics, and session information
# ========================================================================================

suppressPackageStartupMessages(library(dada2))

# ==== 1. Paths and General Settings ====
WORKDIR <- "/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
IN_DIR  <- file.path(WORKDIR, "1_trimming", "3_primer_trim")
OUT_DIR <- file.path(WORKDIR, "2_DADA2")
PLOTS   <- file.path(OUT_DIR, "plots")
LOGS    <- file.path(OUT_DIR, "logs")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(LOGS, showWarnings = FALSE, recursive = TRUE)

THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "12"))
if (is.na(THREADS) || THREADS < 1) THREADS <- 12

# ==== 2. DADA2 Parameters ====
TRUNC_LEN <- c(240, 200)
MAX_EE <- c(2, 2)
TRUNC_Q <- 2
MAX_N <- 0
RM_PHIX <- TRUE

POOL_MODE <- "pseudo"
MIN_OVERLAP <- 100
MAX_MISMATCH <- 0
LEN_ALLOWED <- c(307, 310, 313, 316, 319)

# ==== 3. Input Discovery ====
fnFs <- sort(list.files(IN_DIR, pattern = "_R1_001\\.fastq$", full.names = TRUE))
fnRs <- sort(list.files(IN_DIR, pattern = "_R2_001\\.fastq$", full.names = TRUE))

if (length(fnFs) == 0 || length(fnRs) == 0)
  stop("No input FASTQ files found in: ", IN_DIR)

if (length(fnFs) != length(fnRs))
  stop("R1 and R2 file counts differ. Check: ", IN_DIR)

get_samplename <- function(x)
  sub("_R1_001\\.fastq$", "", basename(x))

sample.names <- get_samplename(fnFs)

expectedR2 <- file.path(IN_DIR, paste0(sample.names, "_R2_001.fastq"))
if (!all(expectedR2 %in% fnRs))
  stop("R1/R2 pairing mismatch")

if (anyDuplicated(sample.names))
  stop("Duplicate sample IDs detected")

valid <- grepl("^[A-Z]{2}[0-9]{2}$", sample.names)
if (!all(valid))
  stop("Unexpected sample ID format: ",
       paste(sample.names[!valid], collapse = ", "))

cat("[INFO] Samples:", length(sample.names), "\n")
cat("[INFO] truncLen(F,R):", paste(TRUNC_LEN, collapse = ","), "\n")
cat("[INFO] maxEE(F,R):", paste(MAX_EE, collapse = ","), "\n")
cat("[INFO] pool:", POOL_MODE, "\n")
cat("[INFO] mergePairs:", MIN_OVERLAP, "bp overlap /",
    MAX_MISMATCH, "mismatches\n")
cat("[INFO] allowed lengths:", paste(LEN_ALLOWED, collapse = ","), "bp\n")
cat("[INFO] threads:", THREADS, "\n\n")

# ==== 4. Quality Profiles ====
nplot <- min(6, length(fnFs))

pdf(file.path(PLOTS, "01_qualityProfile_R1.pdf"), 10, 7)
print(plotQualityProfile(fnFs[1:nplot]))
dev.off()

pdf(file.path(PLOTS, "02_qualityProfile_R2.pdf"), 10, 7)
print(plotQualityProfile(fnRs[1:nplot]))
dev.off()

# ==== 5. Quality Filtering ====
filt_path <- file.path(OUT_DIR, "filtered")
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

cat("[STEP] filterAndTrim...\n")

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = TRUNC_LEN,
  maxN = MAX_N,
  maxEE = MAX_EE,
  truncQ = TRUNC_Q,
  rm.phix = RM_PHIX,
  compress = TRUE,
  multithread = THREADS,
  verbose = TRUE
)

write.csv(out, file.path(LOGS, "filterAndTrim_readCounts.csv"), quote = FALSE)

cat(sprintf("[INFO] filter retention: %.2f%%\n",
            100 * sum(out[, "reads.out"]) / sum(out[, "reads.in"])))

# ==== 6. Error Learning ====
cat("[STEP] Learning error rates...\n")

errF <- learnErrors(filtFs, multithread = THREADS, verbose = TRUE)
errR <- learnErrors(filtRs, multithread = THREADS, verbose = TRUE)

pdf(file.path(PLOTS, "03_errorModel_F.pdf"), 10, 7)
print(plotErrors(errF, nominalQ = TRUE))
dev.off()

pdf(file.path(PLOTS, "04_errorModel_R.pdf"), 10, 7)
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

saveRDS(errF, file.path(OUT_DIR, "errF.rds"))
saveRDS(errR, file.path(OUT_DIR, "errR.rds"))

# ==== 7. Denoising ====
cat("[STEP] DADA denoising (pool=", POOL_MODE, ")...\n", sep = "")

dadaFs <- dada(filtFs, err = errF, pool = POOL_MODE, multithread = THREADS)
dadaRs <- dada(filtRs, err = errR, pool = POOL_MODE, multithread = THREADS)

saveRDS(dadaFs, file.path(OUT_DIR, "dadaFs.rds"))
saveRDS(dadaRs, file.path(OUT_DIR, "dadaRs.rds"))

# ==== 8. Pair Merging ====
cat("[STEP] Merging paired reads...\n")

mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  minOverlap = MIN_OVERLAP,
  maxMismatch = MAX_MISMATCH,
  trimOverhang = FALSE,
  verbose = TRUE
)

saveRDS(mergers, file.path(OUT_DIR, "mergers.rds"))

# ==== 9. Sequence Table Construction ====
cat("[STEP] Creating pre-chimera sequence table...\n")

seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- sample.names
saveRDS(seqtab, file.path(OUT_DIR, "seqtab.prechim.rds"))

cat("[INFO] ASVs pre-chimera:", ncol(seqtab), "\n")

len.pre <- as.data.frame(table(nchar(colnames(seqtab))))
colnames(len.pre) <- c("length", "n_ASVs")
write.csv(len.pre,
          file.path(LOGS, "ASV_length_distribution_prechim.csv"),
          row.names = FALSE, quote = FALSE)

abund.pre <- aggregate(
  colSums(seqtab),
  by = list(length = nchar(colnames(seqtab))),
  FUN = sum
)
colnames(abund.pre) <- c("length", "reads")
write.csv(abund.pre,
          file.path(LOGS, "ASV_length_distribution_prechim_weightedReads.csv"),
          row.names = FALSE, quote = FALSE)

# ==== 10. Chimera Removal ====
cat("[STEP] Removing chimeras...\n")

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = THREADS,
  verbose = TRUE
)

saveRDS(seqtab.nochim, file.path(OUT_DIR, "seqtab.nochim.rds"))
cat("[INFO] ASVs post-chimera:", ncol(seqtab.nochim), "\n")

len.post <- as.data.frame(table(nchar(colnames(seqtab.nochim))))
colnames(len.post) <- c("length", "n_ASVs")
write.csv(len.post,
          file.path(LOGS, "ASV_length_distribution_postchim.csv"),
          row.names = FALSE, quote = FALSE)

abund.post <- aggregate(
  colSums(seqtab.nochim),
  by = list(length = nchar(colnames(seqtab.nochim))),
  FUN = sum
)
colnames(abund.post) <- c("length", "reads")
write.csv(abund.post,
          file.path(LOGS, "ASV_length_distribution_postchim_weightedReads.csv"),
          row.names = FALSE, quote = FALSE)

# ==== 11. Final Length Filtering ====
cat("[STEP] Length filter:",
    paste(LEN_ALLOWED, collapse = ","), "bp\n")

keep <- nchar(colnames(seqtab.nochim)) %in% LEN_ALLOWED
seqtab.final <- seqtab.nochim[, keep, drop = FALSE]

saveRDS(seqtab.final, file.path(OUT_DIR, "seqtab.final.rds"))

cat("[INFO] ASVs before length filter :", ncol(seqtab.nochim), "\n")
cat("[INFO] ASVs after length filter  :", ncol(seqtab.final), "\n")
cat("[INFO] Reads before length filter:", sum(seqtab.nochim), "\n")
cat("[INFO] Reads after length filter :", sum(seqtab.final), "\n")

# ==== 12. Read Tracking ====
getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = out[, "reads.in"],
  filtered = out[, "reads.out"],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim),
  lenfilt = rowSums(seqtab.final)
)

rownames(track) <- sample.names
write.csv(track, file.path(LOGS, "track_reads.csv"), quote = FALSE)

summary_df <- data.frame(
  stage = c("input", "filtered", "merged", "nonchim", "lenfilt"),
  reads = c(
    sum(track[, "input"]),
    sum(track[, "filtered"]),
    sum(track[, "merged"]),
    sum(track[, "nonchim"]),
    sum(track[, "lenfilt"])
  )
)

summary_df$percent_of_input <-
  100 * summary_df$reads / summary_df$reads[1]

write.csv(
  summary_df,
  file.path(LOGS, "global_read_retention.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ==== 13. Per-Sample Final ASV Counts ====
asv_per_sample <- data.frame(
  sample = rownames(seqtab.final),
  n_ASVs = rowSums(seqtab.final > 0),
  reads_final = rowSums(seqtab.final)
)

write.csv(
  asv_per_sample,
  file.path(LOGS, "asv_counts_per_sample.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ==== 14. Export Final ASV FASTA ====
asv_seqs <- colnames(seqtab.final)
asv_headers <- paste0(">ASV_", seq_along(asv_seqs))

fasta_out <- file.path(OUT_DIR, "ASVs.nochim.fasta")
writeLines(c(rbind(asv_headers, asv_seqs)), fasta_out)

# ==== 15. Session Information and Summary ====
capture.output(sessionInfo(), file = file.path(LOGS, "sessionInfo.txt"))

cat("\n========================================================\n")
cat("[SUCCESS] DADA2 Pipeline Completed\n")
cat("Samples processed : ", length(sample.names), "\n")
cat("ASVs pre-chimera  : ", ncol(seqtab), "\n")
cat("ASVs post-chimera : ", ncol(seqtab.nochim), "\n")
cat("ASVs post-lenfilt : ", ncol(seqtab.final), "\n")
cat("Final ASV FASTA   : ", fasta_out, "\n")
cat("========================================================\n")
