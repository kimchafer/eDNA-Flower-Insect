#!/usr/bin/env Rscript

# =================================================================================================================================
# PURPOSE
#   Evaluate Swarm clustering behavior across multiple d-values and identify an appropriate d-value for downstream MOTU definition.
#
# ANALYSIS
#   1. Count the number of MOTUs generated at d = 1-10.
#   2. Compare MOTU reduction between consecutive d-values.
#   3. Fit two-segment linear regressions across candidate breakpoints and identify the minimum-SSE breakpoint.
#   4. Fit an exponential-decay model to the MOTU-count curve.
#   5. Estimate the point of maximum normalized curvature.
#   6. Select the primary d-value using the piecewise breakpoint, with first differences and curvature as supporting evidence.
#
# OUTPUT
#   6_swarm_MOTU/d_selection/swarm_d_counts.tsv
#       : baseline ASV count and MOTU counts across d-values
#   6_swarm_MOTU/d_selection/swarm_d_first_difference.tsv
#       : MOTU reduction between consecutive d-values
#   6_swarm_MOTU/d_selection/swarm_d_piecewise.tsv
#       : breakpoint-specific slopes and combined SSE values
#   6_swarm_MOTU/d_selection/swarm_d_exponential.tsv
#       : exponential model parameters and curvature estimate
#   6_swarm_MOTU/d_selection/swarm_d_fitted_curve.tsv
#       : fitted exponential curve and curvature values
#   6_swarm_MOTU/d_selection/swarm_d_selection.svg
#       : observed MOTU counts, fitted curve, and diagnostic points
#   6_swarm_MOTU/d_selection/swarm_d_piecewise_SSE.svg
#       : combined SSE across candidate breakpoints
#   6_swarm_MOTU/d_selection/swarm_d_summary.txt
#       : final d-value evaluation summary
# =================================================================================================================================

# ==== Paths and Output Files ====
WORK_DIR <- "/storage2/flower_eDNA/analysis/5_20260819_new_pipe"
D_VALUES <- 1:10
OUT_DIR <- file.path(WORK_DIR, "6_swarm_MOTU", "d_selection")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

ABUND_RDS <- file.path(WORK_DIR, "5_subsequence", "ASV_abundance_resolved.rds")
INPUT_TSV <- file.path(OUT_DIR, "swarm_d_counts.tsv")
DIFF_TSV <- file.path(OUT_DIR, "swarm_d_first_difference.tsv")
PIECEWISE_TSV <- file.path(OUT_DIR, "swarm_d_piecewise.tsv")
EXP_TSV <- file.path(OUT_DIR, "swarm_d_exponential.tsv")
CURVE_TSV <- file.path(OUT_DIR, "swarm_d_fitted_curve.tsv")
SUMMARY_TXT <- file.path(OUT_DIR, "swarm_d_summary.txt")
PLOT_SVG <- file.path(OUT_DIR, "swarm_d_selection.svg")
SSE_SVG <- file.path(OUT_DIR, "swarm_d_piecewise_SSE.svg")

# ==== Load Swarm MOTU Counts ====
if (!file.exists(ABUND_RDS)) stop("Missing: ", ABUND_RDS)
abund <- as.matrix(readRDS(ABUND_RDS))
baseline_asvs <- ncol(abund)

motu_counts <- sapply(D_VALUES, function(d) {
  f <- file.path(WORK_DIR, "6_swarm_MOTU", paste0("ALL_d", d), paste0("ALL_d", d, ".swarm"))
  if (!file.exists(f)) stop("Missing Swarm result: ", f)
  sum(nzchar(readLines(f, warn=FALSE)))
})

dat <- data.frame(
  d=c(0, D_VALUES),
  label=c("ASV", paste0("d", D_VALUES)),
  count=c(baseline_asvs, motu_counts)
)
write.table(dat, INPUT_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# ==== 1. First Differences ====
first_diff <- data.frame(
  from=dat$label[-nrow(dat)],
  to=dat$label[-1],
  count_from=dat$count[-nrow(dat)],
  count_to=dat$count[-1],
  reduction=dat$count[-nrow(dat)] - dat$count[-1]
)
first_diff$percent_reduction <- 100 * first_diff$reduction / first_diff$count_from
write.table(first_diff, DIFF_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# ==== 2. Two-Segment Regression ====
candidate_bp <- D_VALUES[-length(D_VALUES)]
piecewise <- do.call(rbind, lapply(candidate_bp, function(bp) {
  left <- dat[dat$d <= bp, ]
  right <- dat[dat$d >= bp, ]
  fit_left <- lm(count ~ d, data=left)
  fit_right <- lm(count ~ d, data=right)
  sse_left <- sum(resid(fit_left)^2)
  sse_right <- sum(resid(fit_right)^2)
  data.frame(
    breakpoint=bp,
    slope_before=unname(coef(fit_left)["d"]),
    slope_after=unname(coef(fit_right)["d"]),
    SSE_before=sse_left,
    SSE_after=sse_right,
    total_SSE=sse_left + sse_right
  )
}))
best_i <- which.min(piecewise$total_SSE)
best_bp <- piecewise$breakpoint[best_i]
write.table(piecewise, PIECEWISE_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# ==== 3. Exponential-Decay Model: M(d) = c + a * exp(-b*d) ====
exp_fit <- nls(
  count ~ c + a * exp(-b*d), data=dat,
  start=list(c=min(dat$count)*0.95, a=max(dat$count)-min(dat$count), b=0.5),
  algorithm="port", lower=c(c=0, a=0, b=1e-8)
)
co <- coef(exp_fit)
pred <- predict(exp_fit)
SSE <- sum((dat$count - pred)^2)
SST <- sum((dat$count - mean(dat$count))^2)
R2 <- 1 - SSE/SST

# ==== 4. Normalized Curvature Analysis ====
# Curvature after min-max normalization of observed x and y ranges.
# This makes curvature independent of the original axis units.
x_range <- diff(range(dat$d))
y_range <- diff(range(dat$count))
curvature <- function(x) {
  fp <- -co["a"] * co["b"] * exp(-co["b"]*x)
  fpp <- co["a"] * co["b"]^2 * exp(-co["b"]*x)
  dy_du <- fp * x_range / y_range
  d2y_du2 <- fpp * x_range^2 / y_range
  abs(d2y_du2) / (1 + dy_du^2)^(3/2)
}
curv_opt <- optimize(curvature, interval=range(dat$d), maximum=TRUE)
curv_d <- curv_opt$maximum
nearest_d <- D_VALUES[which.min(abs(D_VALUES - curv_d))]

exp_stats <- data.frame(
  c=unname(co["c"]), a=unname(co["a"]), b=unname(co["b"]),
  R2=R2, SSE=SSE,
  curvature_d=curv_d,
  nearest_tested_d=nearest_d,
  max_curvature=curv_opt$objective
)
write.table(exp_stats, EXP_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# ==== 5. Export Fitted Curve ====
curve_x <- seq(min(dat$d), max(dat$d), length.out=1001)
curve_df <- data.frame(
  d=curve_x,
  fitted=co["c"] + co["a"]*exp(-co["b"]*curve_x),
  curvature=sapply(curve_x, curvature)
)
write.table(curve_df, CURVE_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# ==== 6. Generate Diagnostic Figures ====
svg(PLOT_SVG, width=7, height=5, family="Arial")
plot(dat$d, dat$count, type="b", pch=16, xaxt="n",
     xlab="Swarm d", ylab="Number of ASVs / MOTUs")
axis(1, at=dat$d, labels=c("ASV", D_VALUES))
lines(curve_df$d, curve_df$fitted, lty=2, lwd=2)
abline(v=best_bp, lty=3, lwd=2)
abline(v=curv_d, lty=4, lwd=2)
legend("topright",
       legend=c("Observed", "Exponential fit", paste0("Piecewise breakpoint: d=", best_bp),
                sprintf("Curvature maximum: d=%.2f", curv_d)),
       lty=c(1,2,3,4), pch=c(16,NA,NA,NA), bty="n")
dev.off()

svg(SSE_SVG, width=6, height=4.5, family="Arial")
plot(piecewise$breakpoint, piecewise$total_SSE, type="b", pch=16,
     xlab="Candidate breakpoint (d)", ylab="Combined SSE")
abline(v=best_bp, lty=3, lwd=2)
dev.off()

# ==== 7. Summary and d-Value Selection ====
best <- piecewise[best_i, ]
summary <- c(
  paste(rep("=", 64), collapse=""),
  "Swarm d-value evaluation",
  sprintf("Unclustered ASV baseline      : %d", baseline_asvs),
  sprintf("Tested d values               : %s", paste(D_VALUES, collapse=", ")),
  "",
  sprintf("Piecewise breakpoint          : d=%d", best_bp),
  sprintf("Slope before breakpoint       : %.4f MOTUs/d", best$slope_before),
  sprintf("Slope after breakpoint        : %.4f MOTUs/d", best$slope_after),
  sprintf("Minimum combined SSE          : %.4f", best$total_SSE),
  "",
  sprintf("Exponential model             : M(d) = %.4f + %.4f * exp(-%.6f*d)", co["c"], co["a"], co["b"]),
  sprintf("Exponential R^2               : %.6f", R2),
  sprintf("Normalized curvature maximum  : d=%.4f", curv_d),
  sprintf("Nearest tested d              : d=%d", nearest_d),
  "",
  sprintf("Primary selected d            : d=%d", best_bp),
  "Selection rule: minimum-SSE piecewise breakpoint; exponential curvature and first differences are supporting diagnostics.",
  paste(rep("=", 64), collapse="")
)
cat(paste(summary, collapse="\n"), "\n")
writeLines(summary, SUMMARY_TXT)
