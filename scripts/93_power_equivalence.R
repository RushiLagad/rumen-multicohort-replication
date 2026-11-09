#!/usr/bin/env Rscript
# Per-target power and equivalence for both arms, from the calibrated pooled estimates (89 for RFI, 83 for methane). No new fits.
# Amendment 39 item 3: the Alberta all-animal analysis (classification_all) is the primary RFI classification. Usage: 93_power_equivalence.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "power"); dir.create(od, showWarnings=FALSE)
stopifnot("commit Amendment 39 first" = any(grepl("^\\*\\*Amendment 39 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
R <- rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")); R$arm <- "RFI"
M <- rd(file.path(root, "results", "methane_pooled", "classification_primary.tsv")); M$arm <- "Methane"; M$rank <- NA; M$stratum <- NA
keep <- c("arm", "set", "id", "target", "rank", "published_dir", "k", "r", "ci_lo", "ci_hi", "q_bh", "label")
D <- rbind(R[, keep], M[, keep]); D <- D[!startsWith(D$label, "untestable") & is.finite(D$r), ]
z <- atanh; za <- qnorm(0.975); zb <- qnorm(0.80)
D$se_z <- (z(D$ci_hi) - z(D$ci_lo)) / (2 * za)                      # pooled SE on Fisher's z, from the calibrated CI (random effects where k >= 2)
D$n_eff <- round(1 / D$se_z^2 + 3)                                   # the single-cohort sample size with the same precision
D$mde_r80 <- tanh((za + zb) * D$se_z)                                # minimum detectable |r| at 80% power, two-sided 0.05
D$power_r30 <- pnorm(z(0.3) / D$se_z - za)                           # power we had for a published-size effect r = 0.3
D$power_r15 <- pnorm(z(0.15) / D$se_z - za)
# equivalence: largest effect in the published direction still inside the 95% CI (two-sided); NA direction -> larger of the two bounds
D$bound_pubdir <- ifelse(is.na(D$published_dir), pmax(abs(D$ci_lo), abs(D$ci_hi)), ifelse(D$published_dir > 0, D$ci_hi, -D$ci_lo))
D$excludes_r30 <- D$bound_pubdir < 0.3; D$excludes_r15 <- D$bound_pubdir < 0.15
D$statement <- ifelse(D$label %in% c("replicated"), "replicated",
  ifelse(D$bound_pubdir < 0, "CI excludes any effect in the published direction",
  ifelse(D$excludes_r15, "effects >= 0.15 in the published direction excluded",
  ifelse(D$excludes_r30, "effects >= 0.30 in the published direction excluded", "published-size effect (r = 0.3) not excluded: underpowered"))))
D <- D[order(D$arm, D$set, D$q_bh), ]
write.table(D, file.path(od, "power_equivalence.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== PER-TARGET POWER AND EQUIVALENCE (calibrated pooled estimates; Alberta all-animal primary, Amendment 39)\n")
print(D[, c("arm", "set", "id", "target", "k", "r", "ci_lo", "ci_hi", "q_bh", "n_eff", "mde_r80", "power_r30", "bound_pubdir", "statement")], digits=2, row.names=FALSE)
cat("\n== ARM COMPARISON\n")
A <- do.call(rbind, lapply(split(D, D$arm), function(x) data.frame(arm=x$arm[1], targets_testable=nrow(x), replicated=sum(x$label == "replicated"),
  median_abs_r=round(median(abs(x$r)), 3), max_abs_r=round(max(abs(x$r)), 3), median_n_eff=median(x$n_eff), median_mde_r80=round(median(x$mde_r80), 3),
  n_power80_for_r30=sum(x$power_r30 >= 0.8), n_excl_r30=sum(x$excludes_r30), n_excl_r15=sum(x$excludes_r15), n_excl_any_pubdir=sum(x$bound_pubdir < 0, na.rm=TRUE),
  n_underpowered=sum(!x$excludes_r30 & x$label != "replicated"))))
print(A, row.names=FALSE); write.table(A, file.path(od, "arm_comparison.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n== statements by arm\n"); print(table(D$arm, D$statement))
