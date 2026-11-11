#!/usr/bin/env Rscript
# Amendment 44 items 5 and 6: summary of Hedges' g (Charolais, descriptive) and REML pooling of the Li et al. diversity indices over the beef cohorts
# (two-sided, "direction not stated"; BH over the four indices). Late analysis, run after all primary results were seen. Usage: 107_late_hedges_diversity.R
suppressMessages(library(metafor)); root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")
sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
LATE_LABEL <- "late analysis, run after all primary results were seen"; rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("commit Amendment 44 first" = any(grepl("^\\*\\*Amendment 44 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
hd <- file.path(root, "results", "late_hedges_g"); hf <- list.files(hd, "_hedges_g\\.tsv$", full.names=TRUE)
if (length(hf)) { H <- do.call(rbind, lapply(hf, rd)); H <- H[H$status == "ok", ]
  cat("== Hedges' g, Charolais mRNA (", LATE_LABEL, "; descriptive, not pooled; extreme-group selection inflates g)\n", sep="")
  cat("targets", nrow(H), "| n high", H$n_high[1], "low", H$n_low[1], "| J", round(H$J[1], 4), "\n")
  cat("g unadjusted", round(min(H$g_unadjusted), 2), "to", round(max(H$g_unadjusted), 2), "(median |g|", round(median(abs(H$g_unadjusted)), 2), ") | g covariate-adjusted",
      round(min(H$g_adjusted), 2), "to", round(max(H$g_adjusted), 2), "(median |g|", round(median(abs(H$g_adjusted)), 2), ")\n")
  print(H[order(-abs(H$g_adjusted)), c("id", "target", "rank", "g_unadjusted", "g_adjusted", "r_clr", "r_ancombc2_calibrated")], row.names=FALSE, digits=2)
  write.table(H, file.path(hd, "charolais_hedges_g_all_ranks.tsv"), sep="\t", quote=FALSE, row.names=FALSE) }
dd <- file.path(root, "results", "late_diversity"); f <- file.path(dd, paste0(c("nelore", "charolais_mrna", "alberta_all"), "_diversity.tsv")); stopifnot(all(file.exists(f)))
D <- do.call(rbind, lapply(f, rd)); cat("\n== Li et al. diversity claims (Section 3.1), per cohort\n"); print(D[, c("index", "cohort", "genera_present", "depth", "n", "df", "W", "r", "status")], row.names=FALSE, digits=3)
ix <- c("bacterial_shannon_evenness", "bacterial_shannon", "archaeal_simpson", "archaeal_inverse_simpson")
e <- D[D$status == "ok", ]; e$id <- e$index; tg <- data.frame(id=ix, dir=NA_real_, testable=ix %in% e$id)
P <- cbind(index=ix, claim_in_Li=c("CHAR, p = 0.03", "HYB, p = 0.04", "CHAR, p = 0.03", "CHAR, p = 0.03"), classify_set(tg, e[, c("id", "cohort", "W", "df", "n")]))
for (cc in c("nelore", "charolais_mrna", "alberta_all")) { x <- D[D$cohort == cc, ]; j <- match(ix, x$index); P[[paste0("r_", cc)]] <- x$r[j]; P[[paste0("status_", cc)]] <- x$status[j] }
P$analysis <- LATE_LABEL; write.table(P, file.path(dd, "diversity_pooled.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n== pooled (REML over cohorts, BH over the four indices, two-sided; direction not stated)\n"); print(P[, c("index", "k", "r", "ci_lo", "ci_hi", "p", "q_bh", "tau2", "label")], row.names=FALSE, digits=3)
