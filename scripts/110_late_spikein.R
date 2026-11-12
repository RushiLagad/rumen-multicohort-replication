#!/usr/bin/env Rscript
# Amendment 45 item 2: spike-in summary. Planted effects (+-b natural-log fold per SD of the trait, 10% of eligible genera, half up and half down) in real
# tables whose trait was first permuted within design groups; ANCOM-BC2 refitted. Recovery and false-positive rates at |W| > 1.96, raw and calibrated
# (W / max(1, per-genus null SD) from the existing permutation tables, as in Amendment 37; table-level median null SD where a genus has no null).
# Late analysis, run after all primary results were seen. Usage: 110_late_spikein.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "late_spikein")
stopifnot("commit Amendment 45 first" = any(grepl("^\\*\\*Amendment 45 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
LATE_LABEL <- "late analysis, run after all primary results were seen"; rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
nulls <- list(wallace_prokaryotes=list(pat="results/methane_dryrun/wallace_prokaryotes_seed%d_ancombc2.csv", seeds=1001:1200, term="CH4_z"),
              alberta_genus=list(pat="results/scan_null/alberta_genus_seed%d_ancombc2.csv", seeds=2001:2200, term="RFI_z"),
              nelore_genus=list(pat="results/scan_null/nelore_genus_seed%d_ancombc2.csv", seeds=2001:2200, term="RFI_z"))
out <- list(); sds <- list()
for (tab in names(nulls)) { nl <- nulls[[tab]]; fs <- sprintf(file.path(root, nl$pat), nl$seeds); fs <- fs[file.exists(fs)]
  if (!length(fs) || !dir.exists(file.path(od, tab))) { cat("skip", tab, "| null files", length(fs), "\n"); next }
  N <- do.call(rbind, lapply(fs, function(f) { x <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE); data.frame(genus=x$taxon_name, W=x[[paste0("lfc_", nl$term)]] / x[[paste0("se_", nl$term)]]) }))
  N <- N[is.finite(N$W) & !grepl("^TGT_", N$genus), ]; sd_g <- tapply(N$W, N$genus, sd); sd_tab <- median(sd_g, na.rm=TRUE)
  sds[[tab]] <- data.frame(table=tab, null_permutations=length(fs), genera_with_null=length(sd_g), median_null_sd=round(sd_tab, 2))
  S <- do.call(rbind, lapply(list.files(file.path(od, tab), "^b.*_seed[0-9]+\\.tsv$", full.names=TRUE), rd)); S <- S[is.finite(S$W), ]
  S$sd_null <- unname(sd_g[S$genus]); S$sd_null[is.na(S$sd_null)] <- sd_tab; S$W_cal <- S$W / pmax(1, S$sd_null); S$W_cal_tab <- S$W / max(1, sd_tab)
  S$planted <- S$planted_lfc_per_sd != 0; S$dir_ok <- sign(S$W) == sign(S$planted_lfc_per_sd)
  for (b in sort(unique(S$b))) { x <- S[S$b == b, ]; p <- x[x$planted, ]; u <- x[!x$planted, ]
    out[[paste(tab, b)]] <- data.frame(table=tab, b_lfc_per_sd=b, fits=length(unique(x$seed)), planted=nrow(p), unplanted=nrow(u),
      power_raw=round(mean(abs(p$W) > 1.96 & p$dir_ok), 3), power_calibrated=round(mean(abs(p$W_cal) > 1.96 & p$dir_ok), 3), power_calibrated_table_sd=round(mean(abs(p$W_cal_tab) > 1.96 & p$dir_ok), 3),
      fpr_raw=round(mean(abs(u$W) > 1.96), 3), fpr_calibrated=round(mean(abs(u$W_cal) > 1.96), 3),
      median_est_over_planted=round(median(p$lfc / p$planted_lfc_per_sd), 2), median_r_raw=round(median(abs(p$W) / sqrt(p$W^2 + p$df)), 3), median_r_calibrated=round(median(abs(p$W_cal) / sqrt(p$W_cal^2 + p$df)), 3),
      df=x$df[1], analysis=LATE_LABEL) } }
R <- do.call(rbind, out); Sd <- do.call(rbind, sds)
write.table(R, file.path(od, "spikein_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE); write.table(Sd, file.path(od, "spikein_null_sd.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== Amendment 45 item 2: spike-in (", LATE_LABEL, ")\n", sep=""); print(Sd, row.names=FALSE); print(R[, -ncol(R)], row.names=FALSE)
cat("power = share of planted genera with |W| > 1.96 in the planted direction; fpr = share of unplanted genera with |W| > 1.96; median_est_over_planted = ANCOM-BC2 lfc / planted lfc (1 = unbiased)\n")
