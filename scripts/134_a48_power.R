#!/usr/bin/env Rscript
# Amendment 48.7: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Power, published effects and attenuation (arithmetic, no fits).
# Minimum detectable r at 80% power, two-sided alpha 0.05: per cohort, tanh(2.8016 / sqrt(df - 1)); per pooled configuration, fixed effect tanh(2.8016 / sqrt(sum(df_i - 1)))
# and the expected mKH half-width at tau = 0 (median interval width / 2 at rho = 0, tau = 0 in the 47.1 Check 1 output). Published effects: Nelore set t = coef / SE
# (plan/target_keys.tsv, Supplementary Table 7), r = t / sqrt(t^2 + 48) (df approximate); Li et al. set, extreme-group point-biserial r from the reported p (qval column)
# with t on 14 df (8 vs 8); methane sources not convertible. Each published r against each independent cohort's CLR interval. Attenuation: true r = observed r /
# sqrt(trait reliability x feature reliability), reliabilities 0.5, 0.7, 0.9. Usage: 134_a48_power.R. Writes results/late_a48/power/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "power"); outs <- file.path(od, c("mde_cohort.tsv", "mde_pooled.tsv", "published_effects.tsv", "attenuation.tsv")); a48_fresh(outs)
cat(A48_LABEL, "| 48.7 power, published effects and attenuation\n")
K <- qnorm(0.975) + qnorm(0.80); mde <- function(m) tanh(K / sqrt(m))
ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv")); PI <- rd(file.path(root, "results", "late_a47", "intervals", "pooled_intervals.tsv")); PI <- PI[PI$version == "primary", ]
c1 <- rd(file.path(root, "results", "late_a47", "intervals", "check1_cells.tsv")); tk <- rd(file.path(root, "plan", "target_keys.tsv")); dirs <- a48_dirs(root)
MC <- unique(ce[, c("cohort", "table", "n", "df")]); MC$mde_80 <- mde(MC$df - 1); MC <- MC[order(MC$cohort, MC$table), ]
P2 <- PI[PI$k >= 2, ]
MP <- do.call(rbind, lapply(seq_len(nrow(P2)), function(i) { p <- P2[i, ]; dfs <- as.numeric(strsplit(p$dfs, ";")[[1]])
  w <- c1[c1$role %in% c("primary", grep("primary", unique(c1$role), value=TRUE)) & c1$cfg == p$cfg & c1$method == "mkh" & c1$rho == 0 & c1$tau == 0, "median_width_z"]
  data.frame(id=p$id, target=p$target, k=p$k, cohorts=p$cohorts, dfs=p$dfs, mde_fe=mde(sum(dfs - 1)), mkh_t_crit=qt(0.975, p$k - 1),
    mkh_expected_halfwidth_r_tau0=if (length(w)) tanh(w[1] / 2) else NA, note=if (p$k == 2) "mKH with k = 2 uses a t on 1 df (critical value 12.7): no useful power" else "") }))
# published effects
pub <- do.call(rbind, lapply(seq_len(nrow(tk)), function(i) { x <- tk[i, ]
  if (x$set == "Li" && is.finite(x$qval)) { t <- qt(1 - x$qval / 2, 14); data.frame(id=x$id, set=x$set, source_value=paste("reported p", x$qval), r_published=t / sqrt(t^2 + 14),
    basis="extreme-group point-biserial r, two-sample t on 14 df (8 vs 8); inflated by extreme-group selection; comparable only with the Charolais design") } else
  if (x$set == "Nelore" && is.finite(x$coef) && is.finite(x$stderr) && x$stderr > 0) { t <- x$coef / x$stderr; data.frame(id=x$id, set=x$set, source_value=sprintf("coef %.4g, SE %.4g", x$coef, x$stderr),
    r_published=abs(t) / sqrt(t^2 + 48), basis=paste0("t = coef / SE from Supplementary Table 7; df 48 assumed (approximate: the source model df is not reported)",
      if (!grepl("^RFI", x$phenotype)) paste0("; published phenotype ", x$phenotype, ", tested against CH4 g/d") else "")) } else NULL }))
pub <- pub[pub$id %in% PI$id, ]   # targets with an independent test cohort
PE <- do.call(rbind, lapply(seq_len(nrow(pub)), function(i) { p <- pub[i, ]; e <- ce[ce$id == p$id, ]; pr <- PI[PI$id == p$id, ]
  if (!nrow(pr) || !nrow(e)) return(data.frame(p, cohort=NA, r_clr=NA, lo_clr=NA, hi_clr=NA, comparison="untestable (no independent cohort)"))
  e <- e[e$cohort %in% strsplit(pr$cohorts, ";")[[1]], ]; d <- dirs[[p$id]] %||% NA
  do.call(rbind, lapply(seq_len(nrow(e)), function(j) { lo <- if (is.na(d) || d > 0) e$lo[j] else -e$hi[j]; hi <- if (is.na(d) || d > 0) e$hi[j] else -e$lo[j]
    cmp <- if (is.na(d)) "no published direction" else if (hi < p$r_published) "published magnitude excluded (interval below it)" else if (lo > p$r_published) "interval above the published magnitude" else "published magnitude inside the interval"
    data.frame(p, cohort=e$cohort[j], r_clr=e$r[j], lo_clr=e$lo[j], hi_clr=e$hi[j], comparison=cmp) })) }))
PE <- rbind(PE, data.frame(id="methane targets", set="Methane, MethaneSecondary", source_value="", r_published=NA, basis="no convertible effect in the sources (group contrasts or community types)", cohort=NA, r_clr=NA, lo_clr=NA, hi_clr=NA, comparison="untestable"))
# attenuation
ex <- PI[PI$excl_15_mkh | PI$excl_30_mkh, ]; ex$bound <- a46_bound(ex$lo_mkh, ex$hi_mkh, ex$published_dir)
obs <- rbind(data.frame(quantity="minimum detectable r (80%)", what=paste(MC$cohort, MC$table), r_observed=MC$mde_80),
  data.frame(quantity="benchmark exclusion bound", what=paste(ex$id, ex$target), r_observed=ex$bound))
mv <- rd(file.path(root, "results", "late_mapping", "charolais_mrna_vs_16s_targets.tsv")); anchor <- median(mv$spearman_clr, na.rm=TRUE)
AT <- do.call(rbind, lapply(c(0.5, 0.7, 0.9), function(rt) do.call(rbind, lapply(c(0.5, 0.7, 0.9), function(rf) data.frame(obs, trait_reliability=rt, feature_reliability=rf,
  r_true=pmin(obs$r_observed / sqrt(rt * rf), 1))))))
AT$note <- paste0("true r = observed r / sqrt(trait x feature reliability), capped at 1; empirical anchor for feature reliability (lower bound only, DNA vs mRNA mixed): Charolais GTDB vs SILVA median Spearman ",
  sprintf("%.2f", anchor), " over ", sum(is.finite(mv$spearman_clr)), " targets (Amendment 45)")
a48_write(MC, outs[1]); a48_write(MP, outs[2]); a48_write(PE, outs[3]); a48_write(AT, outs[4])
cat("== minimum detectable r per cohort\n"); print(unique(MC[, c("cohort", "n", "df", "mde_80")]), digits=3, row.names=FALSE)
cat("\n== pooled configurations of the seven claims\n"); print(MP[MP$id %in% a48_seven, c("id", "k", "dfs", "mde_fe", "mkh_expected_halfwidth_r_tau0", "note")], digits=3, row.names=FALSE)
cat("\n== published effects\n"); print(PE[, c("id", "set", "r_published", "cohort", "r_clr", "lo_clr", "hi_clr", "comparison")], digits=3, row.names=FALSE)
