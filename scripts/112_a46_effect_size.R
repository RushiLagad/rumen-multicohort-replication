#!/usr/bin/env Rscript
# Amendment 46.1 and 46.3: late analysis, run after all primary results were seen.
# Usage: 112_a46_effect_size.R            CLR partial correlation per target and cohort, Freedman-Lane p (9,999 permutations, seeds 461001 on), pooling with HKSJ,
#                                          benchmark exclusions; Amendment 37 permutation p recomputed as (b + 1)/(B + 1) (46.5.1).
#        112_a46_effect_size.R ancombc2    summarise the ANCOM-BC2 Freedman-Lane fits for the methane primary targets (written by 81 with A46_FL_SEED, seeds 46101 on).
# Reads the dumps in work/a46_inputs/. Writes results/late_a46/effect_size/ only. No phenotype value is printed.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); a46_gate(root); od <- a46_dir(root, "effect_size"); mode <- c(commandArgs(TRUE), "main")[1]; stopifnot(mode %in% c("main", "ancombc2"))
B_FL <- as.integer(Sys.getenv("A46_B", "9999")); if (B_FL != 9999) cat("NOTE: A46_B =", B_FL, "(test run; the plan fixes 9,999)\n")

if (mode == "ancombc2") {
  fd <- file.path(root, "results", "late_a46", "ancombc2_fl"); dec <- file.path(fd, "permutation_decision.txt"); stopifnot("no decision file: run the ANCOM-BC2 timing step first" = file.exists(dec))
  out <- file.path(od, "ancombc2_fl_methane.tsv"); a46_fresh(out); cat(readLines(dec), sep="\n")
  R <- list()
  for (tab in c("difford_bact", "wallace_prokaryotes")) {
    o <- rd(file.path(root, "results", "methane", paste0(tab, "_targets.tsv"))); o <- o[o$status == "ok", ]
    if (tab == "difford_bact") o <- o[o$id != "M1", ]; if (tab == "wallace_prokaryotes") o <- o[o$id != "M2", ]
    fs <- list.files(fd, paste0("^", tab, "_fl4[0-9]{4}_targets\\.tsv$"), full.names=TRUE); if (!length(fs)) { cat("no Freedman-Lane fits for", tab, "\n"); next }
    N <- do.call(rbind, lapply(fs, function(f) { x <- rd(f); x[x$status == "ok", c("id", "W")] }))
    for (j in seq_len(nrow(o))) { w <- N$W[N$id == o$id[j]]; b <- sum(abs(w) >= abs(o$W[j]))
      R[[length(R) + 1]] <- data.frame(table=tab, id=o$id[j], W_obs=o$W[j], B=length(w), b=b, p_fl=(b + 1) / (length(w) + 1), sd_null_fl=sd(w), analysis=A46_LABEL) } }
  R <- do.call(rbind, R); if (is.null(R)) { cat("ANCOM-BC2 Freedman-Lane permutations were not run (see decision); nothing to summarise\n"); quit(save="no") }
  a46_write(R, out); cat("== ANCOM-BC2, Freedman-Lane permutation p (b + 1)/(B + 1), methane primary targets\n"); print(R[, -ncol(R)], digits=3, row.names=FALSE); quit(save="no") }

outs <- file.path(od, c("cohort_estimates.tsv", "pooled.tsv", "exclusion_counts.tsv", "a37_perm_p_recomputed.tsv")); a46_fresh(outs)
U <- a46_units(root); tk0 <- attr(U, "tk")
cat(A46_LABEL, "\n46.1 units:", length(U), "|", paste(sapply(U, function(u) paste0(u$prefix, " (", length(u$ids), ")")), collapse=", "), "\n")
# ---- per target and cohort: CLR OLS partial r, its interval, and the Freedman-Lane p. One seed per unit, in this fixed order.
C <- do.call(rbind, lapply(seq_along(U), function(k) { u <- U[[k]]; s <- 461000 + k; t0 <- proc.time()[[3]]
  f <- a46_fl(u$E, u$D$rx, u$D$qz, u$g, B_FL, seed=s); df <- u$D$df; tt <- r2t(f$r, df); ci <- a46_ci(f$r, df)
  cat(sprintf("  %-28s n %4d df %4d targets %2d | Freedman-Lane %d permutations within %s, seed %d, %.0f s\n", u$prefix, u$D$n, df, length(u$ids), B_FL, paste(u$perm_by, collapse=" x "), s, proc.time()[[3]] - t0))
  data.frame(arm=u$arm, id=u$ids, cohort=u$cohort, table=u$prefix, outcome=u$kind, n=u$D$n, df=df, t=tt, r=f$r, lo=ci[, "lo"], hi=ci[, "hi"], p_t=2 * pt(-abs(tt), df),
    b_fl=f$b, B_fl=f$B, failed_fl=f$failed, p_fl=f$p_fl, mcse_p_fl=f$mcse, seed_fl=s, perm_within=paste(u$perm_by, collapse=" x ")) }))
# cross-check against the second estimator on file (results/clr, and the M1 linear model in results/methane): same model, so t must agree
chk <- do.call(rbind, lapply(unique(C$table), function(p) { x <- C[C$table == p, ]
  f <- if (endsWith(p, "_M1")) file.path(root, "results", "methane", paste0(sub("_M1$", "", p), "_targets.tsv")) else file.path(root, "results", "clr", paste0(p, "_clr_targets.tsv"))
  o <- rd(f); j <- match(x$id, o$id); data.frame(table=p, file=basename(f), n=nrow(x), max_abs_dt=max(abs(x$t - o$W[j])), df_same=all(x$df == o$df[j])) }))
cat("== check against the second estimator on file (t must equal its W):\n"); print(chk, row.names=FALSE, digits=3)
stopifnot("46.1 estimates differ from the CLR / M1 estimates on file" = all(chk$max_abs_dt < 1e-6), all(chk$df_same))
j <- match(C$id, tk0$id); C$set <- tk0$set[j]; C$target <- tk0$target[j]; C$site <- tk0$site[j]; C$published_dir <- a46_dir_of(tk0, C$id, C$arm); C$analysis <- A46_LABEL
# ---- pooled per target (all testable targets of each arm)
TT <- unique(C[, c("arm", "id")]); P <- do.call(rbind, lapply(seq_len(nrow(TT)), function(i) { x <- C[C$arm == TT$arm[i] & C$id == TT$id[i], ]
  cbind(arm=TT$arm[i], id=TT$id[i], a46_pool(x$r, x$df), cohorts=paste(x$cohort, collapse=";")) }))
j <- match(P$id, tk0$id); P$set <- tk0$set[j]; P$target <- tk0$target[j]; P$site <- tk0$site[j]; P$published_dir <- a46_dir_of(tk0, P$id, P$arm)
P$stratum <- ifelse(P$arm == "Methane" & P$set == "Nelore", "Nelore RME_GC genus, tested against daily output", ifelse(P$site == "feces", "fecal source",
  ifelse(is.na(P$published_dir), "no published direction", "directed rumen")))
P$k_class <- ifelse(P$k >= 2, "two or more cohorts", "one cohort")
for (m in c("", "_fe")) { b <- a46_bound(P[[paste0("lo", m)]], P[[paste0("hi", m)]], P$published_dir); P[[paste0("bound", m)]] <- b
  P[[paste0("excl_r30", m)]] <- b < 0.30; P[[paste0("excl_r15", m)]] <- b < 0.15; P[[paste0("excl_any", m)]] <- b < 0 }
# previous (ANCOM-BC2-derived, calibrated) equivalence for comparison only
pe <- rd(file.path(root, "results", "power", "power_equivalence.tsv")); j <- match(paste(P$arm, P$id), paste(pe$arm, pe$id))
P$r_ancom <- pe$r[j]; P$lo_ancom <- pe$ci_lo[j]; P$hi_ancom <- pe$ci_hi[j]; P$excl_r30_ancom <- pe$excludes_r30[j]; P$excl_r15_ancom <- pe$excludes_r15[j]; P$label_a37 <- pe$label[j]
P$analysis <- A46_LABEL
cat("testable targets: RFI", sum(P$arm == "RFI"), "(plan: 27) | Methane", sum(P$arm == "Methane"), "(plan: 11) | missing from the ANCOM-BC2 equivalence table:", sum(is.na(j)), "\n")
# ---- exclusion counts by arm and stratum (46.5.2), HKSJ (k >= 2) or the cohort interval (k = 1); fixed effect and ANCOM-BC2 beside
cnt <- function(x, by) do.call(rbind, lapply(split(x, by), function(y) data.frame(arm=y$arm[1], group=NA, targets=nrow(y),
  excl_r30=sum(y$excl_r30), excl_r15=sum(y$excl_r15), excl_any_pubdir=sum(y$excl_any), excl_r30_fe=sum(y$excl_r30_fe), excl_r15_fe=sum(y$excl_r15_fe),
  excl_r30_ancom=sum(y$excl_r30_ancom, na.rm=TRUE), excl_r15_ancom=sum(y$excl_r15_ancom, na.rm=TRUE))))
X <- rbind(transform(cnt(P, P$arm), group="all"), transform(cnt(P, paste(P$arm, P$stratum)), group=sub("^(RFI|Methane) ", "", names(split(P, paste(P$arm, P$stratum))))),
  transform(cnt(P, paste(P$arm, P$k_class)), group=sub("^(RFI|Methane) ", "", names(split(P, paste(P$arm, P$k_class))))))
X$analysis <- A46_LABEL
# ---- 46.5.1: Amendment 37 permutation p recomputed as (b + 1)/(B + 1) from the null fits on file. The Amendment 37 labels use W_cal = W / max(1, SD_null)
# and not these p values, so no label can change; the count is computed rather than assumed: labels are recomputed from W_cal unchanged.
nulls <- list(rfi=list(est=file.path(root, "results", "replication_rfi3_calibrated", "estimates_calibrated_all.tsv"), dir=file.path(root, "results", "rfi_null")),
  methane=list(est=file.path(root, "results", "methane_pooled", "estimates_calibrated.tsv"), dir=file.path(root, "results", "methane_dryrun")))
Q <- do.call(rbind, lapply(names(nulls), function(a) { e <- rd(nulls[[a]]$est); do.call(rbind, lapply(unique(e$table), function(tab) {
  fs <- list.files(nulls[[a]]$dir, paste0("^", tab, "_seed1[0-9]{3}_targets\\.tsv$"), full.names=TRUE)
  N <- do.call(rbind, lapply(fs, function(f) { d <- rd(f); d[d$status == "ok", c("id", "W")] })); x <- e[e$table == tab, ]
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) { w <- N$W[N$id == x$id[i]]; b <- sum(abs(w) >= abs(x$W_raw[i]))
    data.frame(arm=a, table=tab, id=x$id[i], cohort=x$cohort[i], W_raw=x$W_raw[i], B=length(w), b=b, p_perm_a37=x$p_perm[i], p_perm_b1=(b + 1) / (length(w) + 1)) })) })) }))
Q$analysis <- A46_LABEL
cat("== 46.5.1: Amendment 37 permutation p values that were 0:", sum(Q$p_perm_a37 == 0, na.rm=TRUE), "| smallest p as (b + 1)/(B + 1):", signif(min(Q$p_perm_b1), 3),
  "| check b/B equals the old p:", isTRUE(all.equal(Q$b / Q$B, Q$p_perm_a37)), "| primary labels affected: 0 (labels are computed from W_cal; p_perm is reported only)\n")
a46_write(C, outs[1]); a46_write(P, outs[2]); a46_write(X, outs[3]); a46_write(Q, outs[4])
cat("\n== 46.1 POOLED (HKSJ where k >= 2; cohort interval where k = 1). Exclusion statements are conditional on the 46.2 coverage rule (113).\n")
print(P[order(P$arm, P$stratum, P$id), c("arm", "id", "target", "stratum", "k", "r", "lo", "hi", "lo_fe", "hi_fe", "I2", "bound", "excl_r30", "excl_r15", "r_ancom", "excl_r30_ancom")], digits=2, row.names=FALSE)
cat("\n== exclusion counts\n"); print(X[, -ncol(X)], row.names=FALSE)
cat("\n== Freedman-Lane p < 0.05 by cohort table (CLR):", paste(names(tapply(C$p_fl < 0.05, C$table, sum)), tapply(C$p_fl < 0.05, C$table, sum), sep="=", collapse=", "), "\n")
cat("smallest Freedman-Lane p:", signif(min(C$p_fl, na.rm=TRUE), 3), "| failed permutation fits:", sum(C$failed_fl), "\n")
