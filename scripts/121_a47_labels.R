#!/usr/bin/env Rscript
# Amendment 47.2 and 47.4: late analysis, run after all primary, late and Amendment 46 results were seen.
# 47.2 (with the Amendment 47 addendum): replication labels re-derived under the 46.1 CLR estimand: classification as Amendments 14b and 14d (classify_set, lib_replication.R)
# with pooling (k >= 2, and every leave-one-out pool with k >= 2) by modified Knapp-Hartung (descriptive if 47.1 did not validate it) and the OLS t p for k = 1; BH within each target set with the Section 8 family sizes (untestable p = 1); sets,
# testability and families as 96_clr_compare.R. Check: REML-Wald pooling must reproduce the Amendment 42 CLR labels on file. Sensitivity: Freedman-Lane p (46.3) for k = 1.
# 47.4: four-cohort RFI pooling with Wallace RFI NRC under the CLR estimand (intervals from 120), labels as 47.2; depth-adjusted scan counts stated by estimator.
# Usage: 121_a47_labels.R (after 120). Writes results/late_a47/labels/ and results/late_a47/rfi4/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); a47_gate(root)
ld <- a47_dir(root, "labels"); fd <- a47_dir(root, "rfi4"); idir <- file.path(root, "results", "late_a47", "intervals")
outs <- c(file.path(ld, c("labels_rfi.tsv", "labels_methane.tsv", "label_summary.tsv", "check_a42.tsv")), file.path(fd, c("labels.tsv", "pooled.tsv", "canonical.tsv", "depth_statement.tsv"))); a47_fresh(outs)
cm <- readLines(file.path(idir, "chosen_method.txt")); chosen <- cm[length(cm)]; stopifnot(chosen %in% c(A47_METHODS, "none"))
METHOD <- if (chosen == "none") "mkh" else chosen; DESC <- chosen == "none"
cat(A47_LABEL, "| 47.2 and 47.4 | pooling method:", METHOD, if (DESC) "(not validated in 47.1: labels are descriptive)" else "(validated in 47.1)", "\n")
CV <- rd(file.path(idir, "configuration_validation.tsv")); PIv <- rd(file.path(idir, "pooled_intervals.tsv"))
# validation flags per target: the main pool's statement (47.1) and whether every leave-one-out pool with k >= 2 it uses passes at tau <= 0.10
flags <- function(version, arm, ids) { x <- PIv[PIv$version == version & PIv$arm == arm, ]; st <- x$statement[match(ids, x$id)]
  lv <- sapply(ids, function(i) { y <- CV[CV$role == paste0(version, "_loo") & CV$arm == arm & CV$id == i, ]; if (!nrow(y)) NA else all(y$pass_tau010) })
  data.frame(pooled_statement=st, loo_pools_validated=lv) }
tk0 <- rd(file.path(root, "plan", "target_keys.tsv"))

# ---------------------------------------------------------------- classification with a given pooling method (copy of classify_set with the pooling replaced)
pool_with <- function(t, df, pfl, method, kp) { k <- length(t)
  if (k == 0) return(list(k=0, r=NA, lo=NA, hi=NA, p=NA))
  r <- t / sqrt(t^2 + df)
  if (k == 1) { ci <- a46_ci(r, df); p <- if (kp == "fl") pfl else 2 * pt(-abs(t), df); return(list(k=1, r=r, lo=ci[1, "lo"], hi=ci[1, "hi"], p=p)) }
  P <- a47_pool(atanh(r), 1 / (df - 1)); mu <- if (method == "fe") P[["mu_fe"]] else P[["mu"]]
  list(k=k, r=tanh(mu), lo=tanh(P[[paste0("lo_", method)]]), hi=tanh(P[[paste0("hi_", method)]]), p=P[[paste0("p_", method)]]) }
classify_with <- function(tg, est, method, kp="t", alpha=0.05) {
  N <- nrow(tg); cohorts <- sort(unique(est$cohort))
  get <- function(i, drop=NULL) { e <- est[est$id == tg$id[i] & !(est$cohort %in% drop) & is.finite(est$W), ]; pool_with(e$W, e$df, e$p_fl, method, kp) }
  P <- lapply(seq_len(N), function(i) if (tg$testable[i]) get(i) else pool_with(numeric(0), 1, NA, method, kp))
  p <- sapply(P, function(x) if (is.na(x$p)) 1 else x$p); q <- p.adjust(p, "BH")
  has <- function(i, cc) cc %in% est$cohort[est$id == tg$id[i] & is.finite(est$W)]
  loo_q <- sapply(cohorts, function(cc) { pc <- sapply(seq_len(N), function(i) { x <- if (!tg$testable[i]) NA else if (has(i, cc) && P[[i]]$k >= 2) get(i, cc)$p else if (has(i, cc)) NA else P[[i]]$p
    ifelse(is.na(x), 1, x) }); p.adjust(pc, "BH") })
  loo_r <- sapply(cohorts, function(cc) sapply(seq_len(N), function(i) if (P[[i]]$k >= 2 && has(i, cc)) get(i, cc)$r else NA))
  loo_q <- matrix(loo_q, N, dimnames=list(NULL, cohorts)); loo_r <- matrix(loo_r, N, dimnames=list(NULL, cohorts))
  do.call(rbind, lapply(seq_len(N), function(i) { x <- P[[i]]; k <- x$k; nodir <- is.na(tg$dir[i]); s <- if (nodir) sign(x$r) else tg$dir[i]
    conc <- !is.na(x$r) && sign(x$r) == s; sig <- q[i] < alpha; loo_ok <- NA
    if (k >= 2) { used <- cohorts[sapply(cohorts, function(cc) has(i, cc))]; loo_ok <- all(sapply(used, function(cc) !is.na(loo_r[i, cc]) && sign(loo_r[i, cc]) == s && loo_q[i, cc] < alpha)) }
    lab <- if (!tg$testable[i] || k == 0) "untestable" else if (sig && !conc) "opposite sign" else if (k == 1) { if (sig && conc) "single-cohort support" else "not replicated" } else
      if (sig && conc && isTRUE(loo_ok)) "replicated" else if (sig && conc) "single-cohort support" else "not replicated"
    if (nodir && lab != "untestable") lab <- paste(lab, "(direction not stated)")
    data.frame(id=tg$id[i], k=k, r=x$r, ci_lo=x$lo, ci_hi=x$hi, p=p[i], q_bh=q[i], family_size=N, concordant=conc, loo_all_sig=loo_ok, label=lab) })) }

# ---------------------------------------------------------------- inputs: 46.1 CLR estimates (W = OLS t), Wallace RFI NRC from 120; sets and testability as 96
ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv")); ce$W <- ce$t
wr <- rd(file.path(idir, "wallace_rfi_estimates.tsv")); wr$W <- wr$t; wr$p_fl <- NA
cols <- c("id", "cohort", "W", "df", "n", "p_fl")
pR <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
pM <- function(x) ifelse(grepl("lower with higher", x), -1, ifelse(grepl("higher with higher", x), 1, NA))
tk <- tk0[!tk0$set %in% c("Methane", "MethaneSecondary"), ]
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available", ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
tm <- tk0[tk0$set %in% c("Methane", "MethaneSecondary") | (tk0$set == "Nelore" & tk0$phenotype == "RME_GC"), ]
ER <- ce[ce$arm == "RFI", cols]; EM <- ce[ce$arm == "Methane", cols]; E4 <- rbind(ER, wr[, cols])
run_rfi <- function(E, method, kp="t") do.call(rbind, lapply(unique(tk$set), function(st) { d <- tk[tk$set == st, ]
  tg <- data.frame(id=d$id, dir=pR(d$direction), testable=d$stratum %in% c("replication", "cross-site")); cbind(set=st, target=d$target, published_dir=tg$dir, classify_with(tg, E[E$id %in% d$id, ], method, kp)) }))
run_met <- function(E, method, kp="t") do.call(rbind, lapply(unique(tm$set), function(st) { d <- tm[tm$set == st, ]
  tg <- data.frame(id=d$id, dir=pM(d$direction), testable=d$status %in% c("mapped", "special")); cbind(set=st, target=d$target, published_dir=tg$dir, classify_with(tg, E[E$id %in% d$id, ], method, kp)) }))

# ---------------------------------------------------------------- check: REML-Wald pooling reproduces the Amendment 42 CLR labels
a42r <- rd(file.path(root, "results", "clr_pooled", "classification_rfi_clr.tsv")); a42m <- rd(file.path(root, "results", "clr_pooled", "classification_methane_clr.tsv"))
wR <- run_rfi(ER, "wald"); wM <- run_met(EM, "wald")
ck <- function(new, old, arm) { j <- match(paste(new$set, new$id), paste(old$set, old$id)); f <- is.finite(new$r) & is.finite(old$r[j])
  data.frame(arm=arm, targets=nrow(new), matched=sum(!is.na(j)), labels_identical=sum(new$label == old$label[j], na.rm=TRUE), max_abs_dr=max(abs(new$r[f] - old$r[j][f])),
    max_abs_dq=max(abs(new$q_bh - old$q_bh[j]), na.rm=TRUE)) }
CK <- rbind(ck(wR, a42r, "RFI"), ck(wM, a42m, "Methane")); cat("== check: REML-Wald pooling against the Amendment 42 CLR labels on file\n"); print(CK, row.names=FALSE, digits=3)
a47_write(CK, outs[4])
stopifnot("REML-Wald classification does not reproduce the Amendment 42 CLR labels" = all(CK$matched == CK$targets), all(CK$labels_identical == CK$targets), all(CK$max_abs_dr < 1e-4), all(CK$max_abs_dq < 1e-4))   # REML optimiser settings differ slightly (lib_replication default vs stepadj 0.5)

# ---------------------------------------------------------------- 47.2 labels, beside Amendment 37 (primary) and Amendment 42 (CLR, REML-Wald)
a37r <- rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")); a37m <- rd(file.path(root, "results", "methane_pooled", "classification_primary.tsv"))
lab <- function(new, fl, old37, old42, arm) { j37 <- match(paste(new$set, new$id), paste(old37$set, old37$id)); if (anyNA(j37)) j37 <- match(new$id, old37$id); j42 <- match(paste(new$set, new$id), paste(old42$set, old42$id))
  out <- data.frame(arm=arm, new[, c("set", "id", "target", "published_dir", "k", "r", "ci_lo", "ci_hi", "p", "q_bh", "family_size", "concordant", "loo_all_sig")], method=METHOD,
    descriptive=DESC, label_a47=new$label, label_a47_fl_k1=fl$label, label_a42_clr=old42$label[j42], label_a37_primary=old37$label[j37])
  out <- cbind(out, flags("primary", arm, out$id)); out$label_note <- ifelse(out$k >= 2 & (out$pooled_statement != "validated" | out$loo_pools_validated %in% FALSE),
    "configuration not validated (descriptive)", ifelse(DESC & out$k >= 2, "procedure not validated (descriptive)", ""))
  out$differs_from_a37 <- out$label_a47 != out$label_a37_primary; out }
LR <- lab(run_rfi(ER, METHOD), run_rfi(ER, METHOD, "fl"), a37r, a42r, "RFI"); LM <- lab(run_met(EM, METHOD), run_met(EM, METHOD, "fl"), a37m, a42m, "Methane")
LS <- do.call(rbind, lapply(list(LR, LM), function(x) { x <- x[!startsWith(x$label_a47, "untestable") | !startsWith(x$label_a37_primary, "untestable"), ]
  l <- sort(unique(c(x$label_a47, x$label_a37_primary, x$label_a42_clr, x$label_a47_fl_k1)))
  data.frame(arm=x$arm[1], label=l, a47=sapply(l, function(v) sum(x$label_a47 == v)), a47_fl_k1=sapply(l, function(v) sum(x$label_a47_fl_k1 == v)),
    a42_clr=sapply(l, function(v) sum(x$label_a42_clr == v, na.rm=TRUE)), a37_primary=sapply(l, function(v) sum(x$label_a37_primary == v, na.rm=TRUE))) }))
cat("\n== 47.2 labels under the CLR estimand (", METHOD, "), beside Amendment 42 and the Amendment 37 primary labels\n"); print(LS, row.names=FALSE)
for (x in list(LR, LM)) { y <- x[x$label_a47 != x$label_a37_primary & !(startsWith(x$label_a47, "untestable") & startsWith(x$label_a37_primary, "untestable")), ]
  cat(x$arm[1], "targets whose 47.2 label differs from Amendment 37:", nrow(y), "\n"); if (nrow(y)) print(y[, c("id", "target", "k", "r", "q_bh", "label_a47", "label_a42_clr", "label_a37_primary")], digits=3, row.names=FALSE) }

# ---------------------------------------------------------------- 47.4 (a) four-cohort RFI pooling under CLR
L4 <- run_rfi(E4, METHOD); a43 <- rd(file.path(root, "results", "replication_rfi4_wallace_sensitivity", "classification.tsv"))
gain <- unique(wr$id); L4 <- L4[L4$id %in% gain, ]; j3 <- match(paste(L4$set, L4$id), paste(LR$set, LR$id)); j43 <- match(paste(L4$set, L4$id), paste(a43$set, a43$id))
L4 <- data.frame(L4[, c("set", "id", "target", "published_dir", "k", "r", "ci_lo", "ci_hi", "q_bh", "label")], method=METHOD, descriptive=DESC, label_3cohort_a47=LR$label_a47[j3], label_a43_ancombc2=a43$label[j43])
L4 <- cbind(L4, flags("rfi4", "RFI", L4$id))
PI <- rd(file.path(idir, "pooled_intervals.tsv")); P4 <- PI[PI$version == "rfi4", ]; P3 <- PI[PI$version == "primary" & PI$arm == "RFI", ]; j <- match(P4$id, P3$id)
P4$r_3cohort <- P3$r[j]; P4$k_3cohort <- P3$k[j]; for (m in A47_METHODS) { P4[[paste0("excl_30_", m, "_3cohort")]] <- P3[[paste0("excl_30_", m)]][j] }
j <- match(P4$id, a43$id); P4$r_a43_ancombc2 <- a43$r[j]; P4$lo_a43_ancombc2 <- a43$ci_lo[j]; P4$hi_a43_ancombc2 <- a43$ci_hi[j]
P4$excl_30_a43_ancombc2 <- a47_excl(P4$lo_a43_ancombc2, P4$hi_a43_ancombc2, P4$published_dir, 0.30)
CAN <- P4[P4$id %in% paste0("C", 1:6), ]
CS <- data.frame(statement="canonical targets C1 to C6, four cohorts (sensitivity)", method=c(A47_METHODS, "ANCOM-BC2-derived (Amendment 43)"), chosen=c(A47_METHODS == METHOD, FALSE),
  targets=nrow(CAN), excl_r30=c(sapply(A47_METHODS, function(m) sum(CAN[[paste0("excl_30_", m)]])), sum(CAN$excl_30_a43_ancombc2)),
  excl_r15=c(sapply(A47_METHODS, function(m) sum(CAN[[paste0("excl_15_", m)]])), NA), excl_r30_3cohort=c(sapply(A47_METHODS, function(m) sum(CAN[[paste0("excl_30_", m, "_3cohort")]])), NA))
cat("\n== 47.4 (a) four-cohort RFI pooling with Wallace RFI NRC (CLR), canonical targets\n"); print(CAN[, c("id", "target", "k", "r", paste0(c("lo_", "hi_"), rep(A47_METHODS, each=2)), "r_a43_ancombc2", "lo_a43_ancombc2", "hi_a43_ancombc2")], digits=3, row.names=FALSE)
print(CS[, -1], row.names=FALSE); cat("four-cohort labels (", METHOD, "): replicated", sum(L4$label == "replicated"), "of", nrow(L4), "targets gaining Wallace\n")

# ---------------------------------------------------------------- 47.4 (b) depth adjustment by estimator (no new fit)
ss <- rd(file.path(root, "results", "late_a46", "scan", "scan_summary.tsv")); pl <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md"))
stopifnot("Amendment 41 result (c) not found as recorded" = any(grepl("Scan: 14 genera at per-genus permutation q < 0.05 (was 28)", pl, fixed=TRUE)))
DS <- data.frame(estimator=c("CLR OLS, Freedman-Lane permutation, BH within table (Amendment 46.4)", "ANCOM-BC2, per-genus permutation q (Amendment 41 result, item (c))"),
  wallace_q05_unadjusted=c(ss$n_q05[ss$table == "wallace_prokaryotes"], 28), wallace_q05_depth_adjusted=c(ss$n_q05[ss$table == "wallace_prokaryotes_depth"], 14),
  source=c("results/late_a46/scan/scan_summary.tsv", "plan/ANALYSIS_PLAN.md, Amendment 41 result"))
print(DS, row.names=FALSE)
a47_write(LR, outs[1]); a47_write(LM, outs[2]); a47_write(LS, outs[3]); a47_write(L4, outs[5]); a47_write(P4, outs[6]); a47_write(CS, outs[7]); a47_write(DS, outs[8])
