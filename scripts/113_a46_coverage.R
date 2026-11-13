#!/usr/bin/env Rscript
# Amendment 46.2: late analysis, run after all primary results were seen. Coverage of the 46.1 intervals by simulation of the CLR outcome.
# For each primary target table: t_r = trait residualised on the covariates, standardised; e = the target's CLR residual from the 46.1 model, permuted within
# design groups, standardised. x* = rho t_r + sqrt(1 - rho^2) e (rho signed in the published direction, + where none), refit the 46.1 model.
# rho in {0, 0.10, 0.20, 0.30}; 1,000 replicates, seeds 46001 to 47000 (one per replicate, shared by all tables, so replicate b of every cohort enters the
# pooled interval b). Heterogeneity: cohort rho = tanh(atanh(rho) + tau z), z ~ N(0, 1) per target and cohort, tau in {0, 0.05, 0.10}.
# Usage: 113_a46_coverage.R. Env A46_CORES (default 4); A46_REPS (test runs only). Writes results/late_a46/coverage/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); suppressMessages(library(parallel)); a46_gate(root); od <- a46_dir(root, "coverage")
REPS <- as.integer(Sys.getenv("A46_REPS", "1000")); if (REPS != 1000) cat("NOTE: A46_REPS =", REPS, "(test run; the plan fixes 1,000)\n"); CORES <- as.integer(Sys.getenv("A46_CORES", "4"))
RHO <- c(0, 0.10, 0.20, 0.30); TAU <- c(0, 0.05, 0.10)
outs <- file.path(od, c("cohort_target.tsv", "cohort_summary.tsv", "pooled_target.tsv", "pooled_summary.tsv", "pass.txt")); a46_fresh(outs)
U <- a46_units(root); tk0 <- attr(U, "tk"); cat(A46_LABEL, "| 46.2 coverage | units", length(U), "| replicates", REPS, "| cores", CORES, "\n")
for (k in seq_along(U)) { u <- U[[k]]
  R <- t(qr.resid(qr(u$D$mm), t(u$Y))); U[[k]]$es <- R / apply(R, 1, sd)                      # CLR residual from the full 46.1 model, standardised (sd is permutation invariant)
  U[[k]]$tr <- u$D$rx / sd(u$D$rx); U[[k]]$dir <- { d <- a46_dir_of(tk0, u$ids, u$arm); ifelse(is.na(d), 1, d) }
  U[[k]]$key <- paste(u$ids, u$cohort) }
keys <- sort(unique(unlist(lapply(U, `[[`, "key"))))
one <- function(b) { set.seed(46000 + b); zn <- setNames(rnorm(length(keys)), keys); out <- vector("list", length(U))
  for (k in seq_along(U)) { u <- U[[k]]; P <- a46_perm_idx(u$g, 1)[, 1]; ep <- u$es[, P, drop=FALSE]; df <- u$D$df; res <- list()
    for (tau in TAU) for (rho in RHO) { rs <- u$dir * rho; rc <- tanh(atanh(rs) + tau * zn[u$key])
      Xs <- outer(rc, u$tr) + sqrt(1 - rc^2) * ep; r <- a46_r(a46_resid(Xs, u$D), u$D$rx)
      res[[length(res) + 1]] <- data.frame(b=b, unit=u$prefix, id=u$ids, cohort=u$cohort, arm=u$arm, df=df, rho=rho, tau=tau, rho_s=rs, rho_c=rc, r=r) }
    out[[k]] <- do.call(rbind, res) }
  do.call(rbind, out) }
t0 <- proc.time()[[3]]; S <- do.call(rbind, mclapply(seq_len(REPS), one, mc.cores=CORES)); cat("simulation fits done,", round(proc.time()[[3]] - t0), "s;", nrow(S), "estimates\n")
ci <- a46_ci(S$r, S$df); S$lo <- ci[, "lo"]; S$hi <- ci[, "hi"]; S$cover <- S$lo <= S$rho_c & S$rho_c <= S$hi; S$reject <- abs(r2t(S$r, S$df)) > qt(0.975, S$df)
saveRDS(S, file.path(root, "work", "a46_coverage_replicates.rds"))
# ---- per cohort (rule a; tau = 0, where the cohort's rho is rho_s)
S0 <- S[S$tau == 0, ]
CT <- do.call(rbind, lapply(split(S0, list(S0$unit, S0$id, S0$rho), drop=TRUE), function(x) data.frame(arm=x$arm[1], unit=x$unit[1], cohort=x$cohort[1], id=x$id[1], df=x$df[1],
  rho=x$rho[1], rho_s=x$rho_s[1], reps=nrow(x), coverage=mean(x$cover), mcse=sqrt(mean(x$cover) * (1 - mean(x$cover)) / nrow(x)), mean_r=mean(x$r), bias=mean(x$r) - x$rho_s[1],
  rejection=mean(x$reject))))
CT <- CT[order(CT$arm, CT$cohort, CT$id, CT$rho), ]
CS <- do.call(rbind, lapply(split(CT, list(CT$cohort, CT$rho), drop=TRUE), function(x) data.frame(arm=x$arm[1], cohort=x$cohort[1], rho=x$rho[1], targets=nrow(x),
  coverage_mean=mean(x$coverage), coverage_min=min(x$coverage), coverage_max=max(x$coverage), bias_mean=mean(x$bias), type1_mean=if (x$rho[1] == 0) mean(x$rejection) else NA)))
CS$pass_a <- CS$coverage_mean >= 0.93 & CS$coverage_mean <= 0.97; CS <- CS[order(CS$arm, CS$cohort, CS$rho), ]
# ---- pooled per target with k >= 2: HKSJ (rule b), fixed effect and REML-Wald beside it
K <- aggregate(unit ~ id + arm, unique(S0[, c("id", "arm", "unit")]), length); K <- K[K$unit >= 2, c("id", "arm")]
Sp <- S[S$id %in% K$id, ]; grp <- split(seq_len(nrow(Sp)), list(Sp$id, Sp$rho, Sp$tau, Sp$b), drop=TRUE)
t0 <- proc.time()[[3]]
V <- do.call(rbind, mclapply(grp, function(i) { x <- Sp[i, ]; p <- suppressWarnings(a46_pool(x$r, x$df)); s <- x$rho_s[1]
  c(k=p$k, r=p$r, cover_hksj=p$lo <= s & s <= p$hi, cover_fe=p$lo_fe <= s & s <= p$hi_fe, cover_wald=p$lo_wald <= s & s <= p$hi_wald) }, mc.cores=CORES))
f1 <- Sp[vapply(grp, `[`, 1L, 1), c("id", "arm", "rho", "tau", "b", "rho_s")]; PT <- cbind(f1, as.data.frame(V)); rownames(PT) <- NULL
cat("pooled fits done,", round(proc.time()[[3]] - t0), "s;", nrow(PT), "pooled intervals\n")
m <- function(v) mean(v, na.rm=TRUE)
PTs <- do.call(rbind, lapply(split(PT, list(PT$id, PT$rho, PT$tau), drop=TRUE), function(x) data.frame(arm=x$arm[1], id=x$id[1], k=x$k[1], rho=x$rho[1], tau=x$tau[1], reps=nrow(x),
  coverage_hksj=m(x$cover_hksj), coverage_fe=m(x$cover_fe), coverage_reml_wald=m(x$cover_wald), mean_r=mean(x$r), bias=mean(x$r) - x$rho_s[1], failed=sum(is.na(x$cover_hksj)))))
PTs <- PTs[order(PTs$arm, PTs$id, PTs$tau, PTs$rho), ]
PS <- do.call(rbind, lapply(split(PTs, list(PTs$rho, PTs$tau), drop=TRUE), function(x) data.frame(rho=x$rho[1], tau=x$tau[1], targets=nrow(x), arms=paste(unique(x$arm), collapse="+"),
  coverage_hksj_mean=mean(x$coverage_hksj), coverage_hksj_min=min(x$coverage_hksj), coverage_fe_mean=mean(x$coverage_fe), coverage_reml_wald_mean=mean(x$coverage_reml_wald))))
PS$pass_b <- ifelse(PS$tau == 0, PS$coverage_hksj_mean >= 0.93, NA); PS$below_090 <- PS$coverage_hksj_mean < 0.90; PS <- PS[order(PS$tau, PS$rho), ]
for (x in c("CT", "CS", "PTs", "PS")) { y <- get(x); y$analysis <- A46_LABEL; assign(x, y) }
pa <- all(CS$pass_a); pb <- all(PS$pass_b[PS$tau == 0]); t10 <- PS$coverage_hksj_mean[PS$tau == 0.10]
txt <- c(A46_LABEL, sprintf("replicates %d per table and rho; seeds %d to %d", REPS, 46001, 46000 + REPS),
  sprintf("rule (a) per-cohort coverage averaged over targets in [0.93, 0.97] at every rho: %s (%d of %d cohort x rho cells pass)", if (pa) "PASS" else "FAIL", sum(CS$pass_a), nrow(CS)),
  sprintf("rule (b) pooled HKSJ coverage averaged over the %d targets with k >= 2 at least 0.93 at every rho (tau = 0): %s", length(unique(PT$id)), if (pb) "PASS" else "FAIL"),
  sprintf("46.2 overall: %s", if (pa && pb) "PASS (exclusion statements may use the 46.1 HKSJ intervals)" else "FAIL (intervals reported descriptively, no exclusion claim)"),
  sprintf("tau = 0.10: pooled HKSJ coverage %s; %s", paste(sprintf("%.3f", t10), collapse=", "), if (any(t10 < 0.90)) "below 0.90 at some rho, the exclusion statements say so" else "at least 0.90 at every rho"),
  sprintf("type I error of the per-cohort test at alpha 0.05 (rho = 0), by cohort: %s", paste(CS$cohort[CS$rho == 0], sprintf("%.3f", CS$type1_mean[CS$rho == 0]), collapse="; ")))
a46_write(CT, outs[1]); a46_write(CS, outs[2]); a46_write(PTs, outs[3]); a46_write(PS, outs[4]); writeLines(txt, outs[5])
cat(txt, sep="\n"); cat("\n== per cohort\n"); print(CS[, -ncol(CS)], digits=3, row.names=FALSE); cat("\n== pooled\n"); print(PS[, -ncol(PS)], digits=3, row.names=FALSE)
