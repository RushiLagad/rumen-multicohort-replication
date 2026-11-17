#!/usr/bin/env Rscript
# Amendment 48.4: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Wallace farm random effects.
# Inputs: results/late_a46/wallace_farms/per_farm_intervals.tsv (29 features, per-farm ANCOM-BC2 lfc and SE; uncalibrated, model-based) and a per-farm CLR OLS refit
# of the same features on the Amendment 46 Wallace dump (CLR of the full table as 46.1; within each farm: parity, days in milk, CH4 z; a feature below 10% prevalence
# in a farm is not fitted there, as 46.5.4). Per feature and estimand: REML random effects (DL if REML fails, recorded), modified Knapp-Hartung 95% CI, tau2, I2, Q p,
# 95% prediction interval (t on k - 2 df, k >= 3), leave-one-farm-out (metafor leave1out), farms with |studentized residual| > 2. Reading as 48.4.
# Usage: 131_a48_farms.R. Writes results/late_a48/farms/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "farms"); outs <- file.path(od, c("per_farm_clr.tsv", "random_effects.tsv", "leave_one_farm_out.tsv", "reading.tsv")); a48_fresh(outs)
cat(A48_LABEL, "| 48.4 Wallace farm random effects\n")
pf <- rd(file.path(root, "results", "late_a46", "wallace_farms", "per_farm_intervals.tsv")); sm <- rd(file.path(root, "results", "late_a46", "wallace_farms", "summary.tsv"))
FEAT <- sm$genus; stopifnot("expected 29 features" = length(FEAT) == 29 || nzchar(Sys.getenv("A48_OUT")))
# ---------------------------------------------------------------- per-farm CLR OLS
d <- a46_load(root, "wallace_prokaryotes"); C <- a46_clr(d$A); md <- d$md; tr <- d$terms[length(d$terms)]; stopifnot("farm" %in% names(md))
miss <- setdiff(FEAT, rownames(C)); if (length(miss)) stop("features not in the Wallace CLR table: ", paste(miss, collapse=", "))
PFC <- do.call(rbind, lapply(sort(unique(as.character(md$farm))), function(fm) { i <- which(as.character(md$farm) == fm); mf <- md[i, , drop=FALSE]
  terms <- setdiff(d$terms, "farm"); Dz <- a46_design(mf, terms); prev <- rowMeans(d$A[, rownames(mf), drop=FALSE] > 0)
  do.call(rbind, lapply(FEAT, function(f) { p <- if (f %in% names(prev)) prev[[f]] else NA
    if (is.na(p) || p < 0.10) return(data.frame(genus=f, farm=fm, n=nrow(mf), df=Dz$df, r=NA, lo=NA, hi=NA, prevalence=p, fitted=FALSE))
    r <- a46_r(a46_resid(C[f, rownames(mf), drop=FALSE], Dz), Dz$rx); ci <- a46_ci(r, Dz$df)
    data.frame(genus=f, farm=fm, n=nrow(mf), df=Dz$df, r=r, lo=ci[, "lo"], hi=ci[, "hi"], prevalence=p, fitted=TRUE) })) }))
# ---------------------------------------------------------------- random effects per feature and estimand
re_one <- function(yi, vi, farms) { k <- length(yi); out <- list(k=k, mu=NA, lo=NA, hi=NA, p=NA, tau2=NA, I2=NA, Q_p=NA, pi_lo=NA, pi_hi=NA, method=NA, l1o=NULL, out_farms="")
  if (k < 2) return(out)
  f <- tryCatch(suppressWarnings(rma(yi, vi, method="REML", test="adhoc", control=list(stepadj=0.5, maxiter=1000))), error=function(e) NULL); m <- "REML"
  if (is.null(f)) { f <- tryCatch(suppressWarnings(rma(yi, vi, method="DL", test="adhoc")), error=function(e) NULL); m <- "DL" }
  if (is.null(f)) { out$method <- "failed"; return(out) }
  se <- f$se; pil <- pih <- NA; if (k >= 3) { tq <- qt(0.975, k - 2); s <- sqrt(f$tau2 + se^2); pil <- f$beta[1] - tq * s; pih <- f$beta[1] + tq * s }
  l1 <- tryCatch(leave1out(f), error=function(e) NULL); rs <- tryCatch(rstudent(f)$z, error=function(e) rep(NA, k))
  out[c("mu", "lo", "hi", "p", "tau2", "I2", "Q_p", "pi_lo", "pi_hi", "method")] <- list(f$beta[1], f$ci.lb, f$ci.ub, f$pval, f$tau2, f$I2, f$QEp, pil, pih, m)
  if (!is.null(l1)) out$l1o <- data.frame(left_out=farms, estimate=l1$estimate, lo=l1$ci.lb, hi=l1$ci.ub)
  out$out_farms <- paste(farms[which(abs(rs) > 2)], collapse=";"); out }
RE <- list(); L1 <- list()
for (f in FEAT) for (est in c("ANCOM-BC2 lfc (per farm, uncalibrated)", "CLR partial r (Fisher z)")) {
  if (startsWith(est, "ANCOM")) { x <- pf[pf$genus == f & is.finite(pf$lfc) & is.finite(pf$se), ]; yi <- x$lfc; vi <- x$se^2; fm <- x$farm } else {
    x <- PFC[PFC$genus == f & PFC$fitted, ]; yi <- atanh(x$r); vi <- 1 / (x$df - 1); fm <- x$farm }
  o <- re_one(yi, vi, fm); s <- sm[sm$genus == f, ]; full_sign <- sign(s$W_full)
  bt <- function(v) if (startsWith(est, "CLR")) tanh(v) else v
  RE[[length(RE) + 1]] <- data.frame(genus=f, estimand=est, k_farms=o$k, method=o$method, estimate=bt(o$mu), lo_mkh=bt(o$lo), hi_mkh=bt(o$hi), p_mkh=o$p, tau2=o$tau2, I2=o$I2, Q_p=o$Q_p,
    pi_lo=bt(o$pi_lo), pi_hi=bt(o$pi_hi), ci_excludes_0_full_side=!is.na(o$lo) && ((full_sign > 0 && o$lo > 0) || (full_sign < 0 && o$hi < 0)),
    pi_excludes_0=!is.na(o$pi_lo) && (o$pi_lo > 0 || o$pi_hi < 0), pi_excludes_0_full_side=!is.na(o$pi_lo) && ((full_sign > 0 && o$pi_lo > 0) || (full_sign < 0 && o$pi_hi < 0)),
    outlying_farms=o$out_farms, W_full=s$W_full, lofo_sign_kept=s$lofo_sign_kept, farms_same_sign=s$farms_same_sign, farms_fit_a46=s$farms_fit)
  if (!is.null(o$l1o)) { l <- o$l1o; l$estimate <- bt(l$estimate); l$lo <- bt(l$lo); l$hi <- bt(l$hi); L1[[length(L1) + 1]] <- cbind(genus=f, estimand=est, l,
    sign_changes=sign(l$estimate) != sign(bt(o$mu)), shift=abs(l$estimate - bt(o$mu))) } }
RE <- do.call(rbind, RE); L1 <- do.call(rbind, L1)
infl <- do.call(rbind, lapply(split(L1, paste(L1$genus, L1$estimand)), function(x) data.frame(genus=x$genus[1], estimand=x$estimand[1], most_influential_farm=x$left_out[which.max(x$shift)],
  any_omission_changes_sign=any(x$sign_changes))))
RE <- merge(RE, infl, by=c("genus", "estimand"), all.x=TRUE)
RD <- do.call(rbind, lapply(unique(RE$estimand), function(e) { x <- RE[RE$estimand == e, ]
  data.frame(estimand=e, features=nrow(x), estimable_pi=sum(!is.na(x$pi_lo)), consistent_across_farms=sum(x$pi_excludes_0), consistent_ids=paste(x$genus[x$pi_excludes_0], collapse=";"),
    share_pi_includes_0=mean(!x$pi_excludes_0[!is.na(x$pi_lo)]),
    reading=if (any(x$pi_excludes_0)) "features whose 95% prediction interval excludes 0 are 'consistent across farms'; the rest keep 'stable to leaving out any farm'" else
      "no feature is consistent across farms; the text keeps 'stable to leaving out any farm' and reports the share of prediction intervals that include 0") }))
a48_write(PFC, outs[1]); a48_write(RE, outs[2]); a48_write(L1, outs[3]); a48_write(RD, outs[4])
cat("== random effects (CLR)\n"); print(RE[RE$estimand == "CLR partial r (Fisher z)", c("genus", "k_farms", "estimate", "lo_mkh", "hi_mkh", "tau2", "I2", "pi_lo", "pi_hi", "pi_excludes_0", "most_influential_farm")], digits=3, row.names=FALSE)
cat("\n== reading\n"); print(RD, digits=3, row.names=FALSE)
