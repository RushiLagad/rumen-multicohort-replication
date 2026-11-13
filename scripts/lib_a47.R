# Amendment 47: late analysis, run after all primary, late and Amendment 46 results were seen. Shared functions for 120 (47.1), 121 (47.2, 47.4) and 122 (47.3).
# Sourced after lib_a46.R (estimator, design, Freedman-Lane and unit loaders are reused unchanged). Outputs only in results/late_a47/; nothing is overwritten.
A47_LABEL <- "late analysis, run after all primary, late and Amendment 46 results were seen (Amendment 47)"
A47_METHODS <- c("hksj", "mkh", "wald", "fe")
A47_NAMES <- c(hksj="HKSJ", mkh="modified Knapp-Hartung", wald="REML-Wald", fe="fixed effect")
a47_gate <- function(root) {
  pl <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")); stopifnot("commit Amendment 47 first" = any(grepl("^\\*\\*Amendment 47 \\(", pl)))
  stopifnot("Amendment 47 never runs on shuffled data" = !nzchar(Sys.getenv("SHUFFLE_SEED")))
  ok <- tryCatch({ rma(c(0.1, 0.2), c(0.01, 0.02), method="REML", test="adhoc"); TRUE }, error=function(e) FALSE)
  stopifnot("metafor must support test = 'adhoc'" = ok) }
a47_dir <- function(root, ...) { d <- file.path(root, "results", "late_a47", ...); dir.create(d, recursive=TRUE, showWarnings=FALSE); d }
a47_fresh <- function(files) { ex <- files[file.exists(files)]; if (length(ex)) stop("output exists, Amendment 47 never overwrites: ", paste(ex, collapse=", ")) }
a47_write <- function(x, f) { stopifnot("output exists" = !file.exists(f)); x$analysis <- A47_LABEL; tmp <- paste0(f, ".tmp", Sys.getpid())
  write.table(x, tmp, sep="\t", quote=FALSE, row.names=FALSE); invisible(file.rename(tmp, f)) }

# ---------------------------------------------------------------- pooling (47.1): one REML fit (DL if REML fails), four intervals in closed form, on the z scale
# returns c(mu, mu_fe, tau2, lo_/hi_/p_ for each method, failed, fallback). failed = 1 if both REML and DL fail (then only fe is filled).
A47_POOL_COLS <- c("mu", "mu_fe", "tau2", as.vector(t(outer(c("lo_", "hi_", "p_"), A47_METHODS, paste0))), "failed", "fallback")
a47_pool <- function(z, v) {
  k <- length(z); stopifnot(k >= 2); qz <- qnorm(0.975); tq <- qt(0.975, k - 1); out <- setNames(rep(NA_real_, length(A47_POOL_COLS)), A47_POOL_COLS)
  wf <- 1 / v; mf <- sum(wf * z) / sum(wf); sf <- sqrt(1 / sum(wf))
  out[c("mu_fe", "lo_fe", "hi_fe", "p_fe")] <- c(mf, mf - qz * sf, mf + qz * sf, 2 * pnorm(-abs(mf / sf)))
  fb <- 0; f <- tryCatch(suppressWarnings(rma(z, v, method="REML", control=list(stepadj=0.5, maxiter=1000))), error=function(e) NULL)
  if (is.null(f)) { fb <- 1; f <- tryCatch(suppressWarnings(rma(z, v, method="DL")), error=function(e) NULL) }
  out["fallback"] <- fb; if (is.null(f)) { out["failed"] <- 1; return(out) }
  t2 <- f$tau2; w <- 1 / (v + t2); mu <- sum(w * z) / sum(w); se <- sqrt(1 / sum(w)); q <- sum(w * (z - mu)^2) / (k - 1)
  s_kh <- sqrt(q) * se; s_ah <- sqrt(max(1, q)) * se
  out[c("mu", "tau2", "failed")] <- c(mu, t2, 0)
  out[c("lo_wald", "hi_wald", "p_wald")] <- c(mu - qz * se, mu + qz * se, 2 * pnorm(-abs(mu / se)))
  out[c("lo_hksj", "hi_hksj", "p_hksj")] <- c(mu - tq * s_kh, mu + tq * s_kh, if (s_kh > 0) 2 * pt(-abs(mu / s_kh), k - 1) else NA)
  out[c("lo_mkh", "hi_mkh", "p_mkh")] <- c(mu - tq * s_ah, mu + tq * s_ah, 2 * pt(-abs(mu / s_ah), k - 1))
  out }
# check the closed forms against metafor itself (real data only); returns the largest absolute difference over bounds and p
a47_pool_check <- function(z, v) { p <- a47_pool(z, v); if (p["failed"] == 1) return(NA_real_)
  m <- if (p["fallback"] == 1) "DL" else "REML"; ctl <- list(stepadj=0.5, maxiter=1000)
  kh <- rma(z, v, method=m, test="knha", control=ctl); ah <- rma(z, v, method=m, test="adhoc", control=ctl); wa <- rma(z, v, method=m, test="z", control=ctl); fe <- rma(z, v, method="FE")
  max(abs(c(p["lo_hksj"] - kh$ci.lb, p["hi_hksj"] - kh$ci.ub, p["p_hksj"] - kh$pval, p["lo_mkh"] - ah$ci.lb, p["hi_mkh"] - ah$ci.ub, p["p_mkh"] - ah$pval,
    p["lo_wald"] - wa$ci.lb, p["hi_wald"] - wa$ci.ub, p["p_wald"] - wa$pval, p["lo_fe"] - fe$ci.lb, p["hi_fe"] - fe$ci.ub, p["p_fe"] - fe$pval)), na.rm=TRUE) }

# ---------------------------------------------------------------- exclusion logic (47.1; identical to a46_bound): r-scale bounds, published direction (+1, -1 or NA)
a47_excl <- function(lo, hi, dir, delta) ifelse(is.na(dir), lo > -delta & hi < delta, ifelse(dir > 0, hi < delta, -lo < delta))
a47_cp <- function(x, n) { if (n == 0) return(c(NA, NA)); as.numeric(binom.test(x, n)$conf.int) }

# ---------------------------------------------------------------- Wallace RFI NRC unit (47.4), from the dump made by 81 with A47=1 TRAIT=rfi_nrc
a47_wallace_rfi_unit <- function(root) {
  d <- a46_load(root, "wallace_prokaryotes_rfi_nrc"); C <- a46_clr(d$A); i <- match(d$feat, rownames(C)); ok <- !is.na(i)
  D <- a46_design(d$md, d$terms); stopifnot("df differs from the dumped design" = D$df == d$df)
  Y <- C[i[ok], , drop=FALSE]; rownames(Y) <- d$tk$id[ok]
  list(prefix="wallace_prokaryotes_rfi_nrc", arm="RFI", cohort="wallace", ids=d$tk$id[ok], Y=Y, E=a46_resid(Y, D), D=D,
    g=interaction(d$md[, d$perm_by, drop=FALSE], drop=TRUE), perm_by=d$perm_by, kind="CLR") }
