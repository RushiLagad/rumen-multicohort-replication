# Replication pooling + classification (Amendments 14, 14b). Pure functions.
suppressMessages(library(metafor))
Q <- 2; ALPHA <- 0.05
dir_sign <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(grepl("^(\\+|pos|up|high|incr|more)", x), 1,
  ifelse(grepl("^(-|neg|down|low|decr|less)", x), -1, NA_real_))
}
partial_r <- function(W, df) W / sqrt(W^2 + df)
pool_one <- function(W, df, n) {
  k <- length(W); r <- partial_r(W, df); z <- atanh(r); v <- 1/(df - 1)  # = 1/(n - 3 - q), q = covariate columns (Amendment 25)
  if (k == 0) return(list(k=0, r=NA, lo=NA, hi=NA, p=NA, tau2=NA, p_hk=NA, p_fe=NA))
  if (k == 1) {
    p <- 2 * pt(-abs(W), df); s <- sqrt(v)
    return(list(k=1, r=r, lo=tanh(z-1.96*s), hi=tanh(z+1.96*s), p=p, tau2=NA, p_hk=NA, p_fe=NA))
  }
  m <- tryCatch(rma(z, v, method="REML"), error=function(e) rma(z, v, method="DL"))
  hk <- tryCatch(rma(z, v, method="REML", test="knha")$pval, error=function(e) NA)
  fe <- rma(z, v, method="FE")$pval
  list(k=k, r=tanh(m$beta[1]), lo=tanh(m$ci.lb), hi=tanh(m$ci.ub), p=m$pval,
       tau2=m$tau2, p_hk=hk, p_fe=fe)
}
# tg: data.frame(id, dir, testable) for ONE set (all rows = BH family)
# est: data.frame(id, cohort, W, df, n)
classify_set <- function(tg, est, alpha = ALPHA) {
  N <- nrow(tg); cohorts <- sort(unique(est$cohort))
  get <- function(i, drop = NULL) {
    e <- est[est$id == tg$id[i] & !(est$cohort %in% drop) & is.finite(est$W), ]
    pool_one(e$W, e$df, e$n)
  }
  P <- lapply(seq_len(N), function(i) if (tg$testable[i]) get(i) else pool_one(numeric(0),1,1))
  p <- sapply(P, function(x) if (is.na(x$p)) 1 else x$p)
  q <- p.adjust(p, "BH")
  loo_q <- sapply(cohorts, function(cc) {
    pc <- sapply(seq_len(N), function(i) {
      has <- cc %in% est$cohort[est$id == tg$id[i] & is.finite(est$W)]
      x <- if (!tg$testable[i]) NA else if (has && P[[i]]$k >= 2) get(i, cc)$p else if (has) NA else P[[i]]$p
      ifelse(is.na(x), 1, x) })
    p.adjust(pc, "BH") })
  loo_r <- sapply(cohorts, function(cc) sapply(seq_len(N), function(i) {
    if (P[[i]]$k >= 2 && cc %in% est$cohort[est$id == tg$id[i]]) get(i, cc)$r else NA }))
  out <- lapply(seq_len(N), function(i) {
    x <- P[[i]]; k <- x$k; nodir <- is.na(tg$dir[i]); s <- if (nodir) sign(x$r) else tg$dir[i]
    conc <- !is.na(x$r) && sign(x$r) == s; sig <- q[i] < alpha
    loo_ok <- NA; loo_nom <- NA
    if (k >= 2) {
      used <- cohorts[sapply(cohorts, function(cc) cc %in% est$cohort[est$id == tg$id[i] & is.finite(est$W)])]
      loo_ok <- all(sapply(used, function(cc) !is.na(loo_r[i, cc]) && sign(loo_r[i, cc]) == s && loo_q[i, cc] < alpha))
      loo_nom <- all(sapply(used, function(cc) !is.na(loo_r[i, cc]) && sign(loo_r[i, cc]) == s &&
                    get(i, cc)$p < alpha))
    }
    lab <- if (!tg$testable[i] || k == 0) "untestable"
      else if (sig && !conc) "opposite sign"
      else if (k == 1) { if (sig && conc) "single-cohort support" else "not replicated" }
      else if (sig && conc && isTRUE(loo_ok)) "replicated"
      else if (sig && conc) "single-cohort support" else "not replicated"
    if (nodir && lab != "untestable") lab <- paste(lab, "(direction not stated)")
    data.frame(id=tg$id[i], k=k, r=x$r, ci_lo=x$lo, ci_hi=x$hi, p=p[i], q_bh=q[i],
      family_size=N, tau2=x$tau2, p_hk=x$p_hk, p_fe=x$p_fe, concordant=conc,
      loo_all_sig=loo_ok, loo_nominal=loo_nom, label=lab)
  })
  do.call(rbind, out)
}
