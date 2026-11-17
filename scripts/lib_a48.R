# Amendment 48: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Shared functions for 128 to 136.
# Sourced after lib_a46.R and lib_a47.R (estimator, design, Freedman-Lane, unit loaders, mKH pooling and exclusion logic are reused unchanged).
# Outputs only under results/late_a48/ (A48_OUT overrides the base, for test runs only); nothing is overwritten.
A48_LABEL <- "late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen (Amendment 48)"
a48_gate <- function(root) {
  pl <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")); stopifnot("commit Amendment 48 first" = any(grepl("^\\*\\*Amendment 48 \\(", pl)))
  stopifnot("Amendment 48 never runs on shuffled data" = !nzchar(Sys.getenv("SHUFFLE_SEED"))) }
a48_base <- function(root) { b <- Sys.getenv("A48_OUT"); if (nzchar(b)) { cat("NOTE: A48_OUT =", b, "(test run; not a result)\n"); b } else file.path(root, "results", "late_a48") }
a48_dir <- function(root, ...) { d <- file.path(a48_base(root), ...); dir.create(d, recursive=TRUE, showWarnings=FALSE); d }
a48_fresh <- function(files) { ex <- files[file.exists(files)]; if (length(ex)) stop("output exists, Amendment 48 never overwrites: ", paste(ex, collapse=", ")) }
a48_write <- function(x, f) { stopifnot("output exists" = !file.exists(f)); x$analysis <- A48_LABEL; tmp <- paste0(f, ".tmp", Sys.getpid())
  write.table(x, tmp, sep="\t", quote=FALSE, row.names=FALSE); invisible(file.rename(tmp, f)) }
a48_seven <- c("T2", "T4", "T5", "T6", "N15", "S1", "M1")

# ---------------------------------------------------------------- CLR variants (48.5). zero: "pc" (pseudo-count pc) or "czm" (zCompositions::cmultRepl, CZM);
# ref: R0 = non-target features after the 10% prevalence filter (= a46_clr), R1 = every filtered feature, R2 = non-target features present in >= 90% of animals.
a48_clr <- function(A, zero="pc", pc=0.5, ref="R0", drop_ids=character(0)) {
  B <- A[rowMeans(A > 0) >= 0.10, , drop=FALSE]; tg <- grepl("^TGT_", rownames(B))
  L <- if (zero == "pc") log(B + pc) else {
    P <- suppressMessages(suppressWarnings(zCompositions::cmultRepl(t(B), method="CZM", output="prop", z.warning=1, z.delete=FALSE, suppress.print=TRUE)))
    P <- t(as.matrix(P)); stopifnot("cmultRepl changed the table shape" = identical(dim(P), dim(B)), all(is.finite(P)), all(P > 0)); dimnames(P) <- dimnames(B); log(P) }
  rf <- switch(ref, R0=!tg, R1=rep(TRUE, nrow(B)), R2=!tg & rowMeans(B > 0) >= 0.90)
  stopifnot("empty reference" = sum(rf) >= 1); out <- sweep(L, 2, colMeans(L[rf, , drop=FALSE])); attr(out, "n_ref") <- sum(rf); out }
A48_GRID <- expand.grid(zero=c("pc0.1", "pc0.5", "pc1", "czm"), ref=c("R0", "R1", "R2"), stringsAsFactors=FALSE)
a48_clr_cell <- function(A, zero, ref) if (zero == "czm") a48_clr(A, "czm", ref=ref) else a48_clr(A, "pc", as.numeric(sub("pc", "", zero)), ref)

# ---------------------------------------------------------------- estimation units with a given CLR (copy of a46_units with the transform as an argument; M1 kept as 46.1)
a48_units <- function(root, clrfun=a46_clr, drop_animals=NULL) {
  U <- list()
  add <- function(prefix, arm, cohort, ids, Y, md, terms, perm_by, kind, nref=NA) {
    D <- a46_design(md, terms); g <- interaction(md[, perm_by, drop=FALSE], drop=TRUE); rownames(Y) <- ids
    U[[length(U) + 1]] <<- list(prefix=prefix, arm=arm, cohort=cohort, ids=ids, Y=Y, E=a46_resid(Y, D), D=D, g=g, perm_by=perm_by, kind=kind, n_ref=nref) }
  for (p in c(A46_RFI, A46_MET)) { d <- a46_load(root, p); C <- clrfun(d$A); i <- match(d$feat, rownames(C)); ok <- !is.na(i)
    stopifnot("df differs from the dumped design" = a46_design(d$md, d$terms)$df == d$df)
    if (p %in% A46_RFI) { coh <- d$cohort; arm <- "RFI" } else {
      coh <- if (startsWith(p, "difford")) "difford" else "wallace"; arm <- "Methane"
      drop <- d$tk$id %in% c("M1", if (d$difford) "M3", if (coh == "wallace") "M2"); if (p == "difford_arch") drop <- rep(TRUE, length(ok)); ok <- ok & !drop }
    if (any(ok)) add(p, arm, coh, d$tk$id[ok], C[i[ok], , drop=FALSE], d$md, d$terms, d$perm_by, "CLR", attr(C, "n_ref") %||% NA)
    if (!is.null(d$m1)) add(paste0(p, "_M1"), "Methane", coh, "M1", matrix(d$m1$y, 1), d$m1$md, d$terms, d$perm_by, "M1 log ratio (lm)") }
  U }
`%||%` <- function(a, b) if (is.null(a)) b else a
# per-cohort estimates of a unit list, in the layout of cohort_estimates.tsv (W = OLS t, as 121)
a48_est <- function(U) do.call(rbind, lapply(U, function(u) { r <- a46_r(u$E, u$D$rx); df <- u$D$df; ci <- a46_ci(r, df); t <- r2t(r, df)
  data.frame(arm=u$arm, id=u$ids, cohort=u$cohort, table=u$prefix, n=u$D$n, df=df, t=t, W=t, r=r, lo=ci[, "lo"], hi=ci[, "hi"], p_t=2 * pt(-abs(t), df), p_fl=NA, n_ref=u$n_ref %||% NA) }))

# ---------------------------------------------------------------- 47.2 classification, copied from 121_a47_labels.R (pool_with, classify_with, run_rfi, run_met; parametric p)
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
# labels for both arms from per-cohort estimates (columns id, cohort, W, df, n, p_fl), sets and testability as 121
a48_labels <- function(root, est, method="mkh", kp="t") {
  tk0 <- rd(file.path(root, "plan", "target_keys.tsv")); cols <- c("id", "cohort", "W", "df", "n", "p_fl")
  tk <- tk0[!tk0$set %in% c("Methane", "MethaneSecondary"), ]
  tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available", ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
  tm <- tk0[tk0$set %in% c("Methane", "MethaneSecondary") | (tk0$set == "Nelore" & tk0$phenotype == "RME_GC"), ]
  ER <- est[est$arm == "RFI", cols]; EM <- est[est$arm == "Methane", cols]
  R <- do.call(rbind, lapply(unique(tk$set), function(st) { d <- tk[tk$set == st, ]
    tg <- data.frame(id=d$id, dir=pR(d$direction), testable=d$stratum %in% c("replication", "cross-site")); cbind(arm="RFI", set=st, target=d$target, published_dir=tg$dir, classify_with(tg, ER[ER$id %in% d$id, ], method, kp)) }))
  M <- do.call(rbind, lapply(unique(tm$set), function(st) { d <- tm[tm$set == st, ]
    tg <- data.frame(id=d$id, dir=pM(d$direction), testable=d$status %in% c("mapped", "special")); cbind(arm="Methane", set=st, target=d$target, published_dir=tg$dir, classify_with(tg, EM[EM$id %in% d$id, ], method, kp)) }))
  rbind(R, M) }
# pooled mKH interval and benchmark exclusions per target from per-cohort estimates (47.1 rule; k = 1 keeps the cohort interval)
a48_pool_targets <- function(est, dirs) do.call(rbind, lapply(split(est, paste(est$arm, est$id)), function(e) { e <- e[is.finite(e$r), ]; k <- nrow(e); id <- e$id[1]; arm <- e$arm[1]; dr <- dirs[[id]] %||% NA
  if (k == 0) return(NULL)
  if (k == 1) { lo <- e$lo; hi <- e$hi; r <- e$r; tau2 <- NA } else { P <- a47_pool(atanh(e$r), 1 / (e$df - 1)); r <- tanh(P[["mu"]]); lo <- tanh(P[["lo_mkh"]]); hi <- tanh(P[["hi_mkh"]]); tau2 <- P[["tau2"]] }
  data.frame(arm=arm, id=id, k=k, cohorts=paste(e$cohort, collapse=";"), r=r, lo=lo, hi=hi, tau2=tau2, published_dir=dr, excl_0=lo > 0 | hi < 0,
    excl_15=a47_excl(lo, hi, dr, 0.15), excl_30=a47_excl(lo, hi, dr, 0.30), bound=a46_bound(lo, hi, dr)) }))
a48_dirs <- function(root) { tk0 <- rd(file.path(root, "plan", "target_keys.tsv")); m <- tk0$set %in% c("Methane", "MethaneSecondary") | (tk0$set == "Nelore" & tk0$phenotype == "RME_GC")
  d <- ifelse(m, pM(tk0$direction), pR(tk0$direction)); setNames(as.list(d), paste0(tk0$id)) }
