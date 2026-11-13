# Amendment 46: late analysis, run after all primary results were seen. Shared functions for 112 (46.1, 46.3), 113 (46.2) and 114 (46.4).
# The design objects (count table A with summed target features, target map feat, metadata md, model terms, design groups) are dumped by the
# A46=1 hooks in 21, 65 and 81 at the point where those scripts call lib_clr.R, so the CLR outcome and model are those of the second estimator
# (Amendment 42). Dumps go to work/a46_inputs/ (never committed: they hold phenotype values). Outputs only in results/late_a46/; nothing is overwritten.
suppressMessages(library(metafor))
A46_LABEL <- "late analysis, run after all primary results were seen (Amendment 46)"
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
a46_gate <- function(root) {
  pl <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")); stopifnot("commit Amendment 46 first" = any(grepl("^\\*\\*Amendment 46 \\(", pl)))
  stopifnot("Amendment 46 never runs on shuffled data" = !nzchar(Sys.getenv("SHUFFLE_SEED"))) }
a46_dir <- function(root, ...) { d <- file.path(root, "results", "late_a46", ...); dir.create(d, recursive=TRUE, showWarnings=FALSE); d }
# refuse to start if any output of a script already exists (outputs are written together at the end of each script)
a46_fresh <- function(files) { ex <- files[file.exists(files)]; if (length(ex)) stop("output exists, Amendment 46 never overwrites: ", paste(ex, collapse=", ")) }
a46_write <- function(x, f) { stopifnot("output exists" = !file.exists(f)); tmp <- paste0(f, ".tmp", Sys.getpid()); write.table(x, tmp, sep="\t", quote=FALSE, row.names=FALSE); invisible(file.rename(tmp, f)) }

# ---------------------------------------------------------------- dump (called from 21, 65, 81)
a46_dump <- function(A, feat, tk, md, terms, df, cohort, rank, prefix, root, perm_by, m1=NULL, difford=FALSE) {
  a46_gate(root); A <- as.matrix(A); stopifnot(all(colnames(A) %in% rownames(md))); md <- md[colnames(A), , drop=FALSE]
  if (!all(perm_by %in% names(md))) stop("design-group column(s) missing from the model data: ", paste(setdiff(perm_by, names(md)), collapse=", "))
  d <- file.path(root, "work", "a46_inputs"); dir.create(d, recursive=TRUE, showWarnings=FALSE); f <- file.path(d, paste0(prefix, ".rds"))
  if (file.exists(f)) { cat("exists, not overwritten:", f, "\n"); return(invisible()) }
  saveRDS(list(A=A, feat=feat, tk=tk, md=md, terms=terms, df=df, cohort=cohort, rank=rank, prefix=prefix, perm_by=perm_by, m1=m1, difford=difford, made=format(Sys.time())), f)
  cat(A46_LABEL, "| dump", prefix, "| animals", ncol(A), "| features", nrow(A), "| df", df, "| design groups", paste(perm_by, collapse=" x "),
    if (!is.null(m1)) paste("| M1 cows", nrow(m1$md)) else "", "\n") }
# M1 clade log ratio exactly as in 81_estimate_methane.R (cows with >= 100 clade-assigned reads)
a46_m1 <- function(root, coh, md) {
  m1 <- rd(file.path(root, "profiles_16s", coh, paste0(coh, "_m1_counts.tsv"))); m1$cow <- as.character(m1$cow); m1 <- m1[match(rownames(md), m1$cow), ]
  ok <- !is.na(m1$cow) & (m1$gottschalkii + m1$ruminantium) >= 100
  list(md=md[ok, , drop=FALSE], y=log((m1$gottschalkii[ok] + 0.5) / (m1$ruminantium[ok] + 0.5))) }
# Freedman-Lane trait for one ANCOM-BC2 permutation (46.3; used by 81 with A46_FL_SEED)
a46_fl_md <- function(md, terms, perm_by, seed) {
  D <- a46_design(md, terms); g <- interaction(md[, perm_by, drop=FALSE], drop=TRUE); set.seed(seed); P <- a46_perm_idx(g, 1)[, 1]
  v <- terms[length(terms)]; stopifnot(is.numeric(md[[v]])); md[[v]] <- D$fx + D$rx[P]; md }

# ---------------------------------------------------------------- estimator (46.1): CLR as lib_clr.R, OLS partial correlation via Frisch-Waugh-Lovell
a46_clr <- function(A) { B <- A[rowMeans(A > 0) >= 0.10, , drop=FALSE]; L <- log(B + 0.5); sweep(L, 2, colMeans(L[!grepl("^TGT_", rownames(L)), , drop=FALSE])) }
a46_design <- function(md, terms) {
  md <- droplevels(md); mm <- model.matrix(reformulate(terms), md); v <- terms[length(terms)]; j <- grep(paste0("^", v), colnames(mm)); stopifnot(length(j) == 1, qr(mm)$rank == ncol(mm))
  Z <- mm[, -j, drop=FALSE]; qz <- qr(Z); x <- mm[, j]; rx <- qr.resid(qz, x)
  list(mm=mm, qz=qz, x=x, rx=rx, fx=x - rx, df=nrow(mm) - ncol(mm), vcol=colnames(mm)[j], n=nrow(mm)) }
# rows of Y (features x animals) -> covariate-residualised rows
a46_resid <- function(Y, D) { Y <- matrix(Y, ncol=D$n, dimnames=list(rownames(Y), NULL)); E <- t(qr.resid(D$qz, t(Y))); rownames(E) <- rownames(Y); E }
# partial r of each row with the trait given the covariates (= OLS t of the trait in the full model, r = t / sqrt(t^2 + df))
a46_r <- function(E, rx) as.vector(E %*% rx) / (sqrt(rowSums(E^2)) * sqrt(sum(rx^2)))
r2t <- function(r, df) r * sqrt(df / (1 - r^2))
a46_ci <- function(r, df) { z <- atanh(r); s <- 1 / sqrt(df - 1); cbind(lo=tanh(z - 1.96 * s), hi=tanh(z + 1.96 * s)) }

# ---------------------------------------------------------------- permutations within design groups
a46_perm_idx <- function(g, B) { g <- as.integer(factor(g)); n <- length(g); P <- matrix(seq_len(n), n, B)
  for (lv in sort(unique(g))) { i <- which(g == lv); m <- length(i); if (m > 1) P[i, ] <- i[apply(matrix(runif(m * B), m), 2, order)] }; P }
# Freedman-Lane: permuted trait = fitted + permuted residual; its covariate residual is M_Z P rx, so the refitted partial r is cor(E, M_Z P rx).
# Exceedance counts accumulated per feature (no null matrix stored). seed=NULL continues the current RNG stream.
a46_fl <- function(E, rx, qz, g, B, seed=NULL, batch=2000) {
  if (!is.null(seed)) set.seed(seed); nE <- sqrt(rowSums(E^2)); r0 <- as.vector(E %*% rx) / (nE * sqrt(sum(rx^2))); bad0 <- !is.finite(r0)
  b <- numeric(nrow(E)); nf <- numeric(nrow(E)); done <- 0
  while (done < B) { k <- min(batch, B - done); P <- a46_perm_idx(g, k); U <- qr.resid(qz, matrix(rx[P], nrow(P))); U <- sweep(U, 2, sqrt(colSums(U^2)), "/")
    R <- (E %*% U) / nE; ok <- is.finite(R); ex <- abs(R) >= abs(r0) - 1e-12; ex[!ok | is.na(ex)] <- FALSE
    b <- b + rowSums(ex); nf <- nf + rowSums(!ok); done <- done + k }
  Bv <- B - nf; p <- (b + 1) / (Bv + 1); p[bad0] <- NA
  data.frame(r=r0, b=b, B=Bv, failed=nf + bad0 * B, p_fl=p, mcse=sqrt(p * (1 - p) / Bv)) }

# ---------------------------------------------------------------- pooling (46.1): REML with HKSJ, fixed effect and REML-Wald beside it
a46_pool <- function(r, df) {
  k <- sum(is.finite(r)); ok <- is.finite(r); r <- r[ok]; df <- df[ok]; z <- atanh(r); v <- 1 / (df - 1)
  na <- data.frame(k=k, r=NA, lo=NA, hi=NA, p=NA, interval=NA, r_fe=NA, lo_fe=NA, hi_fe=NA, lo_wald=NA, hi_wald=NA, tau2=NA, I2=NA)
  if (k == 0) return(na)
  if (k == 1) { ci <- a46_ci(r, df); p <- 2 * pt(-abs(r2t(r, df)), df); na[c("r", "lo", "hi", "p", "r_fe", "lo_fe", "hi_fe", "lo_wald", "hi_wald")] <- c(r, ci, p, r, ci, ci); na$interval <- "cohort (k = 1)"; return(na) }
  fit <- function(test, method="REML") rma(z, v, method=method, test=test, control=list(stepadj=0.5, maxiter=1000))
  hk <- tryCatch(fit("knha"), error=function(e) fit("knha", "DL")); wa <- tryCatch(fit("z"), error=function(e) fit("z", "DL")); fe <- rma(z, v, method="FE")
  data.frame(k=k, r=tanh(hk$beta[1]), lo=tanh(hk$ci.lb), hi=tanh(hk$ci.ub), p=hk$pval, interval=paste0("HKSJ (", hk$method, ", t on ", k - 1, " df)"),
    r_fe=tanh(fe$beta[1]), lo_fe=tanh(fe$ci.lb), hi_fe=tanh(fe$ci.ub), lo_wald=tanh(wa$ci.lb), hi_wald=tanh(wa$ci.ub), tau2=hk$tau2, I2=hk$I2) }
# largest effect in the published direction inside the interval (both bounds where no direction was published), as 93_power_equivalence.R
a46_bound <- function(lo, hi, dir) ifelse(is.na(dir), pmax(abs(lo), abs(hi)), ifelse(dir > 0, hi, -lo))
pR <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))   # as 96_clr_compare.R
pM <- function(x) ifelse(grepl("lower with higher", x), -1, ifelse(grepl("higher with higher", x), 1, NA))

# ---------------------------------------------------------------- the primary target tables (as 96_clr_compare.R) as estimation units
A46_RFI <- c(paste0(rep(c("nelore", "charolais_mrna"), each=4), "_", c("species", "genus", "family", "phylum")), "alberta_genus_all", "alberta_family_all")
A46_MET <- c("difford_bact", "wallace_prokaryotes", "difford_arch")
a46_load <- function(root, prefix) { f <- file.path(root, "work", "a46_inputs", paste0(prefix, ".rds")); if (!file.exists(f)) stop("dump missing (run the A46=1 step first): ", f); readRDS(f) }
a46_units <- function(root) {
  tk0 <- rd(file.path(root, "plan", "target_keys.tsv")); U <- list()
  add <- function(prefix, arm, cohort, ids, Y, md, terms, perm_by, kind) {
    D <- a46_design(md, terms); g <- interaction(md[, perm_by, drop=FALSE], drop=TRUE); rownames(Y) <- ids
    U[[length(U) + 1]] <<- list(prefix=prefix, arm=arm, cohort=cohort, ids=ids, Y=Y, E=a46_resid(Y, D), D=D, g=g, perm_by=perm_by, kind=kind) }
  for (p in c(A46_RFI, A46_MET)) { d <- a46_load(root, p); C <- a46_clr(d$A); i <- match(d$feat, rownames(C)); ok <- !is.na(i)
    stopifnot("df differs from the dumped design" = a46_design(d$md, d$terms)$df == d$df)
    if (p %in% A46_RFI) { coh <- d$cohort; arm <- "RFI" } else {
      coh <- if (startsWith(p, "difford")) "difford" else "wallace"; arm <- "Methane"
      drop <- d$tk$id %in% c("M1", if (d$difford) "M3", if (coh == "wallace") "M2")   # M1: linear model below; M3 untestable in Difford; Wallace M2 re-analysis (Amendment 32)
      if (p == "difford_arch") drop <- rep(TRUE, length(ok))                         # difford_arch supplies M1 only (as 83 and 96)
      ok <- ok & !drop }
    if (any(ok)) add(p, arm, coh, d$tk$id[ok], C[i[ok], , drop=FALSE], d$md, d$terms, d$perm_by, "CLR")
    if (!is.null(d$m1)) add(paste0(p, "_M1"), "Methane", coh, "M1", matrix(d$m1$y, 1), d$m1$md, d$terms, d$perm_by, "M1 log ratio (lm)") }
  pairs <- do.call(rbind, lapply(U, function(u) data.frame(id=u$ids, cohort=u$cohort, prefix=u$prefix))); stopifnot("a target appears twice for one cohort" = !anyDuplicated(pairs[, 1:2]))
  attr(U, "tk") <- tk0; U }
# arm is recycled to the length of id (fix after the first 113 run: with a single arm value, ifelse() returned the first target's direction for every target of a unit)
a46_dir_of <- function(tk0, id, arm) { j <- match(id, tk0$id); arm <- rep_len(arm, length(id)); ifelse(arm == "RFI", pR(tk0$direction[j]), pM(tk0$direction[j])) }
