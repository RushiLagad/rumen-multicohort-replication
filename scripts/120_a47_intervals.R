#!/usr/bin/env Rscript
# Amendment 47.1 (with the Amendment 47 addendum): late analysis, run after all primary, late and Amendment 46 results were seen.
# Validation of the designated pooled-interval procedure, REML with modified Knapp-Hartung (adhoc); HKSJ, REML-Wald and fixed effect beside it as comparators.
# One REML fit (DL if REML fails, recorded) gives HKSJ, modified KH and REML-Wald in closed form (checked against metafor on the real data); fixed effect uses its own weights.
# Check 1 (pooling alone): theta_i ~ N(mu, tau^2), z_i ~ N(theta_i, 1/(df_i - 1)), mu = atanh(rho), for every configuration; seeds 4710000 + b.
# Check 2 (end to end): the 46.2 simulation on the real traits and covariates with the corrected generator (CLR residual centred within the permutation groups
# (addendum 2), permuted, residualised on the covariates and rescaled to the length of t_r), verified first over 20,000 permutations per table (seeds 4729000 + unit); seeds 4720000 + b, shared by all tables.
# Grid rho in {0, +-0.10, +-0.15, +-0.20, +-0.30}, tau in {0, 0.05, 0.10, 0.20, 0.30} (45 cells); configurations: primary pooling, its leave-one-out pools, four-cohort RFI
# pooling (47.4) and its leave-one-out pools. B = 1,000; a cell where modified KH lies within 2 MCSE of 0.93 is extended to 5,000 and decided there.
# Validated: modified KH passes (c + 2 MCSE >= 0.93) every cell with tau <= 0.10 of the primary configurations (Check 1) and the 18 primary targets (Check 2).
# Per pooled target: tau_req = smallest grid value >= max(0.10, tau-hat), capped at 0.30; its statement stands only if its configuration and the target pass up to tau_req.
# Usage: 120_a47_intervals.R. Env A47_CORES (default 4); A47_REPS, A47_REPS_MAX, A47_GEN_PERMS (test runs only). Writes results/late_a47/intervals/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); suppressMessages(library(parallel)); a47_gate(root); od <- a47_dir(root, "intervals")
CORES <- as.integer(Sys.getenv("A47_CORES", "4")); B0 <- as.integer(Sys.getenv("A47_REPS", "1000")); BX <- as.integer(Sys.getenv("A47_REPS_MAX", "5000")); NGEN <- as.integer(Sys.getenv("A47_GEN_PERMS", "20000"))
TEST <- nzchar(Sys.getenv("A47_TEST")); if (B0 != 1000 || BX != 5000 || NGEN != 20000) cat("NOTE: test run (A47_REPS", B0, "A47_REPS_MAX", BX, "A47_GEN_PERMS", NGEN, "); the plan fixes 1,000, 5,000 and 20,000\n")
RHO <- c(0, 0.10, -0.10, 0.15, -0.15, 0.20, -0.20, 0.30, -0.30); TAU <- c(0, 0.05, 0.10, 0.20, 0.30); THR <- 0.93; DELTA <- c(0.15, 0.30); CAND <- "mkh"
CELLS <- expand.grid(rho=RHO, tau=TAU); CELLS$mu <- atanh(CELLS$rho); CELLS$cell <- seq_len(nrow(CELLS)); NC <- nrow(CELLS)
outs <- file.path(od, c("generator_check.tsv", "check1_cells.tsv", "check2_cells.tsv", "check2_cohort.tsv", "validation.tsv", "chosen_method.txt", "pooled_intervals.tsv",
  "exclusion_counts.tsv", "closed_form_check.tsv", "wallace_rfi_estimates.tsv", "projection.txt", "configuration_validation.tsv")); a47_fresh(outs)
cat(A47_LABEL, "| 47.1 | designated procedure: modified Knapp-Hartung | cores", CORES, "| B", B0, "->", BX, "| cells", NC, "\n")

# ---------------------------------------------------------------- units: the 14 Amendment 46 estimation units plus Wallace RFI NRC (47.4)
U <- a46_units(root); tk0 <- attr(U, "tk"); U[[length(U) + 1]] <- a47_wallace_rfi_unit(root)
for (k in seq_along(U)) { u <- U[[k]]; U[[k]]$r <- a46_r(u$E, u$D$rx); U[[k]]$df <- u$D$df; U[[k]]$key <- paste(u$arm, u$ids, u$cohort); U[[k]]$dir <- a46_dir_of(tk0, u$ids, u$arm) }
ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv"))
A <- do.call(rbind, lapply(U, function(u) data.frame(arm=u$arm, id=u$ids, cohort=u$cohort, table=u$prefix, n=u$D$n, df=u$df, r=u$r, dir=u$dir, key=u$key)))
j <- match(paste(A$arm, A$id, A$cohort), paste(ce$arm, ce$id, ce$cohort)); a46 <- !is.na(j)
stopifnot("Amendment 46 estimates not reproduced" = sum(a46) == nrow(ce), max(abs(A$r[a46] - ce$r[j[a46]])) < 1e-10, all(A$df[a46] == ce$df[j[a46]]))
W <- A[A$table == "wallace_prokaryotes_rfi_nrc", ]; stopifnot("Wallace RFI NRC df must equal Amendment 43 (993)" = all(W$df == 993) || TEST)
W$t <- r2t(W$r, W$df); W$p_t <- 2 * pt(-abs(W$t), W$df); ci <- a46_ci(W$r, W$df); W$lo <- ci[, "lo"]; W$hi <- ci[, "hi"]
cat("Wallace RFI NRC (CLR, 46.1 estimand): n", W$n[1], "df", W$df[1], "targets", nrow(W), "\n")
# pooling versions: primary = the Amendment 46 estimates; rfi4 = RFI targets that Wallace adds to, with Wallace
grp <- function(ver) { x <- if (ver == "primary") A[A$table != "wallace_prokaryotes_rfi_nrc", ] else A[A$arm == "RFI" & A$id %in% W$id, ]
  s <- split(seq_len(nrow(x)), paste(x$arm, x$id)); lapply(s, function(i) x[i, ]) }
G <- list(primary=grp("primary"), rfi4=grp("rfi4"))
TG <- do.call(rbind, lapply(names(G), function(v) do.call(rbind, lapply(G[[v]], function(x) data.frame(version=v, arm=x$arm[1], id=x$id[1], k=nrow(x), dir=x$dir[1],
  cohorts=paste(x$cohort, collapse=";"), dfs=paste(x$df, collapse=";"), cfg=paste(sort(x$df), collapse=";"), keys=paste(x$key, collapse="|"))))))
# pools to evaluate: every main pool with k >= 2 and every leave-one-out pool with k >= 2
TG2 <- do.call(rbind, lapply(seq_len(nrow(TG)), function(i) { tg <- TG[i, ]; if (tg$k < 2) return(NULL); x <- G[[tg$version]][[paste(tg$arm, tg$id)]]
  main <- data.frame(role=tg$version, version=tg$version, arm=tg$arm, id=tg$id, k=tg$k, dir=tg$dir, left_out=NA, cfg=tg$cfg, keys=tg$keys)
  loo <- if (tg$k >= 3) do.call(rbind, lapply(seq_len(nrow(x)), function(l) { y <- x[-l, ]
    data.frame(role=paste0(tg$version, "_loo"), version=tg$version, arm=tg$arm, id=tg$id, k=nrow(y), dir=tg$dir, left_out=x$cohort[l], cfg=paste(sort(y$df), collapse=";"), keys=paste(y$key, collapse="|")) }))
  rbind(main, loo) })); rownames(TG2) <- NULL
cat("pools with k >= 2:", paste(names(table(TG2$role)), table(TG2$role), collapse=", "), "\n")
stopifnot("plan: 18 primary targets with k >= 2" = sum(TG2$role == "primary") == 18 || TEST)
CF <- data.frame(cfg=unique(TG2$cfg)); CF$k <- sapply(strsplit(CF$cfg, ";"), length)
CF$roles <- sapply(CF$cfg, function(c) paste(sort(unique(TG2$role[TG2$cfg == c])), collapse=",")); CF$primary <- sapply(CF$cfg, function(c) any(TG2$cfg == c & TG2$role == "primary"))
CF$targets <- sapply(CF$cfg, function(c) paste(unique(paste(TG2$arm, TG2$id)[TG2$cfg == c]), collapse=", "))
cat("distinct cohort configurations:", nrow(CF), "(primary", sum(CF$primary), ")\n")

# ---------------------------------------------------------------- timing probe and projection
t0 <- proc.time()[[3]]; for (i in 1:300) a47_pool(rnorm(3, 0, 0.05), c(0.02, 0.03, 0.002)); tp <- 1.3 * (proc.time()[[3]] - t0) / 300   # 1.3: task overhead seen in test runs
ntask <- (nrow(CF) + nrow(TG2)) * NC; h <- ntask * B0 * tp / CORES / 3600
proj <- c(A47_LABEL, sprintf("seconds per pooled fit (with overhead): %.4f", tp),
  sprintf("pooling tasks: (%d configurations + %d pools) x %d cells = %d, at B = %d", nrow(CF), nrow(TG2), NC, ntask, B0),
  sprintf("projected wall time at B = %d on %d cores: %.2f h; extension of borderline cells to B = %d adds up to %.2f h if every cell were borderline (typically a small fraction)", B0, CORES, h, BX, h * (BX - B0) / B0))
writeLines(proj, outs[11]); cat(proj, sep="\n")

# ---------------------------------------------------------------- Check 2 inputs and generator verification
# Addendum 2: the residual is centred within the permutation groups before it is permuted (in Charolais the groups, diet x year, are not model terms, so the residual of the
# additive model keeps cell means that a within-cell permutation carries along, correlated with the trait); es_a46 is the uncentred residual of the 46.2 generator, for comparison
for (k in seq_along(U)) { u <- U[[k]]; R <- t(qr.resid(qr(u$D$mm), t(u$Y))); U[[k]]$es_a46 <- R / apply(R, 1, sd)
  Rc <- R - t(apply(R, 1, function(x) ave(x, u$g))); U[[k]]$es <- Rc / apply(Rc, 1, sd); U[[k]]$tr <- u$D$rx / sd(u$D$rx) }
gen_noise <- function(u, P) { et <- a46_resid(u$es[, P, drop=FALSE], u$D); et / apply(et, 1, sd) }   # corrected: residualised on the covariates, length of t_r
t0 <- proc.time()[[3]]
GC <- do.call(rbind, mclapply(seq_along(U), function(k) { u <- U[[k]]; set.seed(4729000 + k); nr <- length(u$ids); T2 <- sum(u$tr^2)
  an <- bn <- ao <- bo <- an2 <- numeric(nr); done <- 0
  while (done < NGEN) { m <- min(500, NGEN - done); PP <- a46_perm_idx(u$g, m)
    for (p in seq_len(m)) { P <- PP[, p]; Mo <- a46_resid(u$es_a46[, P, drop=FALSE], u$D); Mn <- a46_resid(u$es[, P, drop=FALSE], u$D); en <- Mn / apply(Mn, 1, sd)
      a1 <- as.vector(en %*% u$tr); an <- an + a1; an2 <- an2 + a1^2; bn <- bn + rowSums(en^2); ao <- ao + as.vector(Mo %*% u$tr); bo <- bo + rowSums(Mo^2) }
    done <- done + m }
  sa <- sqrt(pmax(0, an2 / NGEN - (an / NGEN)^2)); an <- an / NGEN; bn <- bn / NGEN; ao <- ao / NGEN; bo <- bo / NGEN
  do.call(rbind, lapply(c(0.15, -0.15, 0.30, -0.30), function(rc) { s <- sqrt(1 - rc^2)
    pop <- function(a, b) (rc * T2 + s * a) / sqrt(T2 * (rc^2 * T2 + s^2 * b + 2 * rc * s * a))   # <M_Z x*, t_r> / sqrt(|M_Z x*|^2 |t_r|^2), expectations over permutations
    data.frame(unit=u$prefix, arm=u$arm, id=u$ids, cohort=u$cohort, n=u$D$n, df=u$df, rho_c=rc, rho_pop_corrected=pop(an, bn), rho_pop_a46=pop(ao, bo),
      mcse_corrected=s * sa / T2 / sqrt(NGEN), noise_length_ratio_a46=sqrt(bo / T2), permutations=NGEN, seed=4729000 + k) })) }, mc.cores=CORES))
GC$diff_corrected <- GC$rho_pop_corrected - GC$rho_c; GC$diff_a46 <- GC$rho_pop_a46 - GC$rho_c
cat(sprintf("generator check (%d permutations per table, %.0f s): corrected |diff| max %.4f (largest MCSE %.4f); Amendment 46 generator diff %.4f to %.4f (largest |diff| in %s)\n", NGEN,
  proc.time()[[3]] - t0, max(abs(GC$diff_corrected)), max(GC$mcse_corrected), min(GC$diff_a46), max(GC$diff_a46), GC$cohort[which.max(abs(GC$diff_a46))]))
a47_write(GC, outs[1])
stopifnot("47.1 generator check failed: population partial correlation of the corrected generator differs from rho_c by more than 0.005" = all(abs(GC$diff_corrected) <= 0.005))

# ---------------------------------------------------------------- Check 2 simulation: r for every unit, target and cell, replicates 1..BX (cheap; pooling is the cost)
keys <- sort(unique(A$key)); NK <- length(keys); kidx <- setNames(seq_len(NK), keys)
one2 <- function(b) { set.seed(4720000 + b); zn <- setNames(rnorm(NK), keys); out <- matrix(NA_real_, NK, NC)
  for (u in U) { P <- a46_perm_idx(u$g, 1)[, 1]; et <- gen_noise(u, P); rows <- kidx[u$key]
    for (c in seq_len(NC)) { rc <- tanh(CELLS$mu[c] + CELLS$tau[c] * zn[u$key]); X <- outer(rc, u$tr) + sqrt(1 - rc^2) * et
      out[rows, c] <- a46_r(a46_resid(X, u$D), u$D$rx) } }
  out }
t0 <- proc.time()[[3]]; S2 <- mclapply(seq_len(BX), one2, mc.cores=CORES); S2 <- array(unlist(S2), c(NK, NC, BX)); cat("Check 2 simulation:", BX, "replicates,", round(proc.time()[[3]] - t0), "s\n")
dfk <- A$df[match(keys, A$key)]
CO <- do.call(rbind, lapply(which(CELLS$tau == 0), function(c) { r <- matrix(S2[, c, seq_len(B0)], NK); s <- 1 / sqrt(dfk - 1)   # k = 1 intervals, tau = 0 (cohort rho = rho)
  cov <- abs(atanh(r) - CELLS$mu[c]) <= 1.96 * s; data.frame(key=keys, df=dfk, rho=CELLS$rho[c], reps=B0, coverage=rowMeans(cov), bias_r=rowMeans(r) - CELLS$rho[c],
    rejection=rowMeans(abs(r2t(r, dfk)) > qt(0.975, dfk))) }))
CO <- cbind(A[match(CO$key, A$key), c("arm", "id", "cohort", "table")], CO); CO$mcse <- sqrt(CO$coverage * (1 - CO$coverage) / CO$reps)

# ---------------------------------------------------------------- Check 1 draws: replicates 1..BX, every configuration and cell in a fixed order
cfv <- lapply(strsplit(CF$cfg, ";"), function(x) 1 / (as.numeric(x) - 1)); KM <- max(CF$k); Z1 <- array(NA_real_, c(KM, nrow(CF), NC, BX))
for (b in seq_len(BX)) { set.seed(4710000 + b); for (f in seq_len(nrow(CF))) { v <- cfv[[f]]; k <- length(v)
  for (c in seq_len(NC)) { th <- rnorm(k, CELLS$mu[c], CELLS$tau[c]); Z1[seq_len(k), f, c, b] <- rnorm(k, th, sqrt(v)) } } }

# ---------------------------------------------------------------- pooling tasks
pool_reps <- function(zs, v) t(vapply(seq_len(ncol(zs)), function(b) a47_pool(zs[, b], v), numeric(length(A47_POOL_COLS))))
task_z <- function(tsk, reps) { if (tsk$check == 1) { f <- tsk$f; k <- CF$k[f]; list(z=matrix(Z1[seq_len(k), f, tsk$cell, reps], k), v=cfv[[f]]) } else {
  kk <- kidx[strsplit(TG2$keys[tsk$t], "|", fixed=TRUE)[[1]]]; list(z=atanh(matrix(S2[kk, tsk$cell, reps], length(kk))), v=1 / (dfk[kk] - 1)) } }
run_task <- function(tsk, reps) { x <- task_z(tsk, reps); cbind(b=reps, pool_reps(x$z, x$v)) }
summ <- function(M, rho, dir_list) { mu <- atanh(rho); out <- list()
  for (m in A47_METHODS) { lo <- M[, paste0("lo_", m)]; hi <- M[, paste0("hi_", m)]; est <- if (m == "fe") M[, "mu_fe"] else M[, "mu"]
    cov <- !is.na(lo) & lo <= mu & mu <= hi; cv <- mean(cov); n <- nrow(M)
    row <- data.frame(method=m, reps=n, coverage=cv, mcse=sqrt(cv * (1 - cv) / n), bias_z=mean(est - mu, na.rm=TRUE), bias_r=mean(tanh(est), na.rm=TRUE) - rho,
      median_width_z=median(hi - lo, na.rm=TRUE), failed=sum(is.na(lo)), fallback=sum(M[, "fallback"]))
    for (dn in c("dir", "mag")) { d <- dir_list[[dn]]; for (de in DELTA) { cn <- sprintf("false_excl_%s_%02d", dn, round(100 * de))   # dir: published (Check 1: positive) direction; mag: none
      if (is.null(d)) { row[[cn]] <- NA; next }; reach <- if (is.na(d)) abs(rho) >= de - 1e-12 else d * rho >= de - 1e-12
      row[[cn]] <- if (reach) mean(a47_excl(tanh(lo), tanh(hi), d, de) %in% TRUE) else NA } }
    out[[m]] <- row }
  do.call(rbind, out) }
TASKS <- c(unlist(lapply(seq_len(nrow(CF)), function(f) lapply(seq_len(NC), function(c) list(check=1, f=f, cell=c))), recursive=FALSE),
  unlist(lapply(seq_len(nrow(TG2)), function(t) lapply(seq_len(NC), function(c) list(check=2, t=t, cell=c))), recursive=FALSE))
t0 <- proc.time()[[3]]; RES <- mclapply(TASKS, run_task, reps=seq_len(B0), mc.cores=CORES, mc.preschedule=FALSE)
bad <- sapply(RES, function(x) inherits(x, "try-error") || !is.matrix(x)); if (any(bad)) stop("pooling tasks failed: ", sum(bad)); cat("pooling at B =", B0, ":", round(proc.time()[[3]] - t0), "s\n")
summarise <- function(idx=seq_along(TASKS)) do.call(rbind, lapply(idx, function(i) { tsk <- TASKS[[i]]; rho <- CELLS$rho[tsk$cell]; tau <- CELLS$tau[tsk$cell]
  if (tsk$check == 1) { s <- summ(RES[[i]], rho, list(dir=1, mag=NA)); cbind(task=i, check=1, role=CF$roles[tsk$f], version=NA, arm=NA, id=NA, left_out=NA, cfg=CF$cfg[tsk$f], k=CF$k[tsk$f],
    cell=tsk$cell, rho=rho, tau=tau, rule_cell=tau <= 0.10 & CF$primary[tsk$f], s) } else {
    tg <- TG2[tsk$t, ]; s <- summ(RES[[i]], rho, if (is.na(tg$dir)) list(mag=NA) else list(dir=tg$dir)); cbind(task=i, check=2, role=tg$role, version=tg$version, arm=tg$arm, id=tg$id,
    left_out=tg$left_out, cfg=tg$cfg, k=tg$k, cell=tsk$cell, rho=rho, tau=tau, rule_cell=tau <= 0.10 & tg$role == "primary", s) } }))
SM <- summarise()
# extension: every cell (any configuration, any tau) where modified KH lies within 2 MCSE of 0.93 goes to BX replicates (all four methods together)
xi <- sort(unique(SM$task[SM$method == CAND & abs(SM$coverage - THR) < 2 * SM$mcse]))
cat("cells where modified Knapp-Hartung lies within 2 MCSE of", THR, "at B =", B0, ":", length(xi), "-> extended to B =", BX, "\n")
if (length(xi) && BX > B0) { t0 <- proc.time()[[3]]
  EX <- mclapply(TASKS[xi], run_task, reps=(B0 + 1):BX, mc.cores=CORES, mc.preschedule=FALSE); if (any(sapply(EX, function(x) !is.matrix(x)))) stop("extension tasks failed")
  for (j in seq_along(xi)) RES[[xi[j]]] <- rbind(RES[[xi[j]]], EX[[j]]); cat("extension:", round(proc.time()[[3]] - t0), "s\n"); SM <- summarise() }
SM$pass <- SM$coverage + 2 * SM$mcse >= THR
C1 <- SM[SM$check == 1, setdiff(names(SM), c("check", "version", "arm", "id", "left_out"))]; C2 <- SM[SM$check == 2, setdiff(names(SM), "check")]

# ---------------------------------------------------------------- validation of the designated procedure (comparators reported with the same rule)
VAL <- do.call(rbind, lapply(A47_METHODS, function(m) { a <- SM[SM$rule_cell & SM$method == m, ]; w <- a[a$check == 2, ][which.min(a$coverage[a$check == 2]), ]
  data.frame(method=m, name=A47_NAMES[[m]], designated=m == CAND, cells_check1=sum(a$check == 1), pass_check1=sum(a$pass & a$check == 1), cells_check2=sum(a$check == 2),
    pass_check2=sum(a$pass & a$check == 2), min_coverage_check1=min(a$coverage[a$check == 1]), min_coverage_check2=min(a$coverage[a$check == 2]),
    worst_check2=sprintf("%s %s rho %.2f tau %.2f: %.3f (B %d)", w$arm, w$id, w$rho, w$tau, w$coverage, w$reps), median_width_z_check2=median(a$median_width_z[a$check == 2]),
    passes_rule=all(a$pass)) }))
validated <- VAL$passes_rule[VAL$method == CAND]; chosen <- if (validated) CAND else "none"
cat("\n== 47.1 validation (cells: tau <= 0.10, primary configurations and targets; pass if coverage + 2 MCSE >=", THR, "); comparators under the same rule\n"); print(VAL, row.names=FALSE, digits=3)
# every pool (main, leave-one-out, four-cohort) under modified KH: passes up to tau 0.10, 0.20, 0.30 in Check 1 (its configuration) and Check 2 (the pool itself)
upto <- function(x, t) all(x$pass[x$tau <= t + 1e-9])
CV <- do.call(rbind, lapply(seq_len(nrow(TG2)), function(i) { tg <- TG2[i, ]; a <- C1[C1$method == CAND & C1$cfg == tg$cfg, ]
  b <- C2[C2$method == CAND & C2$role == tg$role & C2$id == tg$id & C2$arm == tg$arm & ((is.na(tg$left_out) & is.na(C2$left_out)) | (!is.na(C2$left_out) & C2$left_out %in% tg$left_out)), ]
  data.frame(tg[, c("role", "version", "arm", "id", "k", "left_out", "cfg")], pass_tau010=upto(a, 0.10) && upto(b, 0.10), pass_tau020=upto(a, 0.20) && upto(b, 0.20),
    pass_tau030=upto(a, 0.30) && upto(b, 0.30), min_coverage_tau010=min(c(a$coverage[a$tau <= 0.10], b$coverage[b$tau <= 0.10]))) }))
txt <- c(A47_LABEL, sprintf("designated procedure: modified Knapp-Hartung (REML, adhoc); validated: %s", if (validated) "yes" else "no (every pooled statement is descriptive; 47.2 uses modified Knapp-Hartung, labelled descriptive)"),
  sprintf("comparators passing the same rule: %s", paste(VAL$method[VAL$passes_rule & !VAL$designated], collapse=", ")), sprintf("B = %d per cell, %d cells extended to %d", B0, length(xi), BX), chosen)

# ---------------------------------------------------------------- the real data: four intervals per target with k >= 2 (both versions), k = 1 the cohort interval
CK <- c(); PI <- do.call(rbind, lapply(seq_len(nrow(TG)), function(i) { tg <- TG[i, ]; x <- G[[tg$version]][[paste(tg$arm, tg$id)]]
  base <- data.frame(version=tg$version, arm=tg$arm, id=tg$id, k=tg$k, cohorts=tg$cohorts, dfs=tg$dfs, cfg=tg$cfg, published_dir=tg$dir)
  if (tg$k == 1) { ci <- a46_ci(x$r, x$df); p <- 2 * pt(-abs(r2t(x$r, x$df)), x$df)
    for (m in A47_METHODS) { base[[paste0("lo_", m)]] <- ci[1, "lo"]; base[[paste0("hi_", m)]] <- ci[1, "hi"]; base[[paste0("p_", m)]] <- p }
    base$r <- x$r; base$r_fe <- x$r; base$tau2 <- NA; base$failed <- 0; base$fallback <- 0; base$interval <- "cohort (k = 1)"; return(base) }
  z <- atanh(x$r); v <- 1 / (x$df - 1); p <- a47_pool(z, v); CK <<- c(CK, a47_pool_check(z, v))
  for (m in A47_METHODS) { base[[paste0("lo_", m)]] <- tanh(p[[paste0("lo_", m)]]); base[[paste0("hi_", m)]] <- tanh(p[[paste0("hi_", m)]]); base[[paste0("p_", m)]] <- p[[paste0("p_", m)]] }
  base$r <- tanh(p[["mu"]]); base$r_fe <- tanh(p[["mu_fe"]]); base$tau2 <- p[["tau2"]]; base$failed <- p[["failed"]]; base$fallback <- p[["fallback"]]; base$interval <- "pooled"; base }))
cfc <- data.frame(targets=length(CK), max_abs_difference=max(CK, na.rm=TRUE), tolerance=1e-8); a47_write(cfc, outs[9])
stopifnot("closed-form intervals differ from metafor" = max(CK, na.rm=TRUE) < 1e-8)
pa <- rd(file.path(root, "results", "late_a46", "effect_size", "pooled.tsv")); pp <- PI[PI$version == "primary", ]; j <- match(paste(pp$arm, pp$id), paste(pa$arm, pa$id))
stopifnot("primary pooling must reproduce the Amendment 46 HKSJ intervals" = !anyNA(j), nrow(pp) == nrow(pa), max(abs(pp$lo_hksj - pa$lo[j]), abs(pp$hi_hksj - pa$hi[j]), abs(pp$r - pa$r[j])) < 1e-8)
j <- match(PI$id, tk0$id); PI$set <- tk0$set[j]; PI$target <- tk0$target[j]; PI$site <- tk0$site[j]
PI$stratum <- ifelse(PI$arm == "Methane" & PI$set == "Nelore", "Nelore RME_GC genus, tested against daily output", ifelse(PI$site == "feces", "fecal source", ifelse(is.na(PI$published_dir), "no published direction", "directed rumen")))
PI$k_class <- ifelse(PI$k >= 2, "pooled (k >= 2)", "single cohort (k = 1)"); PI$claim <- ifelse(is.na(PI$published_dir), "magnitude", "directional")
for (m in A47_METHODS) { for (de in DELTA) PI[[sprintf("excl_%02d_%s", round(100 * de), m)]] <- a47_excl(PI[[paste0("lo_", m)]], PI[[paste0("hi_", m)]], PI$published_dir, de)
  PI[[paste0("excl0_", m)]] <- PI[[paste0("lo_", m)]] > 0 | PI[[paste0("hi_", m)]] < 0 }
# per pooled target: tau-hat, tau_req and whether its statement stands under the designated procedure (addendum item 3)
PI$tau_hat <- sqrt(PI$tau2); PI$tau_req <- ifelse(PI$k >= 2, sapply(pmax(0.10, PI$tau_hat), function(t) if (is.na(t)) NA else c(TAU[TAU >= t - 1e-9], 0.30)[1]), NA)
PI$beyond_range <- PI$k >= 2 & !is.na(PI$tau_hat) & PI$tau_hat > 0.30
PI$pool_passes_to_tau_req <- sapply(seq_len(nrow(PI)), function(i) { if (PI$k[i] < 2) return(NA); x <- CV[CV$role == PI$version[i] & CV$id == PI$id[i] & CV$arm == PI$arm[i] & is.na(CV$left_out), ]
  if (is.na(PI$tau_req[i])) return(FALSE); if (PI$tau_req[i] <= 0.10) x$pass_tau010 else if (PI$tau_req[i] <= 0.20) x$pass_tau020 else x$pass_tau030 })
PI$statement <- ifelse(PI$k < 2, "single cohort (cohort interval)", ifelse(rep(!validated, nrow(PI)), "descriptive (procedure not validated)",   # rep: ifelse with a scalar test returns one value
  ifelse(PI$beyond_range, "descriptive (heterogeneity beyond the simulated range)", ifelse(PI$pool_passes_to_tau_req %in% TRUE, "validated", "descriptive (coverage below 0.93 up to tau_req)"))))
stopifnot(length(PI$statement) == nrow(PI), !any(PI$statement == "validated" & PI$beyond_range %in% TRUE))
PI$designated_method <- CAND; PI$validated_procedure <- validated
# counts (pointwise), by version, arm, single vs pooled, directional vs magnitude, method and benchmark; for modified KH also counting only statements that stand
EC <- do.call(rbind, lapply(split(PI, list(PI$version, PI$arm), drop=TRUE), function(x) do.call(rbind, lapply(c("all", "single cohort (k = 1)", "pooled (k >= 2)"), function(kc)
  do.call(rbind, lapply(c("all", "directional", "magnitude"), function(cl) { y <- x[(kc == "all" | x$k_class == kc) & (cl == "all" | x$claim == cl), ]; ok <- y$statement %in% c("validated", "single cohort (cohort interval)")
    do.call(rbind, lapply(A47_METHODS, function(m) data.frame(version=x$version[1], arm=x$arm[1], k_class=kc, claim=cl, method=m, designated=m == CAND, targets=nrow(y),
      excl_r30=sum(y[[paste0("excl_30_", m)]]), excl_r15=sum(y[[paste0("excl_15_", m)]]), interval_excludes_0=sum(y[[paste0("excl0_", m)]]),
      excl_r30_standing=if (m == CAND) sum(y[[paste0("excl_30_", m)]] & ok) else NA, excl_r15_standing=if (m == CAND) sum(y[[paste0("excl_15_", m)]] & ok) else NA,
      note="pointwise; BH does not adjust exclusion claims"))) }))))))
a47_write(C1, outs[2]); a47_write(C2, outs[3]); a47_write(CO, outs[4]); a47_write(VAL, outs[5]); writeLines(txt, outs[6]); a47_write(PI, outs[7]); a47_write(EC, outs[8])
a47_write(W[, c("arm", "id", "cohort", "table", "n", "df", "r", "t", "p_t", "lo", "hi")], outs[10]); a47_write(CV, outs[12])
cat("\n", txt[1:4], sep="\n")
cat("\n== exclusion counts, primary pooling, modified Knapp-Hartung (pointwise)\n"); print(EC[EC$version == "primary" & EC$designated, c("arm", "k_class", "claim", "targets", "excl_r30", "excl_r15", "excl_r30_standing", "excl_r15_standing", "interval_excludes_0")], row.names=FALSE)
cat("\n== pooled targets (k >= 2), primary: r, the four intervals, tau-hat and the statement\n")
print(PI[PI$version == "primary" & PI$k >= 2, c("arm", "id", "target", "k", "r", "lo_mkh", "hi_mkh", "lo_hksj", "hi_hksj", "lo_wald", "hi_wald", "lo_fe", "hi_fe", "tau_hat", "tau_req", "statement")], digits=3, row.names=FALSE)
cat("\n== Check 2, minimum coverage over rho by target, tau and method (primary)\n"); x <- C2[C2$role == "primary", ]
print(reshape(aggregate(coverage ~ arm + id + tau + method, x, min), idvar=c("arm", "id", "tau"), timevar="method", direction="wide"), digits=3, row.names=FALSE)
cat("\n== pools under modified Knapp-Hartung not passing at tau <= 0.10 (any role):", sum(!CV$pass_tau010), "\n"); if (any(!CV$pass_tau010)) print(CV[!CV$pass_tau010, ], digits=3, row.names=FALSE)
