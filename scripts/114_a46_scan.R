#!/usr/bin/env Rscript
# Amendment 46.4: late analysis, run after all primary results were seen. Genus scan by CLR OLS with Freedman-Lane permutation (exploratory, as Amendment 38).
# Tables and features as the 91 scan: every non-target genus of the genus table after the 10 percent prevalence filter (summed target features excluded).
# Alberta is the all-animal table (Amendment 39). 99,999 permutations for the shotgun and metatranscriptome tables, 9,999 for the 16S tables; BH across genera
# within each table. Seed 470000 + table index for the scan; null check: 20 Freedman-Lane-permuted traits per table, seed 480000 + 20 (table index - 1) + replicate,
# each followed by the full scan from the same RNG stream. Usage: 114_a46_scan.R. Env A46_CORES (default 4); A46_SCAN_B, A46_NULL_REPS (test runs only).
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); suppressMessages(library(parallel)); a46_gate(root); od <- a46_dir(root, "scan"); CORES <- as.integer(Sys.getenv("A46_CORES", "4"))
TABS <- c(nelore_genus=99999, charolais_mrna_genus=99999, alberta_genus_all=9999, difford_bact=9999, wallace_prokaryotes=9999, wallace_prokaryotes_depth=9999)
tb <- Sys.getenv("A46_SCAN_B"); if (nzchar(tb)) { TABS[] <- as.integer(tb); cat("NOTE: A46_SCAN_B =", tb, "(test run)\n") }
NREP <- as.integer(Sys.getenv("A46_NULL_REPS", "20")); if (NREP != 20) cat("NOTE: A46_NULL_REPS =", NREP, "(test run; the plan fixes 20)\n")
outs <- c(file.path(od, paste0(names(TABS), "_scan.tsv")), file.path(od, c("scan_summary.tsv", "null_check.tsv", "wallace_difford_correlation.tsv"))); a46_fresh(outs)
cat(A46_LABEL, "| 46.4 scan | tables", length(TABS), "| null replicates", NREP, "| cores", CORES, "\n")
S <- lapply(names(TABS), function(p) { d <- a46_load(root, p); C <- a46_clr(d$A); C <- C[!grepl("^TGT_", rownames(C)), , drop=FALSE]; D <- a46_design(d$md, d$terms)
  list(prefix=p, E=a46_resid(C, D), D=D, g=interaction(d$md[, d$perm_by, drop=FALSE], drop=TRUE), perm_by=d$perm_by) }); names(S) <- names(TABS)
jobs <- expand.grid(rep=0:NREP, t=seq_along(TABS)); jobs <- jobs[order(-TABS[jobs$t], jobs$t, jobs$rep), ]   # longest first
run <- function(j) { t <- jobs$t[j]; i <- jobs$rep[j]; s <- S[[t]]; B <- TABS[[t]]; t0 <- proc.time()[[3]]
  if (i == 0) { seed <- 470000 + t; set.seed(seed); rx <- s$D$rx } else {
    seed <- 480000 + 20 * (t - 1) + i; set.seed(seed); P <- a46_perm_idx(s$g, 1)[, 1]; rx <- qr.resid(s$D$qz, s$D$rx[P]) }   # Freedman-Lane-permuted trait: its covariate residual
  f <- a46_fl(s$E, rx, s$D$qz, s$g, B); f$q <- p.adjust(f$p_fl, "BH")
  list(t=t, rep=i, seed=seed, f=f, sec=proc.time()[[3]] - t0) }
t0 <- proc.time()[[3]]; R <- mclapply(seq_len(nrow(jobs)), run, mc.cores=CORES, mc.preschedule=FALSE); cat("scan done,", round(proc.time()[[3]] - t0), "s\n")
bad <- sapply(R, function(x) inherits(x, "try-error") || is.null(x$f)); if (any(bad)) stop("scan jobs failed: ", paste(which(bad), collapse=", "))
key <- sapply(R, function(x) paste(x$t, x$rep))
SUM <- list(); NC <- list(); OUT <- list()
for (t in seq_along(TABS)) { s <- S[[t]]; x <- R[[match(paste(t, 0), key)]]; f <- x$f; B <- TABS[[t]]; G <- nrow(s$E); df <- s$D$df
  d <- data.frame(genus=rownames(s$E), n=s$D$n, df=df, t=r2t(f$r, df), r=f$r, p_t=2 * pt(-abs(r2t(f$r, df)), df), b=f$b, B=f$B, failed=f$failed, p_fl=f$p_fl, mcse_p_fl=f$mcse, q_bh=f$q,
    seed=x$seed, perm_within=paste(s$perm_by, collapse=" x "), analysis=A46_LABEL); d <- d[order(d$p_fl, -abs(d$t)), ]; OUT[[t]] <- d
  nul <- R[match(paste(t, seq_len(NREP)), key)]
  nd <- sapply(nul, function(y) sum(y$f$q < 0.05, na.rm=TRUE)); fp <- sapply(nul, function(y) mean(y$f$p_fl < 0.05, na.rm=TRUE))
  NC[[t]] <- data.frame(table=names(TABS)[t], replicate=seq_len(NREP), seed=sapply(nul, `[[`, "seed"), n_q05=nd, frac_p05=fp, failed=sapply(nul, function(y) sum(y$f$failed)), analysis=A46_LABEL)
  SUM[[t]] <- data.frame(table=names(TABS)[t], n=s$D$n, df=df, n_genera=G, permutations=B, smallest_p=1 / (B + 1),
    min_attainable_discoveries=ceiling(G / (0.05 * (B + 1))), n_p05=sum(d$p_fl < 0.05, na.rm=TRUE), n_q05=sum(d$q_bh < 0.05, na.rm=TRUE), n_q10=sum(d$q_bh < 0.10, na.rm=TRUE),
    failed_fits=sum(d$failed), max_mcse=max(d$mcse_p_fl, na.rm=TRUE), null_mean_q05=mean(nd), null_max_q05=max(nd), null_frac_p05=mean(fp),
    top_genus=d$genus[1], top_t=d$t[1], top_p_fl=d$p_fl[1], top_q=d$q_bh[1], minutes=round(x$sec / 60, 1), analysis=A46_LABEL) }
SUM <- do.call(rbind, SUM); NC <- do.call(rbind, NC)
# 46.5.4: Wallace versus Difford genus statistics, descriptive (no p value: genera are not independent)
W <- do.call(rbind, lapply(c("wallace_prokaryotes", "wallace_prokaryotes_depth"), function(w) { a <- OUT[[match("difford_bact", names(TABS))]]; b <- OUT[[match(w, names(TABS))]]; g <- intersect(a$genus, b$genus)
  a <- a[match(g, a$genus), ]; b <- b[match(g, b$genus), ]
  data.frame(table_a="difford_bact", table_b=w, shared_genera=length(g), spearman_t=cor(a$t, b$t, method="spearman"), pearson_r=cor(a$r, b$r),
    same_sign=sum(sign(a$t) == sign(b$t)), both_p05_same_sign=sum(a$p_fl < 0.05 & b$p_fl < 0.05 & sign(a$t) == sign(b$t)), note="descriptive; no p value (genera are not independent)", analysis=A46_LABEL) }))
for (t in seq_along(TABS)) a46_write(OUT[[t]], outs[t]); a46_write(SUM, outs[length(TABS) + 1]); a46_write(NC, outs[length(TABS) + 2]); a46_write(W, outs[length(TABS) + 3])
cat("== 46.4 SCAN (CLR OLS, Freedman-Lane permutation p = (b + 1)/(B + 1), BH within table)\n"); print(SUM[, setdiff(names(SUM), "analysis")], digits=3, row.names=FALSE)
cat("\n== null check: BH discoveries per Freedman-Lane-permuted trait\n"); print(aggregate(n_q05 ~ table, NC, function(v) c(mean=mean(v), max=max(v))), row.names=FALSE)
cat("\n== Wallace vs Difford (descriptive)\n"); print(W[, setdiff(names(W), "analysis")], digits=3, row.names=FALSE)
for (t in seq_along(TABS)) { cat("--", names(TABS)[t], "top 6\n"); print(head(OUT[[t]][, c("genus", "t", "r", "p_fl", "q_bh")], 6), digits=3, row.names=FALSE) }
