#!/usr/bin/env Rscript
# Amendment 47.3: late analysis, run after all primary, late and Amendment 46 results were seen. Null-scan calibration at an adequate replicate count.
# For each of the six 46.4 tables: 200 outer null replicates. Each is a Freedman-Lane-permuted trait (trait residuals permuted within design groups; covariate residual
# M_Z P r_x), scanned exactly as in 46.4 (CLR OLS, same genera, 99,999 or 9,999 inner Freedman-Lane permutations, BH within table).
# Seeds 4730000 + 1000 x table index + replicate; the inner permutations continue that stream. Replicate 1 of every table is timed first and the projected wall time for 200
# is written to projection.txt before the rest runs; a table projected over 48 h on its cores runs the largest count that fits, written there first.
# Addendum item 5: a table whose Clopper-Pearson 95% CI from 200 replicates contains both 0.05 and 0.10 is extended to 1,000 (replicates 201 to 1,000, same seeding)
# if the projected time is within 48 h; that projection is appended to projection.txt before the extension runs.
# Usage: 122_a47_null_scan.R. Env A47_CORES (default 4); A47_NULL_REPS, A47_SCAN_B (test runs only). Writes results/late_a47/null_scan/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); suppressMessages(library(parallel)); a47_gate(root); od <- a47_dir(root, "null_scan")
CORES <- as.integer(Sys.getenv("A47_CORES", "4")); NREP <- as.integer(Sys.getenv("A47_NULL_REPS", "200")); if (NREP != 200) cat("NOTE: A47_NULL_REPS =", NREP, "(test run; the plan fixes 200)\n")
TABS <- c(nelore_genus=99999, charolais_mrna_genus=99999, alberta_genus_all=9999, difford_bact=9999, wallace_prokaryotes=9999, wallace_prokaryotes_depth=9999)
tb <- Sys.getenv("A47_SCAN_B"); if (nzchar(tb)) { TABS[] <- as.integer(tb); cat("NOTE: A47_SCAN_B =", tb, "(test run)\n") }
outs <- file.path(od, c("projection.txt", "replicates.tsv", "summary.tsv")); a47_fresh(outs)
cat(A47_LABEL, "| 47.3 null-scan calibration | tables", length(TABS), "| outer replicates", NREP, "| cores", CORES, "\n")
S <- lapply(names(TABS), function(p) { d <- a46_load(root, p); C <- a46_clr(d$A); C <- C[!grepl("^TGT_", rownames(C)), , drop=FALSE]; D <- a46_design(d$md, d$terms)
  list(prefix=p, E=a46_resid(C, D), D=D, g=interaction(d$md[, d$perm_by, drop=FALSE], drop=TRUE), perm_by=d$perm_by) }); names(S) <- names(TABS)
# genera must be those of the 46.4 scan
for (t in seq_along(TABS)) { o <- rd(file.path(root, "results", "late_a46", "scan", paste0(names(TABS)[t], "_scan.tsv"))); stopifnot("genera differ from the 46.4 scan" = setequal(o$genus, rownames(S[[t]]$E))) }
run <- function(t, i) { s <- S[[t]]; t0 <- proc.time()[[3]]; seed <- 4730000 + 1000 * t + i; set.seed(seed)
  P <- a46_perm_idx(s$g, 1)[, 1]; rx <- qr.resid(s$D$qz, s$D$rx[P]); f <- a46_fl(s$E, rx, s$D$qz, s$g, TABS[[t]]); q <- p.adjust(f$p_fl, "BH")
  data.frame(table=names(TABS)[t], replicate=i, seed=seed, genera=nrow(s$E), inner_permutations=TABS[[t]], n_q05=sum(q < 0.05, na.rm=TRUE), n_q10=sum(q < 0.10, na.rm=TRUE),
    frac_p05=mean(f$p_fl < 0.05, na.rm=TRUE), min_p=min(f$p_fl, na.rm=TRUE), failed=sum(f$failed), seconds=proc.time()[[3]] - t0) }
# ---- replicate 1 of every table, timed; projection and the 48 h rule before anything else runs
R1 <- do.call(rbind, mclapply(seq_along(TABS), run, i=1, mc.cores=min(CORES, length(TABS))))
R1$projected_h_200_own_cores <- R1$seconds * 200 / CORES / 3600; R1$reps <- ifelse(R1$projected_h_200_own_cores > 48, pmax(1, floor(48 * 3600 * CORES / R1$seconds)), NREP)
tot <- sum(R1$seconds * (R1$reps - 1)) / CORES / 3600
proj <- c(A47_LABEL, sprintf("replicate 1 timed per table (s): %s", paste(R1$table, round(R1$seconds, 1), collapse="; ")),
  sprintf("projected wall time for %d replicates per table on %d cores: %s", NREP, CORES, paste(R1$table, sprintf("%.2f h", R1$projected_h_200_own_cores), collapse="; ")),
  sprintf("replicates to run (48 h rule): %s", paste(R1$table, R1$reps, collapse="; ")), sprintf("projected wall time for the rest, all tables together: %.2f h", tot))
writeLines(proj, outs[1]); cat(proj, sep="\n")
jobs <- do.call(rbind, lapply(seq_along(TABS), function(t) if (R1$reps[t] > 1) data.frame(t=t, i=2:R1$reps[t]))); jobs <- jobs[order(-TABS[jobs$t], jobs$t, jobs$i), ]
t0 <- proc.time()[[3]]; RR <- mclapply(seq_len(nrow(jobs)), function(j) run(jobs$t[j], jobs$i[j]), mc.cores=CORES, mc.preschedule=FALSE)
bad <- sapply(RR, function(x) inherits(x, "try-error") || !is.data.frame(x)); if (any(bad)) stop("null-scan replicates failed: ", sum(bad)); cat("null scans done,", round(proc.time()[[3]] - t0), "s\n")
R <- rbind(R1[, names(RR[[1]])], do.call(rbind, RR)); R <- R[order(R$table, R$replicate), ]
# ---- addendum item 5: extension to 1,000 replicates where the CI from 200 cannot separate 0.05 from 0.10
NEXT <- as.integer(Sys.getenv("A47_NULL_EXT", "1000")); amb <- sapply(names(TABS), function(tn) { x <- R[R$table == tn, ]; ci <- a47_cp(sum(x$n_q05 > 0), nrow(x))
  nrow(x) == NREP && ci[1] <= 0.05 && ci[2] >= 0.10 })
sec <- setNames(R1$seconds, R1$table); hx <- sum(sec[names(TABS)[amb]] * (NEXT - NREP)) / CORES / 3600; ext <- names(TABS)[amb]
if (length(ext) && hx > 48) { ext <- character(0) }
add <- c(sprintf("extension rule (CI from %d contains both 0.05 and 0.10): %s", NREP, if (any(amb)) paste(names(TABS)[amb], collapse=", ") else "no table"),
  if (any(amb)) sprintf("projected wall time to extend to %d: %.2f h -> %s", NEXT, hx, if (length(ext)) "run" else "not run (over 48 h)"))
cat(add, sep="\n"); write(add, outs[1], append=TRUE)
if (length(ext)) { jx <- expand.grid(i=(NREP + 1):NEXT, t=match(ext, names(TABS))); t0 <- proc.time()[[3]]
  RX <- mclapply(seq_len(nrow(jx)), function(j) run(jx$t[j], jx$i[j]), mc.cores=CORES, mc.preschedule=FALSE)
  if (any(sapply(RX, function(x) !is.data.frame(x)))) stop("extension replicates failed"); R <- rbind(R, do.call(rbind, RX)); R <- R[order(R$table, R$replicate), ]
  cat("extension done,", round(proc.time()[[3]] - t0), "s\n") }
SU <- do.call(rbind, lapply(names(TABS), function(tn) { x <- R[R$table == tn, ]; a <- sum(x$n_q05 > 0); ci <- a47_cp(a, nrow(x)); tt <- table(x$n_q05)
  data.frame(table=tn, genera=x$genera[1], inner_permutations=x$inner_permutations[1], replicates=nrow(x), any_rejection=a, p_any_rejection=a / nrow(x), cp_lo=ci[1], cp_hi=ci[2],
    discovery_counts=paste(names(tt), tt, sep=":", collapse="; "), mean_discoveries=mean(x$n_q05), max_discoveries=max(x$n_q05), mean_frac_p05=mean(x$frac_p05), failed=sum(x$failed),
    note="P(any BH rejection at q < 0.05) under the complete null equals the FDR; Clopper-Pearson 95% CI") }))
a47_write(R, outs[2]); a47_write(SU, outs[3])
cat("\n== 47.3 null-scan calibration: probability of any BH rejection (q < 0.05) under the complete null\n"); print(SU[, setdiff(names(SU), "note")], digits=3, row.names=FALSE)
