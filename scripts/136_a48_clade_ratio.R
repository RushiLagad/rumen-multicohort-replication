#!/usr/bin/env Rscript
# Amendment 48.9: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Clade-ratio robustness.
# Clade counts rebuilt from results/16s/<cohort>/m1_asv_clades.tsv (best mismatch to each clade per ASV) and the DADA2 tables results/16s/<cohort>/dada2.rds (st, per cow):
# an ASV is assigned only if |best_gott - best_rum| >= margin (1 = current, 2, 3), to the closer clade. Margin 1 must reproduce profiles_16s/<cohort>/<cohort>_m1_counts.tsv exactly.
# Grid: margin x pseudo-count (0.1, 0.5 current, 1) x minimum assigned share (70% current, 80%; a cohort below it is untestable in that cell). Cows with >= 100 assigned reads.
# Analyses: 1 Wallace CH4 per kg DMI (A48 dump), 2 Wallace CH4 g/d, 3 Difford CH4 g/d (A46 dumps), 4 Wallace CH4 per kg DMI without farm SE1. Model as 46.1 (terms of the dump);
# lm t, partial r with 95% CI, Freedman-Lane p (9,999 within farm or herd; seed 4890000 + 100 x analysis + cell index). Influence for analysis 1 at the current cell:
# Cook's distance and the refit without D > 4/n; MM regression (robustbase::lmrob, defaults). Reading as 48.9. Usage: 136_a48_clade_ratio.R. Env A48_FL_B (tests only).
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "clade_ratio"); outs <- file.path(od, c("assignment.tsv", "grid.tsv", "influence.tsv", "reading.tsv")); a48_fresh(outs)
if (Sys.getenv("A48_TEST_NO_ROBUST") != "1") stopifnot("robustbase is required (48.9)" = requireNamespace("robustbase", quietly=TRUE)) else cat("NOTE: A48_TEST_NO_ROBUST (test run)\n")
B <- as.integer(Sys.getenv("A48_FL_B", "9999")); if (B != 9999) cat("NOTE: A48_FL_B =", B, "(test run; the plan fixes 9,999)\n")
cat(A48_LABEL, "| 48.9 clade ratio | permutations", B, "\n")
COH <- list(wallace=list(dir="wallace_prokaryotes", prof="wallace_prokaryotes"), difford=list(dir="difford_arch", prof="difford_arch"))
counts_at <- function(co, margin) { a <- rd(file.path(root, "results", "16s", COH[[co]]$dir, "m1_asv_clades.tsv")); st <- readRDS(file.path(root, "results", "16s", COH[[co]]$dir, "dada2.rds"))$st
  j <- match(a$sequence, colnames(st)); stopifnot("ASV sequences missing from the DADA2 table" = !anyNA(j)); d <- a$best_gott - a$best_rum
  cl <- ifelse(abs(d) >= margin & d < 0, "gottschalkii", ifelse(abs(d) >= margin & d > 0, "ruminantium", "unassigned"))
  cnt <- sapply(c("gottschalkii", "ruminantium", "unassigned"), function(k) rowSums(st[, j[cl == k], drop=FALSE]))
  rd_ <- colSums(st[, j, drop=FALSE]); list(cnt=data.frame(cow=rownames(st), cnt, check.names=FALSE), share=sum(rd_[cl != "unassigned"]) / sum(rd_), asv_assigned=sum(cl != "unassigned"), asv=length(cl)) }
# margin 1 must reproduce the clade counts on file
AS <- list(); CN <- list()
for (co in names(COH)) for (m in 1:3) { x <- counts_at(co, m); CN[[paste(co, m)]] <- x
  if (m == 1) { f <- rd(file.path(root, "profiles_16s", COH[[co]]$prof, paste0(COH[[co]]$prof, "_m1_counts.tsv"))); f$cow <- as.character(f$cow); k <- match(f$cow, x$cnt$cow)
    ok <- !anyNA(k) && all(f$gottschalkii == x$cnt$gottschalkii[k]) && all(f$ruminantium == x$cnt$ruminantium[k]); cat(co, "margin 1 reproduces the clade counts on file:", ok, "\n")
    if (!ok) stop("48.9: margin 1 does not reproduce profiles_16s/", COH[[co]]$prof, "/", COH[[co]]$prof, "_m1_counts.tsv") }
  AS[[length(AS) + 1]] <- data.frame(cohort=co, margin=m, share_assigned=x$share, asvs_assigned=x$asv_assigned, asvs=x$asv) }
AS <- do.call(rbind, AS)
AN <- list(list(i=1, name="Wallace CH4 per kg DMI", prefix="wallace_prokaryotes_ch4_dmi", co="wallace", noSE1=FALSE), list(i=2, name="Wallace CH4 g/d", prefix="wallace_prokaryotes", co="wallace", noSE1=FALSE),
  list(i=3, name="Difford CH4 g/d", prefix="difford_arch", co="difford", noSE1=FALSE), list(i=4, name="Wallace CH4 per kg DMI without SE1", prefix="wallace_prokaryotes_ch4_dmi", co="wallace", noSE1=TRUE))
GRID <- expand.grid(margin=1:3, pc=c(0.1, 0.5, 1), min_share=c(0.70, 0.80)); GRID$cell <- seq_len(nrow(GRID))
fit_cell <- function(an, g, keep_all=FALSE) { d <- a46_load(root, an$prefix); md <- d$md; if (an$noSE1) md <- md[as.character(md$farm) != "SE1", , drop=FALSE]
  x <- CN[[paste(an$co, g$margin)]]; sh <- x$share; c0 <- x$cnt[match(rownames(md), x$cnt$cow), ]; ok <- !is.na(c0$cow) & (c0$gottschalkii + c0$ruminantium) >= 100
  base <- data.frame(analysis=an$name, cell=g$cell, margin=g$margin, pseudo_count=g$pc, min_share=g$min_share, share_assigned=sh)
  if (sh < g$min_share) return(list(row=cbind(base, testable=FALSE, cows=NA, df=NA, t=NA, r=NA, lo=NA, hi=NA, p_t=NA, p_fl=NA, seed=NA)))
  mdk <- md[ok, , drop=FALSE]; y <- log((c0$gottschalkii[ok] + g$pc) / (c0$ruminantium[ok] + g$pc)); D <- a46_design(mdk, d$terms); E <- a46_resid(matrix(y, 1), D)
  seed <- 4890000 + 100 * an$i + g$cell; f <- a46_fl(E, D$rx, D$qz, interaction(mdk[, d$perm_by, drop=FALSE], drop=TRUE), B, seed=seed); ci <- a46_ci(f$r, D$df)
  list(row=cbind(base, testable=TRUE, cows=D$n, df=D$df, t=r2t(f$r, D$df), r=f$r, lo=ci[, "lo"], hi=ci[, "hi"], p_t=2 * pt(-abs(r2t(f$r, D$df)), D$df), p_fl=f$p_fl, seed=seed),
    y=y, md=mdk, terms=d$terms) }
GR <- do.call(rbind, lapply(AN, function(an) do.call(rbind, lapply(seq_len(nrow(GRID)), function(i) fit_cell(an, GRID[i, ])$row))))
# influence, analysis 1 at the current cell (margin 1, pseudo-count 0.5, share 70%)
cur <- fit_cell(AN[[1]], GRID[GRID$margin == 1 & GRID$pc == 0.5 & GRID$min_share == 0.70, ]); mm <- model.matrix(reformulate(cur$terms), droplevels(cur$md)); tv <- grep(paste0("^", cur$terms[length(cur$terms)]), colnames(mm))
fit <- lm(cur$y ~ mm - 1); cd <- cooks.distance(fit); n <- length(cur$y); dropc <- cd > 4 / n; fit2 <- lm(cur$y[!dropc] ~ mm[!dropc, , drop=FALSE] - 1)
st2 <- summary(fit2)$coefficients[tv, ]; rr <- tryCatch(robustbase::lmrob(cur$y ~ mm - 1), error=function(e) NULL); str <- if (is.null(rr)) rep(NA, 4) else summary(rr)$coefficients[tv, ]
df2 <- fit2$df.residual; r2 <- st2[["t value"]] / sqrt(st2[["t value"]]^2 + df2); ci2 <- a46_ci(r2, df2)
IN <- data.frame(analysis=AN[[1]]$name, fit=c("OLS, current cell", "OLS without Cook's D > 4/n", "MM regression (robustbase::lmrob, defaults)"),
  cows=c(n, sum(!dropc), n), removed=c(0, sum(dropc), 0), max_cooks_d=c(max(cd), NA, NA), estimate=c(coef(fit)[tv], st2[["Estimate"]], str[1]), t=c(summary(fit)$coefficients[tv, "t value"], st2[["t value"]], str[3]),
  p=c(summary(fit)$coefficients[tv, "Pr(>|t|)"], st2[["Pr(>|t|)"]], str[4]), r=c(cur$row$r, r2, NA), lo=c(cur$row$lo, ci2[1, "lo"], NA), hi=c(cur$row$hi, ci2[1, "hi"], NA),
  excludes_0=c(cur$row$lo > 0 | cur$row$hi < 0, ci2[1, "lo"] > 0 | ci2[1, "hi"] < 0, if (is.null(rr)) NA else str[4] < 0.05))
y1 <- GR[GR$analysis == AN[[1]]$name, ]; ex <- y1$testable & (y1$lo > 0 | y1$hi < 0)
RD <- data.frame(analysis=AN[[1]]$name, cells=nrow(y1), cells_testable=sum(y1$testable), cells_interval_excludes_0=sum(ex, na.rm=TRUE), influence_refits_excluding_0=sum(IN$excludes_0[-1], na.rm=TRUE),
  robust=all(ex %in% TRUE) && all(IN$excludes_0[-1] %in% TRUE), failing_cells=paste(y1$cell[!(ex %in% TRUE)], collapse=";"),
  reading=if (all(ex %in% TRUE) && all(IN$excludes_0[-1] %in% TRUE)) "robust within Wallace (still a single-cohort exploratory association)" else "not robust within Wallace: failing cells or refits listed (still a single-cohort exploratory association)",
  note="an untestable cell (assigned share below the minimum) counts as not excluding 0")
a48_write(AS, outs[1]); a48_write(GR, outs[2]); a48_write(IN, outs[3]); a48_write(RD, outs[4])
cat("== assignment\n"); print(AS, digits=3, row.names=FALSE); cat("\n== grid\n"); print(GR[, c("analysis", "cell", "margin", "pseudo_count", "min_share", "testable", "cows", "r", "lo", "hi", "p_fl")], digits=3, row.names=FALSE)
cat("\n== influence\n"); print(IN, digits=3, row.names=FALSE); cat("\n== reading\n"); print(RD, row.names=FALSE)
