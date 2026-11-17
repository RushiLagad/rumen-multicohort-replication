#!/usr/bin/env Rscript
# Amendment 48.5: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. CLR transformation sensitivity.
# 12 cells: zero handling (pseudo-count 0.1, 0.5, 1; multiplicative replacement, zCompositions::cmultRepl CZM) x reference (R0 current = non-target features after the
# 10% prevalence filter; R1 = every filtered feature including targets; R2 = non-target features present in >= 90% of that table's animals). Model, covariates and design
# as 46.1; M1 (a log ratio) unchanged. Cell (pseudo-count 0.5, R0) must reproduce results/late_a46/effect_size/cohort_estimates.tsv to 1e-8, or the run stops.
# Per cell: per-cohort r and CI, mKH pooling (k >= 2), exclusion of 0.15, 0.30 and 0, CLR label (47.2 rules, parametric p; all 38 targets refitted for the BH families).
# Usage: 132_a48_clr_grid.R. Writes results/late_a48/clr_grid/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "clr_grid"); outs <- file.path(od, c("cohort_cells.tsv", "pooled_cells.tsv", "changes.tsv", "robustness.tsv")); a48_fresh(outs)
GRID <- A48_GRID; if (Sys.getenv("A48_TEST_NO_CZM") == "1") { GRID <- GRID[GRID$zero != "czm", ]; cat("NOTE: A48_TEST_NO_CZM (test run; the plan has 12 cells)\n") } else
  stopifnot("zCompositions is required for the multiplicative-replacement cells (48.5)" = requireNamespace("zCompositions", quietly=TRUE))
SCOPE <- unique(c(a48_seven, "C5", "C6", "N02", "N04", "N11", "N13", "N20", "N22", "M2", "M3", "N32", "N38", "N39", "N40", "N06", "C3"))
cat(A48_LABEL, "| 48.5 CLR grid | cells", nrow(GRID), "| targets reported", length(SCOPE), "\n")
dirs <- a48_dirs(root); ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv"))
CO <- list(); PO <- list()
for (i in seq_len(nrow(GRID))) { z <- GRID$zero[i]; rf <- GRID$ref[i]; cell <- paste(z, rf); t0 <- proc.time()[[3]]
  U <- a48_units(root, clrfun=function(A) a48_clr_cell(A, z, rf)); est <- a48_est(U)
  if (z == "pc0.5" && rf == "R0") { j <- match(paste(ce$id, ce$cohort), paste(est$id, est$cohort)); stopifnot("cohort rows differ from cohort_estimates.tsv" = !anyNA(j), nrow(est) == nrow(ce))
    dr <- max(abs(ce$r - est$r[j])); cat("reference cell reproduces cohort_estimates.tsv: max |dr| =", signif(dr, 3), "\n")
    if (dr > 1e-8) stop("48.5: cell (pseudo-count 0.5, R0) does not reproduce the Amendment 46 estimates (max |dr| ", signif(dr, 3), ")") }
  lab <- a48_labels(root, est); po <- a48_pool_targets(est, dirs)
  lab1 <- lab[!duplicated(lab$id), c("id", "label", "q_bh")]; po <- merge(po, lab1, by="id", all.x=TRUE)
  CO[[i]] <- cbind(cell=cell, zero=z, ref=rf, est[est$id %in% SCOPE, c("arm", "id", "cohort", "table", "n", "df", "r", "lo", "hi", "p_t", "n_ref")])
  PO[[i]] <- cbind(cell=cell, zero=z, ref=rf, po[po$id %in% SCOPE, ]); cat(cell, "done,", round(proc.time()[[3]] - t0), "s\n") }
CO <- do.call(rbind, CO); PO <- do.call(rbind, PO)
ref <- PO[PO$cell == "pc0.5 R0", ]; j <- match(PO$id, ref$id)
PO$sign_changed <- sign(PO$r) != sign(ref$r[j]); PO$excl15_changed <- PO$excl_15 != ref$excl_15[j]; PO$excl30_changed <- PO$excl_30 != ref$excl_30[j]
PO$excl0_changed <- PO$excl_0 != ref$excl_0[j]; PO$label_changed <- PO$label != ref$label[j]
CH <- PO[PO$sign_changed | PO$excl15_changed | PO$excl30_changed | PO$excl0_changed | PO$label_changed, c("cell", "id", "k", "r", "lo", "hi", "excl_0", "excl_15", "excl_30", "label",
  "sign_changed", "excl15_changed", "excl30_changed", "excl0_changed", "label_changed")]
if (!nrow(CH)) CH <- data.frame(cell=NA, id=NA, note="no sign, exclusion or label differs from the reference cell in any cell")
RB <- do.call(rbind, lapply(SCOPE, function(id) { x <- PO[PO$id == id, ]; r0 <- ref[ref$id == id, ]; if (!nrow(r0)) return(NULL)
  st <- if (isTRUE(r0$excl_30) || isTRUE(r0$excl_15)) { d <- if (isTRUE(r0$excl_15)) "excl_15" else "excl_30"; c(statement=paste("excludes", if (d == "excl_15") "r = 0.15" else "r = 0.30"), col=d) } else
    if (isTRUE(r0$excl_0)) c(statement="interval excludes 0", col="excl_0") else c(statement="no exclusion", col=NA)
  holds <- if (is.na(st[["col"]])) NA else all(x[[st[["col"]]]] %in% TRUE); fails <- if (is.na(st[["col"]])) "" else paste(x$cell[!(x[[st[["col"]]]] %in% TRUE)], collapse=";")
  data.frame(id=id, k=r0$k, statement_reference=st[["statement"]], cells=nrow(x), holds_in_all_cells=holds, failing_cells=fails,
    reading=if (is.na(holds)) "" else if (holds) "robust to the transformation" else "transformation-sensitive",
    r_min=min(x$r), r_max=max(x$r), bound_min=min(x$bound), bound_max=max(x$bound), labels=paste(unique(x$label), collapse=";")) }))
a48_write(CO, outs[1]); a48_write(PO, outs[2]); a48_write(CH, outs[3]); a48_write(RB, outs[4])
cat("== robustness of the reported statements across", length(unique(PO$cell)), "cells\n"); print(RB, digits=3, row.names=FALSE)
