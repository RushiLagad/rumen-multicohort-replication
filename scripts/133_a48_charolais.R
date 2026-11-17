#!/usr/bin/env Rscript
# Amendment 48.6: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Charolais influence.
# Targets: the six Charolais benchmark exclusions (N02, N04, N11, N13, N20, N22) and the Charolais intervals excluding 0 (N06, N13; C3 added in the reconciliation).
# Leave-one-animal-out: 30 refits, each recomputing the whole 46.1 pipeline on 29 animals (prevalence filter, CLR, covariate residuals, r and CI); range of r and of the bound
# in the published direction; refits in which the exclusion holds; the most influential animal by DFBETAS for the RFI-class term (full-data OLS of the target's CLR value
# on diet, year and RFI class), by position, never by identifier. The 48.5 grid for these targets is reported beside it (after 132). Reading as 48.6.
# Usage: 133_a48_charolais.R (after 132). Writes results/late_a48/charolais/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "charolais"); outs <- file.path(od, c("loao.tsv", "summary.tsv")); a48_fresh(outs)
gd <- file.path(a48_base(root), "clr_grid", "robustness.tsv"); stopifnot("run 132 first (48.5 grid)" = file.exists(gd))
EXCL <- c("N02", "N04", "N11", "N13", "N20", "N22"); ZERO <- c("N06", "N13", "C3"); TG <- unique(c(EXCL, ZERO)); MARGINAL <- c("N11", "N20", "N22")
cat(A48_LABEL, "| 48.6 Charolais influence | targets", paste(TG, collapse=", "), "\n")
dirs <- a48_dirs(root); ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv")); ce <- ce[ce$cohort == "charolais_mrna" & ce$id %in% TG, ]
stopifnot(all(TG %in% ce$id))
LO <- list(); SU <- list()
for (id in TG) { e <- ce[ce$id == id, ]; d <- a46_load(root, e$table); f <- d$feat[d$tk$id == id]; dr <- dirs[[id]] %||% NA; md <- d$md; n <- nrow(md)
  one <- function(keep) { A <- d$A[, rownames(md)[keep], drop=FALSE]; C <- a46_clr(A); if (!f %in% rownames(C)) return(c(r=NA, lo=NA, hi=NA, df=NA))
    D <- a46_design(md[keep, , drop=FALSE], d$terms); r <- unname(a46_r(a46_resid(C[f, , drop=FALSE], D), D$rx)); ci <- a46_ci(r, D$df); c(r=r, lo=unname(ci[1, "lo"]), hi=unname(ci[1, "hi"]), df=D$df) }
  full <- one(rep(TRUE, n)); stopifnot("full-data r differs from cohort_estimates.tsv" = abs(full[["r"]] - e$r) < 1e-8)
  L <- do.call(rbind, lapply(seq_len(n), function(j) { o <- one(seq_len(n) != j); data.frame(id=id, left_out_position=j, r=o[["r"]], lo=o[["lo"]], hi=o[["hi"]], df=o[["df"]]) }))
  L$bound <- a46_bound(L$lo, L$hi, dr); L$excl_15 <- a47_excl(L$lo, L$hi, dr, 0.15); L$excl_30 <- a47_excl(L$lo, L$hi, dr, 0.30); L$excl_0 <- L$lo > 0 | L$hi < 0
  # DFBETAS for the RFI-class term, full-data OLS of the target's CLR value
  C <- a46_clr(d$A[, rownames(md), drop=FALSE]); mm <- model.matrix(reformulate(d$terms), droplevels(md)); y <- as.vector(C[f, ]); fit <- lm(y ~ mm - 1)
  jb <- grep("RFI", colnames(mm)); db <- dfbetas(fit)[, jb]; L$dfbetas_rfi_class <- db
  fl <- c(excl_15=a47_excl(full[["lo"]], full[["hi"]], dr, 0.15), excl_30=a47_excl(full[["lo"]], full[["hi"]], dr, 0.30), excl_0=full[["lo"]] > 0 | full[["hi"]] < 0)
  st <- if (id %in% EXCL) { if (fl[["excl_15"]]) "excl_15" else "excl_30" } else "excl_0"
  hold <- sum(L[[st]] %in% TRUE); g <- rd(gd); g <- g[g$id == id, ]
  SU[[id]] <- data.frame(id=id, table=e$table, published_dir=dr, r_full=full[["r"]], lo_full=full[["lo"]], hi_full=full[["hi"]], bound_full=a46_bound(full[["lo"]], full[["hi"]], dr),
    statement=c(excl_15="excludes r = 0.15", excl_30="excludes r = 0.30", excl_0="interval excludes 0")[[st]], r_min=min(L$r, na.rm=TRUE), r_max=max(L$r, na.rm=TRUE),
    bound_min=min(L$bound, na.rm=TRUE), bound_max=max(L$bound, na.rm=TRUE), refits=nrow(L), refits_holding=hold, robust_to_single_animals=hold == nrow(L),
    most_influential_position=which.max(abs(db)), max_abs_dfbetas=max(abs(db)), grid_holds_in_all_cells=if (nrow(g)) g$holds_in_all_cells else NA, grid_failing_cells=if (nrow(g)) g$failing_cells else NA,
    at_threshold=id %in% MARGINAL, note=if (id %in% MARGINAL) "bound within 0.01 of r = 0.15: described as at the threshold whatever the result" else "")
  LO[[id]] <- L }
LO <- do.call(rbind, LO); SU <- do.call(rbind, SU)
a48_write(LO, outs[1]); a48_write(SU, outs[2])
cat("== 48.6 summary\n"); print(SU[, c("id", "statement", "r_full", "bound_full", "bound_min", "bound_max", "refits_holding", "robust_to_single_animals", "most_influential_position", "grid_holds_in_all_cells", "at_threshold")], digits=3, row.names=FALSE)
