#!/usr/bin/env Rscript
# Amendment 40: leave-one-farm-out and per-farm fits for the Wallace methane genera. Usage: 94_wallace_farms.R. No phenotype value is printed.
suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment); library(metafor)})
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); set.seed(20261110); plan <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md"))
stopifnot("commit Amendment 40 first" = any(grepl("^\\*\\*Amendment 40 \\(", plan)))
od <- file.path(root, "results", paste0("wallace_farms", if (Sys.getenv("DEPTH") == "1") "_depth" else "")); dir.create(od, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote=""); rc <- function(f) read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
coh <- "wallace_prokaryotes"
X <- rd(file.path(root, "profiles_16s", coh, paste0(coh, "_genus_counts.tsv"))); rownames(X) <- X[[1]]; X <- as.matrix(X[, -1, drop=FALSE])
lin <- rd(file.path(root, "results", "16s", coh, "genus_lineage_silva.tsv"))
an <- rd(file.path(root, "harmonized", "wallace_animals.tsv")); an$cow <- as.character(an$cow); an <- an[an$sequenced %in% c("True", TRUE, "TRUE"), ]
an$y <- suppressWarnings(as.numeric(an$ch4_gd)); an$dim <- suppressWarnings(as.numeric(an$dim))
an <- an[!is.na(an$y) & !is.na(an$dim) & an$parity %in% c("1", "2", "3+") & an$cow %in% colnames(X), ]
an$CH4_z <- NA_real_; for (g in unique(an$farm)) { i <- which(an$farm == g); an$CH4_z[i] <- as.numeric(scale(an$y[i])) }
md <- data.frame(farm=factor(an$farm), parity=factor(an$parity), dim=an$dim, CH4_z=an$CH4_z, row.names=an$cow); X <- X[, rownames(md)]; DEP <- Sys.getenv("DEPTH") == "1"; DT <- if (DEP) "log_depth" else NULL; if (DEP) md$log_depth <- log10(colSums(X))
fam <- lin$family[match(rownames(X), lin$genus)]; m2 <- rownames(X)[fam %in% "Succinivibrionaceae"]; A <- rbind(X, TGT_M2=colSums(X[m2, , drop=FALSE]))
# genera fixed by Amendment 40 item 2
sc <- rd(file.path(root, "results", "denovo_scan", "wallace_prokaryotes_scan.tsv")); gen <- c(sc$genus[sc$q_perm < 0.05], "TGT_M2", "Lachnospira"); gen <- unique(gen[gen %in% rownames(A)]); if (Sys.getenv("DEPTH") == "1") { s2 <- rd(file.path(root, "results", "denovo_scan", "wallace_prokaryotes_depth_scan.tsv")); gen <- unique(c(gen, s2$genus[s2$q_perm < 0.05])); gen <- gen[gen %in% rownames(A)] }
cat("genera under test:", length(gen), "| cows", ncol(A), "| farms", nlevels(md$farm), "\n")
fit1 <- function(A, md, terms) { orig <- rownames(A); rownames(A) <- sprintf("F%05d", seq_along(orig))
  tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(A)), colData=S4Vectors::DataFrame(md[, terms, drop=FALSE]))
  r <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH", prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)$res
  r$genus <- orig[match(r$taxon, rownames(A))]; data.frame(genus=r$genus, lfc=r$lfc_CH4_z, se=r$se_CH4_z, W=r$lfc_CH4_z / r$se_CH4_z) }
full <- fit1(A, md, c("farm", "parity", "dim", DT, "CH4_z")); write.table(full, file.path(od, "full.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
farms <- levels(md$farm); lofo <- list(); perfarm <- list()
for (f in farms) { k <- md$farm != f; m <- droplevels(md[k, ]); cat("leave out", f, ": cows", sum(k), "\n"); r <- fit1(A[, k], m, c("farm", "parity", "dim", DT, "CH4_z")); r$drop <- f; lofo[[f]] <- r
  k <- md$farm == f; m <- droplevels(md[k, ]); tm <- c(if (nlevels(m$parity) > 1) "parity", "dim", DT, "CH4_z"); cat("within", f, ": cows", sum(k), "| terms", paste(tm, collapse=" + "), "\n")
  r <- tryCatch(fit1(A[, k], m, tm), error=function(e) { cat("  within-farm fit failed:", conditionMessage(e), "\n"); NULL }); if (!is.null(r)) { r$farm <- f; perfarm[[f]] <- r } }
L <- do.call(rbind, lofo); P <- do.call(rbind, perfarm); write.table(L, file.path(od, "leave_one_farm_out.tsv"), sep="\t", quote=FALSE, row.names=FALSE); write.table(P, file.path(od, "per_farm.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# per-genus null 95th percentile of |W| from the 200 permutations of the full fit (Amendment 37 nulls)
fs <- list.files(file.path(root, "results", "methane_dryrun"), paste0("^wallace_prokaryotes", if (Sys.getenv("DEPTH") == "1") "_depth" else "", "_seed1[0-9]{3}_ancombc2\\.csv$"), full.names=TRUE); stopifnot(length(fs) >= 190)
NW <- do.call(rbind, lapply(fs, function(f) { r <- rc(f); nm <- if ("taxon_name" %in% names(r)) r$taxon_name else r$taxon; data.frame(genus=nm, W=r$lfc_CH4_z / r$se_CH4_z) }))
q95 <- tapply(abs(NW$W), NW$genus, quantile, 0.95, na.rm=TRUE)
out <- do.call(rbind, lapply(gen, function(g) { w0 <- full$W[full$genus == g]; if (!length(w0)) return(NULL); l <- L[L$genus == g, ]; p <- P[P$genus == g, ]
  i2 <- if (nrow(p) >= 3) tryCatch(rma(yi=p$lfc, sei=p$se, method="REML")$I2, error=function(e) NA) else NA
  data.frame(genus=g, W_full=round(w0, 2), null_q95=round(unname(q95[g]), 2), W_lofo_min=round(min(l$W), 2), W_lofo_max=round(max(l$W), 2), lofo_sign_kept=all(sign(l$W) == sign(w0)),
    farm_robust=all(sign(l$W) == sign(w0)) && min(abs(l$W)) > unname(q95[g]), farms_fit=nrow(p), farms_same_sign=sum(sign(p$W) == sign(w0)), I2=round(i2, 0),
    weakest_without=l$drop[which.min(abs(l$W))], stringsAsFactors=FALSE) }))
out <- out[order(-abs(out$W_full)), ]; write.table(out, file.path(od, "farm_consistency.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== WALLACE FARM CONSISTENCY (Amendment 40)\n"); print(out, row.names=FALSE)
cat("== farm-robust (sign kept in all 7 leave-one-farm-out refits and |W| stays above the per-genus null 95th percentile):", sum(out$farm_robust), "of", nrow(out), "\n")
cat("== sign agreement in per-farm fits, genera with >= 5 of 7 farms agreeing:", sum(out$farms_same_sign >= 5), "of", nrow(out), "\n")
