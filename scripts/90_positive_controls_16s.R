#!/usr/bin/env Rscript
# Amendment 38 item 2: positive controls for the 16S cohorts, same animals and model family as the real fits (65, 81), blocking factor replaced by the control contrast.
# Usage: 90_positive_controls_16s.R <alberta|difford_bact|difford_arch|wallace_prokaryotes>. No phenotype value is printed.
suppressMessages({library(vegan); library(zCompositions); library(ANCOMBC); library(TreeSummarizedExperiment)})
a <- commandArgs(TRUE); coh <- a[1]; stopifnot(coh %in% c("alberta", "difford_bact", "difford_arch", "wallace_prokaryotes")); set.seed(20261108)
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); plan <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md"))
stopifnot("commit Amendment 38 first" = any(grepl("^\\*\\*Amendment 38 \\(", plan)))
out <- file.path(root, "results", "positive_controls", "16s"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
X <- rd(file.path(root, "profiles_16s", coh, paste0(coh, "_genus_counts.tsv"))); rownames(X) <- X[[1]]; X <- as.matrix(X[, -1, drop=FALSE])
if (coh == "alberta") {   # selection copied from 65 (primary: overlap IDs excluded, management groups < 5 dropped)
  runs <- rd(file.path(root, "harmonized/alberta_runs.tsv")); an <- rd(file.path(root, "harmonized/alberta_animals.tsv"))
  aid <- runs$animal_id[match(colnames(X), runs$run_accession)]; stopifnot(!anyNA(aid), !anyDuplicated(aid)); colnames(X) <- aid
  an <- an[an$has_rfi %in% c(TRUE, "True") & an$animal_id %in% colnames(X), ]; ov <- readLines(file.path(root, "harmonized", "alberta_overlap_448333.txt")); an <- an[!an$animal_id %in% ov, ]
  an$group <- paste(an$diet, an$sex, an$collection_date, sep="|"); gs <- table(an$group); an <- an[!an$group %in% names(gs)[gs < 5], ]
  md <- data.frame(diet=factor(an$diet), breed=factor(an$breed), RFI_z=as.numeric(scale(as.numeric(an$rfi))), row.names=an$animal_id)
  terms <- c("diet", "breed", "RFI_z"); ctrl <- "diet"; lab <- "diet"
} else {   # selection copied from 81 (trait ch4_gd, all farms)
  src <- if (startsWith(coh, "difford")) "difford" else "wallace"; an <- rd(file.path(root, "harmonized", paste0(src, "_animals.tsv"))); an$cow <- as.character(an$cow)
  if (src == "wallace") an <- an[an$sequenced %in% c("True", TRUE, "TRUE"), ]
  an$y <- suppressWarnings(as.numeric(an$ch4_gd)); an$dim <- suppressWarnings(as.numeric(an$dim))
  an <- an[!is.na(an$y) & !is.na(an$dim) & an$parity %in% c("1", "2", "3+") & an$cow %in% colnames(X), ]
  grp <- if (src == "difford") an$herd else an$farm
  if (src == "difford") an$CH4_z <- as.numeric(scale(an$y)) else { an$CH4_z <- NA_real_; for (g in unique(grp)) { i <- which(grp == g); an$CH4_z[i] <- as.numeric(scale(an$y[i])) } }
  if (src == "difford") {
    md <- data.frame(herd=factor(an$herd), parity=factor(an$parity), dim=an$dim, batch=factor(an$batch), CH4_z=an$CH4_z, row.names=an$cow)
    mm <- model.matrix(~ herd + parity + dim + batch + CH4_z, md)
    if (qr(mm)$rank < ncol(mm)) { g <- paste(an$herd, an$batch, sep="|"); gs <- table(g); keep <- !g %in% names(gs)[gs < 5]; md <- md[keep, ]; cat("same animals as the real fit: dropped", sum(!keep), "cows in herd x batch groups < 5\n") }
    md$batch <- NULL; terms <- c("herd", "parity", "dim", "CH4_z"); ctrl <- "herd"; lab <- "herd-or-batch"
  } else { md <- data.frame(farm=factor(an$farm), parity=factor(an$parity), dim=an$dim, CH4_z=an$CH4_z, row.names=an$cow); terms <- c("farm", "parity", "dim", "CH4_z"); ctrl <- "farm"; lab <- "farm-or-breed" }
}
X <- X[, rownames(md), drop=FALSE]; X <- X[rowMeans(X > 0) >= 0.10, , drop=FALSE]
mm <- model.matrix(reformulate(terms), md); stopifnot("design not of full rank" = qr(mm)$rank == ncol(mm))
cat(coh, ": animals", ncol(X), "| genera (>=10% prevalence)", nrow(X), "| control", ctrl, "(", nlevels(md[[ctrl]]), "levels ) | model", paste(terms, collapse=" + "), "\n")
Z <- suppressMessages(cmultRepl(t(X), method="CZM", label=0, output="prop", z.warning=1, z.delete=FALSE, suppress.print=TRUE))
clr <- t(apply(Z, 1, function(p) log(p) - mean(log(p)))); d <- dist(clr)
pm <- adonis2(as.formula(paste("d ~", paste(terms, collapse=" + "))), data=md, by="margin", permutations=1999)
bd <- permutest(betadisper(d, md[[ctrl]]), permutations=1999)$tab[1, "Pr(>F)"]
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(X)), colData=S4Vectors::DataFrame(md[, terms, drop=FALSE]))
fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH", prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
res <- fit$res; write.csv(res, file.path(out, paste0(coh, "_control_ancombc2.csv")), row.names=FALSE)
qc <- grep(paste0("^q_", ctrl), names(res), value=TRUE); stopifnot(length(qc) >= 1)
qmin <- apply(as.matrix(res[qc]), 1, min, na.rm=TRUE)
s <- data.frame(table=coh, control=ctrl, label=lab, n_animals=ncol(X), n_genera=nrow(res), n_levels=length(qc) + 1L, n_q05=sum(qmin < 0.05), min_q=signif(min(qmin), 3),
  permanova_R2=round(pm[ctrl, "R2"], 3), permanova_p=pm[ctrl, "Pr(>F)"], betadisper_p=bd, pass=sum(qmin < 0.05) >= 1)
write.table(s, file.path(out, paste0(coh, "_summary.tsv")), sep="\t", quote=FALSE, row.names=FALSE); print(s, row.names=FALSE)
cat("POSITIVE CONTROL", if (s$pass) "PASS" else "FAIL", "\n")
