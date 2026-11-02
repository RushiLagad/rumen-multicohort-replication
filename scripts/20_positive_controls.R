#!/usr/bin/env Rscript
# Positive controls (plan Sec. 7): Nelore diet; Charolais diet and year. Usage: 20_positive_controls.R <profile_dir> <nelore|charolais_mrna>
suppressMessages({library(vegan); library(zCompositions); library(ANCOMBC); library(TreeSummarizedExperiment)})
a <- commandArgs(TRUE); prof <- a[1]; coh <- a[2]; set.seed(20261102)
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")
out <- file.path(root, "results", "positive_controls", prof); dir.create(out, recursive=TRUE, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
meta <- rd(file.path(root, "harmonized", "animal_metadata.tsv")); runs <- rd(file.path(root, "harmonized", "runs.tsv"))
cfg <- list(nelore=list(c="Nelore_PRJNA987743", l="DNA_shotgun"), charolais_mrna=list(c="Charolais_PRJEB70895", l="mRNA"))[[coh]]
X <- rd(file.path(root, prof, coh, paste0(coh, "_genus_counts.tsv"))); rownames(X) <- X[[1]]; X <- X[, -1]
colnames(X) <- sub("_clean$", "", colnames(X))
r <- runs[runs$cohort == cfg$c & runs$layer == cfg$l, ]
stopifnot(all(colnames(X) %in% r$run_accession))
aid <- as.character(r$animal_id[match(colnames(X), r$run_accession)])
X <- round(sapply(split(seq_along(aid), aid), function(i) rowSums(X[, i, drop=FALSE])))  # several runs per animal: counts summed
md <- meta[meta$cohort == cfg$c, ]; md <- md[match(colnames(X), as.character(md$animal_id)), ]
stopifnot(!anyNA(md$animal_id)); rownames(md) <- colnames(X)
if (coh == "nelore") {
  md$diet <- factor(md$diet, levels=c("Conventional", "By-products")); md$batch <- factor(md$batch)
  cat("rfi_value sd =", round(sd(md$rfi_value), 3), "\n"); md$RFI_z <- as.numeric(scale(md$rfi_value))
  terms <- c("diet", "batch", "RFI_z"); ctrl <- "diet"
} else {
  md$diet <- factor(md$diet, levels=c("CS", "GS")); md$year <- factor(sub("batch_", "", md$batch))
  md$RFI_class <- factor(md$rfi_class, levels=c("low_RFI", "high_RFI"))
  terms <- c("diet", "year", "RFI_class"); ctrl <- c("diet", "year")
}
X <- X[rowMeans(X > 0) >= 0.10, ]   # plan prevalence filter: present in >=10% of samples
cat(coh, ":", ncol(X), "samples,", nrow(X), "genera after prevalence filter\n")
Z <- suppressMessages(cmultRepl(t(X), method="CZM", label=0, output="prop", z.warning=1, z.delete=FALSE, suppress.print=TRUE))
clr <- t(apply(Z, 1, function(p) log(p) - mean(log(p)))); d <- dist(clr)   # Aitchison
pm <- adonis2(as.formula(paste("d ~", paste(terms, collapse=" + "))), data=md, by="margin", permutations=9999)
bd <- sapply(ctrl, function(v) permutest(betadisper(d, md[[v]]), permutations=9999)$tab[1, "Pr(>F)"])
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(X)), colData=S4Vectors::DataFrame(md[, terms, drop=FALSE]))
fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH",
                prv_cut=0.10, group=ctrl[1], struc_zero=TRUE, verbose=FALSE)
res <- fit$res; write.csv(res, file.path(out, paste0(coh, "_ancombc2.csv")), row.names=FALSE)
s <- do.call(rbind, lapply(ctrl, function(v) {
  qc <- grep(paste0("^q_", v), names(res), value=TRUE); if (!length(qc)) stop("no q_ column for ", v, ": ", paste(names(res), collapse=", "))
  q <- unlist(res[qc])
  data.frame(profile=prof, cohort=coh, term=v, n_genera=nrow(res), n_q05=sum(q < 0.05, na.rm=TRUE), min_q=signif(min(q, na.rm=TRUE), 3),
             permanova_R2=round(pm[v, "R2"], 3), permanova_p=pm[v, "Pr(>F)"], betadisper_p=unname(bd[v]))
}))
write.table(s, file.path(out, paste0(coh, "_summary.tsv")), sep="\t", quote=FALSE, row.names=FALSE); print(s, row.names=FALSE)
cat("POSITIVE CONTROL", if (any(s$n_q05 > 0)) "PASS" else "FAIL", "(>=1 genus at BH q<0.05 for:", paste(ctrl, collapse=" or "), ")\n")
