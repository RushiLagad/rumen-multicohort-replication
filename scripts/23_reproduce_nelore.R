#!/usr/bin/env Rscript
# Level 1 reproducibility (Amendment 15): Nelore-set targets re-tested in Nelore itself. Usage: 23_reproduce_nelore.R <RFI|RME> [dry]
suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment)})
a <- commandArgs(TRUE); phen <- a[1]; dry <- length(a) > 1 && a[2] == "dry"; stopifnot(phen %in% c("RFI", "RME"))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); prof <- Sys.getenv("PROFILE", "profiles_gtdb_full"); set.seed(20261102)
out <- file.path(root, "results", "reproducibility"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
meta <- rd(file.path(root, "harmonized", "animal_metadata.tsv")); runs <- rd(file.path(root, "harmonized", "runs.tsv"))
X <- rd(file.path(root, prof, "nelore", "nelore_genus_counts.tsv")); rownames(X) <- X[[1]]; X <- X[, -1]
colnames(X) <- sub("_clean$", "", colnames(X))
r <- runs[runs$cohort == "Nelore_PRJNA987743" & runs$layer == "DNA_shotgun", ]; stopifnot(all(colnames(X) %in% r$run_accession))
aid <- as.character(r$animal_id[match(colnames(X), r$run_accession)])
X <- round(sapply(split(seq_along(aid), aid), function(i) rowSums(X[, i, drop=FALSE])))
md <- meta[meta$cohort == "Nelore_PRJNA987743", ]; md <- md[match(colnames(X), as.character(md$animal_id)), ]
stopifnot(!anyNA(md$animal_id)); rownames(md) <- colnames(X)
pv <- md[[if (phen == "RFI") "rfi_value" else "methane_value"]]; stopifnot(!anyNA(pv))
md$diet <- factor(md$diet, levels=c("Conventional", "By-products")); md$batch <- factor(md$batch)
term <- paste0(phen, "_z"); md[[term]] <- as.numeric(scale(pv)); terms <- c("diet", "batch", term)
n <- ncol(X); df <- n - ncol(model.matrix(reformulate(terms), md)); A <- as.matrix(X)
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
tk <- tk[tk$set == "Nelore" & tk$status == "mapped" & tk$n_keys > 0 & grepl(paste0("^", phen), tk$phenotype), ]
cat("nelore genus", phen, ": targets", nrow(tk), " features", nrow(A), " n", n, " df", df, "\n")
keys <- lapply(strsplit(tk$gtdb_keys, ";"), function(x) { x <- trimws(x); sort(unique(x[nzchar(x)])) })
pres <- lapply(keys, intersect, rownames(A)); psig <- vapply(pres, paste, "", collapse=";")
ok <- nzchar(psig); feat <- ifelse(ok, paste0("TGT_", tk$id[match(psig, psig)]), NA)
u <- unique(psig[ok]); allk <- unlist(strsplit(u, ";"))
if (anyDuplicated(allk)) stop("Partial key overlap between targets: ", paste(unique(allk[duplicated(allk)]), collapse=", "))
for (s in u) { k <- strsplit(s, ";")[[1]]; f <- paste0("TGT_", tk$id[match(s, psig)])
  A <- rbind(A[!rownames(A) %in% k, , drop=FALSE], colSums(A[k, , drop=FALSE])); rownames(A)[nrow(A)] <- f }
m <- A[match(feat, rownames(A)), , drop=FALSE]; prev <- rowMeans(m > 0)
print(data.frame(id=tk$id, target=tk$target, site=tk$site, n_keys=tk$n_keys, n_present=lengths(pres), reads=rowSums(m), prevalence=round(prev, 3)), row.names=FALSE)
if (dry) quit(save="no")
orig <- rownames(A); rownames(A) <- sprintf("F%05d", seq_along(orig)); sf <- rownames(A)[match(feat, orig)]
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(A)), colData=md[, terms, drop=FALSE])
fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH",
                prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
res <- fit$res; stopifnot(nrow(res) > 0); cat("ANCOM-BC2 kept", nrow(res), "of", nrow(A), "features\n"); i <- match(sf, res$taxon)
res$taxon_name <- orig[match(res$taxon, rownames(A))]; write.csv(res, file.path(out, paste0("nelore_", phen, "_ancombc2.csv")), row.names=FALSE)
tc <- grep(paste0("^lfc_", phen), names(res), value=TRUE); stopifnot(length(tc) == 1); sfx <- sub("^lfc_", "", tc)
est <- data.frame(id=tk$id, set=tk$set, site=tk$site, phenotype=tk$phenotype, cohort="nelore", feature=feat, n=n, df=df,
  lfc=res[[tc]][i], se=res[[paste0("se_", sfx)]][i], p_ancom=res[[paste0("p_", sfx)]][i]); est$W <- est$lfc / est$se
est$status <- ifelse(is.na(feat), "no_feature", ifelse(is.na(i), "filtered", ifelse(is.finite(est$W), "ok", "no_estimate")))
write.table(est, file.path(out, paste0("nelore_", phen, "_targets.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
print(est[, c("id", "site", "status", "lfc", "se", "W", "p_ancom")], row.names=FALSE)
