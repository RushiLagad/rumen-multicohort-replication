#!/usr/bin/env Rscript
# Amendments 19, 23 and 35 (SILVA families): Alberta (PRJNA393057, 16S) target estimates. Usage: 65_estimate_alberta.R <genus|family> [dry]
a <- commandArgs(TRUE); rk <- a[1]; dry <- length(a) > 1 && a[2] == "dry"; stopifnot(rk %in% c("genus", "family"))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); set.seed(20261103)
out <- file.path(root, "results", "alberta"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("commit Amendment 24 first" = any(grepl("^\\*\\*Amendment 24 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
X <- rd(file.path(root, "profiles_16s/alberta/alberta_genus_counts.tsv")); rownames(X) <- X[[1]]; X <- as.matrix(X[, -1, drop=FALSE])
runs <- rd(file.path(root, "harmonized/alberta_runs.tsv")); an <- rd(file.path(root, "harmonized/alberta_animals.tsv"))
aid <- runs$animal_id[match(colnames(X), runs$run_accession)]; stopifnot(!anyNA(aid), !anyDuplicated(aid)); colnames(X) <- aid
an <- an[an$has_rfi %in% c(TRUE, "True") & an$animal_id %in% colnames(X), ]
allA <- Sys.getenv("ALBERTA_ALL") == "1"; tag <- if (allA) "_all" else ""; ovf <- file.path(root, "harmonized", "alberta_overlap_448333.txt")
stopifnot("overlap list missing (Amendment 24)" = file.exists(ovf)); ov <- readLines(ovf); nov <- sum(an$animal_id %in% ov); if (!allA) an <- an[!an$animal_id %in% ov, ]
an$group <- paste(an$diet, an$sex, an$collection_date, sep="|"); gs <- table(an$group); small <- names(gs)[gs < 5]
write.table(data.frame(group=names(gs), n=as.integer(gs), excluded=names(gs) %in% small), file.path(out, paste0("management_groups", tag, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
an <- an[!an$group %in% small, ]; X <- X[, an$animal_id, drop=FALSE]
md <- data.frame(group=factor(an$group), breed=factor(an$breed), RFI_z=as.numeric(scale(as.numeric(an$rfi))), row.names=an$animal_id)
terms <- c("group", "breed", "RFI_z"); mm <- model.matrix(reformulate(terms), md)
stopifnot("Alberta design not of full rank" = qr(mm)$rank == ncol(mm)); n <- ncol(X); df <- n - ncol(mm)
cat("alberta", rk, if (allA) "(sensitivity: all animals)" else "(primary: PRJNA448333 ID numbers excluded)", "| overlap IDs in analysis set", nov, ": animals", n, "| groups", nlevels(md$group), "(excluded", length(small), "with <5 animals:", sum(gs[small]), "animals) | df", df, "\n")
if (rk == "family") {
  lin <- rd(file.path(root, "results", "16s", "alberta", "genus_lineage_silva.tsv")); grp <- lin$family[match(rownames(X), lin$genus)]
  grp[is.na(grp)] <- "unassigned_family"; A <- rowsum(X, grp)
} else A <- X
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
tk <- tk[tk$target_rank == rk & tk$status == "mapped" & tk$n_keys > 0 & grepl("^RFI", tk$phenotype), ]
keys <- lapply(strsplit(tk$gtdb_keys, ";"), function(x) { x <- trimws(x); sort(unique(x[nzchar(x)])) })
pres <- lapply(keys, intersect, rownames(A)); psig <- vapply(pres, paste, "", collapse=";")
ok <- nzchar(psig); feat <- ifelse(ok, paste0("TGT_", tk$id[match(psig, psig)]), NA)
u <- unique(psig[ok]); allk <- unlist(strsplit(u, ";"))
if (anyDuplicated(allk)) stop("Partial key overlap between targets: ", paste(unique(allk[duplicated(allk)]), collapse=", "))
for (s in u) { k <- strsplit(s, ";")[[1]]; f <- paste0("TGT_", tk$id[match(s, psig)])
  A <- rbind(A[!rownames(A) %in% k, , drop=FALSE], colSums(A[k, , drop=FALSE])); rownames(A)[nrow(A)] <- f }
m <- A[match(feat, rownames(A)), , drop=FALSE]; prev <- rowMeans(m > 0)
cov <- data.frame(id=tk$id, target=tk$target, set=tk$set, n_keys=tk$n_keys, n_present=lengths(pres), feature=feat,
  reads=unname(rowSums(m)), prevalence=round(unname(prev), 3), pass10=unname(prev) >= 0.10)
cat("targets", nrow(tk), "| matched a 16S feature", sum(ok), "| pass 10% prevalence", sum(cov$pass10, na.rm=TRUE), "| features", nrow(A), "\n"); print(cov, row.names=FALSE)
if (dry) quit(save="no")
suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment)})
orig <- rownames(A); rownames(A) <- sprintf("F%05d", seq_along(orig)); sf <- rownames(A)[match(feat, orig)]
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(A)), colData=md[, terms, drop=FALSE])
fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH",
                prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
res <- fit$res; stopifnot(nrow(res) > 0, all(res$taxon %in% rownames(A))); cat("ANCOM-BC2 kept", nrow(res), "of", nrow(A), "features\n"); i <- match(sf, res$taxon)
res$taxon_name <- orig[match(res$taxon, rownames(A))]; write.csv(res, file.path(out, paste0("alberta_", rk, tag, "_ancombc2.csv")), row.names=FALSE)
tc <- grep("^lfc_RFI", names(res), value=TRUE); stopifnot(length(tc) == 1); sfx <- sub("^lfc_", "", tc)
est <- data.frame(id=tk$id, set=tk$set, site=tk$site, cohort=paste0("alberta", tag), rank=rk, feature=feat, n=n, df=df,
  lfc=res[[tc]][i], se=res[[paste0("se_", sfx)]][i], p_ancom=res[[paste0("p_", sfx)]][i])
est$W <- est$lfc / est$se
est$status <- ifelse(is.na(feat), "no_feature", ifelse(is.na(i), "filtered", ifelse(is.finite(est$W), "ok", "no_estimate")))
write.table(est, file.path(out, paste0("alberta_", rk, tag, "_targets.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
print(est[, c("id", "feature", "status", "lfc", "se", "W", "p_ancom")], row.names=FALSE)
