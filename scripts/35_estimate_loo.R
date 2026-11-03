#!/usr/bin/env Rscript
# Target estimates (Amendment 14b). Usage: 21_estimate_targets.R <nelore|charolais_mrna> <genus|family|phylum|species> [dry]
suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment)})
a <- commandArgs(TRUE); coh <- a[1]; rk <- a[2]; dry <- length(a) > 2 && a[3] == "dry"
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); prof <- Sys.getenv("PROFILE", "profiles_gtdb_full")
tax <- Sys.getenv("TAXDIR", "/data/refs/k2_gtdb207_full/taxonomy"); set.seed(20261102)
out <- file.path(root, "results", if (coh == "charolais_16s") "charolais16s" else "replication"); if (nzchar(Sys.getenv("DROP"))) out <- file.path(root, "results", "loo", paste0("drop_", Sys.getenv("DROP"))); dir.create(out, recursive=TRUE, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
meta <- rd(file.path(root, "harmonized", "animal_metadata.tsv")); runs <- rd(file.path(root, "harmonized", "runs.tsv"))
cfg <- list(nelore=list(c="Nelore_PRJNA987743", l="DNA_shotgun"), charolais_mrna=list(c="Charolais_PRJEB70895", l="mRNA"), charolais_16s=list(c="Charolais_PRJEB70895", l="16S_amplicon"))[[coh]]
lev <- if (rk == "species") "species" else "genus"
X <- rd(file.path(root, prof, coh, paste0(coh, "_", lev, "_counts.tsv"))); rownames(X) <- X[[1]]; X <- X[, -1]
colnames(X) <- sub("_clean$", "", colnames(X))
r <- runs[runs$cohort == cfg$c & runs$layer == cfg$l, ]; stopifnot(all(colnames(X) %in% r$run_accession))
aid <- as.character(r$animal_id[match(colnames(X), r$run_accession)])
X <- round(sapply(split(seq_along(aid), aid), function(i) rowSums(X[, i, drop=FALSE])))
dr <- Sys.getenv("DROP"); if (nzchar(dr)) { stopifnot(dr %in% colnames(X)); X <- X[, colnames(X) != dr, drop=FALSE] }; md <- meta[meta$cohort == cfg$c, ]; md <- md[match(colnames(X), as.character(md$animal_id)), ]
stopifnot(!anyNA(md$animal_id)); rownames(md) <- colnames(X)
if (coh == "nelore") {
  md$diet <- factor(md$diet, levels=c("Conventional", "By-products")); md$batch <- factor(md$batch)
  md$RFI_z <- as.numeric(scale(md$rfi_value)); terms <- c("diet", "batch", "RFI_z")
} else {
  md$diet <- factor(md$diet, levels=c("CS", "GS")); md$year <- factor(sub("batch_", "", md$batch))
  md$RFI_class <- factor(md$rfi_class, levels=c("low_RFI", "high_RFI")); terms <- c("diet", "year", "RFI_class")
}
n <- ncol(X); df <- n - ncol(model.matrix(reformulate(terms), md))
# --- aggregate to rank (genus->family/phylum via taxonomy)
if (rk %in% c("family", "phylum")) {
  lf <- file.path(root, "results", "replication", "genus_lineage.tsv")
  if (!file.exists(lf)) {
    sp <- function(l) strsplit(l, "\t|\t", fixed=TRUE); col <- function(v, i) vapply(v, `[`, "", i)
    nd <- sp(readLines(file.path(tax, "nodes.dmp"))); id <- col(nd, 1)
    par <- setNames(col(nd, 2), id); rnk <- setNames(col(nd, 3), id)
    nl <- readLines(file.path(tax, "names.dmp")); nl <- sp(nl[grepl("scientific name", nl, fixed=TRUE)])
    nm <- setNames(col(nl, 2), col(nl, 1)); g <- names(rnk)[rnk == "genus"]
    climb <- function(t, tr) { for (i in 1:40) { go <- !is.na(t) & unname(rnk[t]) != tr & t != "1"
      if (!any(go, na.rm=TRUE)) break; t[go] <- unname(par[t[go]]) }; ifelse(!is.na(t) & unname(rnk[t]) == tr, unname(nm[t]), NA) }
    lin <- data.frame(genus=unname(nm[g]), family=climb(g, "family"), phylum=climb(g, "phylum"))
    write.table(lin[!duplicated(lin$genus), ], lf, sep="\t", quote=FALSE, row.names=FALSE)
  }
  lin <- rd(lf); grp <- lin[[rk]][match(rownames(X), lin$genus)]; grp[is.na(grp)] <- paste0("unassigned_", rk)
  A <- rowsum(X, grp)
} else A <- as.matrix(X)
# --- targets for this cohort and rank
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
sel <- tk$target_rank == rk & tk$status == "mapped" & tk$n_keys > 0 & grepl("^RFI", tk$phenotype)
if (coh == "nelore") sel <- sel & tk$set != "Nelore"
tk <- tk[sel, ]; cat(coh, rk, ": targets", nrow(tk), " features", nrow(A), " n", n, " df", df, "\n")
keys <- lapply(strsplit(tk$gtdb_keys, ";"), function(x) { x <- trimws(x); sort(unique(x[nzchar(x)])) })
pres <- lapply(keys, intersect, rownames(A)); psig <- vapply(pres, paste, "", collapse=";")
ok <- nzchar(psig); feat <- ifelse(ok, paste0("TGT_", tk$id[match(psig, psig)]), NA)
u <- unique(psig[ok]); allk <- unlist(strsplit(u, ";"))
if (anyDuplicated(allk)) stop("Partial key overlap between targets: ", paste(unique(allk[duplicated(allk)]), collapse=", "))
for (s in u) { k <- strsplit(s, ";")[[1]]; f <- paste0("TGT_", tk$id[match(s, psig)])
  A <- rbind(A[!rownames(A) %in% k, , drop=FALSE], colSums(A[k, , drop=FALSE])); rownames(A)[nrow(A)] <- f }
m <- A[match(feat, rownames(A)), , drop=FALSE]; prev <- rowMeans(m > 0)
cov <- data.frame(id=tk$id, target=tk$target, set=tk$set, n_keys=tk$n_keys, n_present=lengths(pres), feature=feat,
  reads=rowSums(m), prevalence=round(prev, 3), pass10=prev >= 0.10)
print(cov, row.names=FALSE)
if (dry) quit(save="no")
# --- fit ANCOM-BC2 (package defaults; Amendment 14b item 2). Features get safe IDs so no package renames them.
orig <- rownames(A); rownames(A) <- sprintf("F%05d", seq_along(orig)); sf <- rownames(A)[match(feat, orig)]
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(A)), colData=md[, terms, drop=FALSE])
fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH",
                prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
res <- fit$res; stopifnot(nrow(res) > 0, all(res$taxon %in% rownames(A)))
cat("ANCOM-BC2 kept", nrow(res), "of", nrow(A), "features\n"); i <- match(sf, res$taxon)
res$taxon_name <- orig[match(res$taxon, rownames(A))]; write.csv(res, file.path(out, paste0(coh, "_", rk, "_ancombc2.csv")), row.names=FALSE)
tc <- grep("^lfc_RFI", names(res), value=TRUE); stopifnot(length(tc) == 1); sfx <- sub("^lfc_", "", tc)
est <- data.frame(id=tk$id, set=tk$set, site=tk$site, cohort=coh, rank=rk, feature=feat, n=n, df=df,
  lfc=res[[tc]][i], se=res[[paste0("se_", sfx)]][i], p_ancom=res[[paste0("p_", sfx)]][i])
est$W <- est$lfc / est$se
est$status <- ifelse(is.na(feat), "no_feature", ifelse(is.na(i), "filtered", ifelse(is.finite(est$W), "ok", "no_estimate")))
write.table(est, file.path(out, paste0(coh, "_", rk, "_targets.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
print(est[, c("id", "feature", "status", "lfc", "se", "W", "p_ancom")], row.names=FALSE)
