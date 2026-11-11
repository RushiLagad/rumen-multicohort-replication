#!/usr/bin/env Rscript
# Amendments 19, 28, 32, 34, 35, 36: methane target estimates for one 16S cohort table.
# Usage: 81_estimate_methane.R <difford_bact|difford_arch|wallace_prokaryotes|wallace_bacteria|wallace_archaea> <shuffle|real>
# Env: NO_SE1=1 (Wallace without farm SE1, Amendment 32 item 2); TRAIT=ch4_dmi (Wallace sensitivity trait, Amendment 34). No phenotype value is printed.
a <- commandArgs(TRUE); coh <- a[1]; mode <- a[2]; stopifnot(length(a) == 2, mode %in% c("shuffle", "real"))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); SEED <- as.integer(Sys.getenv("SEED", "20261105")); set.seed(SEED)
plan <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")); stopifnot("commit Amendment 36 first" = any(grepl("^\\*\\*Amendment 36 \\(", plan)))
if (mode == "real") stopifnot("real mode is locked until the dry-run result is recorded in the plan (Amendment 36 item 6)" = any(grepl("^\\*\\*Amendment 36 dry-run result", plan)))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
src <- if (startsWith(coh, "difford")) "difford" else "wallace"; trait <- Sys.getenv("TRAIT", "ch4_gd"); noSE1 <- Sys.getenv("NO_SE1") == "1"
if (src == "difford") stopifnot(trait == "ch4_gd", !noSE1)
tag <- paste0(if (trait != "ch4_gd") paste0("_", trait) else "", if (noSE1) "_noSE1" else "", if (Sys.getenv("DEPTH") == "1") "_depth" else "", if (mode == "shuffle" && SEED != 20261105) paste0("_seed", SEED) else "")
out <- file.path(root, "results", if (mode == "shuffle") "methane_dryrun" else "methane"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
X <- rd(file.path(root, "profiles_16s", coh, paste0(coh, "_genus_counts.tsv"))); rownames(X) <- X[[1]]; X <- as.matrix(X[, -1, drop=FALSE])
lin <- rd(file.path(root, "results", "16s", coh, "genus_lineage_silva.tsv"))
an <- rd(file.path(root, "harmonized", paste0(src, "_animals.tsv"))); an$cow <- as.character(an$cow)
if (src == "wallace") { an <- an[an$sequenced %in% c("True", TRUE, "TRUE"), ]; if (noSE1) an <- an[an$farm != "SE1", ] }
an$y <- suppressWarnings(as.numeric(an[[trait]])); an$dim <- suppressWarnings(as.numeric(an$dim))
an <- an[!is.na(an$y) & !is.na(an$dim) & an$parity %in% c("1", "2", "3+") & an$cow %in% colnames(X), ]
grp <- if (src == "difford") an$herd else an$farm
if (mode == "shuffle") { for (g in unique(grp)) { i <- which(grp == g); an$y[i] <- sample(an$y[i]) }; cat("DRY RUN: methane permuted within", if (src == "difford") "herd" else "farm", "\n") }
if (src == "difford") an$CH4_z <- as.numeric(scale(an$y)) else { an$CH4_z <- NA_real_; for (g in unique(grp)) { i <- which(grp == g); an$CH4_z[i] <- as.numeric(scale(an$y[i])) } }
md <- if (src == "difford") data.frame(herd=factor(an$herd), parity=factor(an$parity), dim=an$dim, batch=factor(an$batch), CH4_z=an$CH4_z, row.names=an$cow) else
  data.frame(farm=factor(an$farm), parity=factor(an$parity), dim=an$dim, CH4_z=an$CH4_z, row.names=an$cow)
if (Sys.getenv("DEPTH") == "1") { stopifnot("DEPTH=1 is for Wallace only (Amendment 41)" = src == "wallace"); md$log_depth <- log10(colSums(X[, rownames(md), drop=FALSE])); md <- md[, c("farm", "parity", "dim", "log_depth", "CH4_z")] }
terms <- c(setdiff(names(md), "CH4_z"), "CH4_z"); mm <- model.matrix(reformulate(terms), md); small <- character(0)
if (qr(mm)$rank < ncol(mm)) { stopifnot(src == "difford"); cat("herd + batch not of full rank: management group = herd x batch (Amendments 23, 28)\n")
  md$group <- factor(paste(an$herd, an$batch, sep="|")); gs <- table(md$group); small <- names(gs)[gs < 5]; keep <- !md$group %in% small
  md <- md[keep, ]; an <- an[keep, ]; md$group <- droplevels(md$group); md$herd <- NULL; md$batch <- NULL
  terms <- c("group", "parity", "dim", "CH4_z"); mm <- model.matrix(reformulate(terms), md); cat("groups", nlevels(md$group), "| excluded", length(small), "small groups\n") }
stopifnot("design not of full rank" = qr(mm)$rank == ncol(mm))
X <- X[, rownames(md), drop=FALSE]; n <- ncol(X); df <- n - ncol(mm)
cat(coh, mode, if (nzchar(tag)) tag else "", "| cows", n, "| df", df, "| terms:", paste(terms, collapse=" + "), "\n")
tk <- rd(file.path(root, "plan", "target_keys.tsv")); tk <- if (trait == "rfi_nrc") tk[grepl("^RFI", tk$phenotype) & tk$status == "mapped" & tk$n_keys > 0 & tk$target_rank %in% c("genus", "family"), ] else tk[(tk$set %in% c("Methane", "MethaneSecondary") | (tk$set == "Nelore" & tk$phenotype == "RME_GC")) & tk$status %in% c("mapped", "special"), ]   # Amendment 43
fam <- lin$family[match(rownames(X), lin$genus)]; fam[is.na(fam)] <- "unassigned_family"; kng <- lin$kingdom[match(rownames(X), lin$genus)]
A <- X; feat <- rep(NA_character_, nrow(tk)); npres <- integer(nrow(tk))
for (j in seq_len(nrow(tk))) { rk <- tk$target_rank[j]; if (rk == "clade_ratio") next
  keys <- trimws(strsplit(tk$gtdb_keys[j], ";")[[1]]); keys <- keys[nzchar(keys)]
  rows <- if (rk == "genus") intersect(keys, rownames(X)) else if (rk == "family") rownames(X)[fam %in% keys] else if (rk == "kingdom") rownames(X)[kng %in% keys] else character(0)
  npres[j] <- length(rows); if (!length(rows)) next
  if (rk == "genus" && length(rows) == 1) { feat[j] <- rows } else { f <- paste0("TGT_", tk$id[j]); A <- rbind(A, colSums(X[rows, , drop=FALSE])); rownames(A)[nrow(A)] <- f; feat[j] <- f } }   # aggregate targets are added alongside the genera (Amendment 36 item 3)
m <- A[match(feat, rownames(A)), , drop=FALSE]; prev <- rowMeans(m > 0, na.rm=TRUE)
cov <- data.frame(id=tk$id, set=tk$set, target=tk$target, rank=tk$target_rank, n_present=npres, feature=feat, prevalence=round(unname(prev), 3), pass10=unname(prev) >= 0.10)
cat("targets", nrow(tk), "| matched", sum(!is.na(feat)), "| pass 10% prevalence", sum(cov$pass10, na.rm=TRUE), "| features", nrow(A), "\n"); print(cov[cov$set != "Nelore" | !is.na(cov$feature), ], row.names=FALSE)
if (Sys.getenv("CLR") == "1") { source(file.path(root, "scripts", "lib_clr.R")); clr_write(A, feat, tk, md, terms, n, df, cohort=paste0(coh, tag), rank=tk$target_rank, prefix=paste0(coh, tag), root=root, m1=TRUE, difford=startsWith(coh, "difford")); quit(save="no") }
if (nzchar(Sys.getenv("LATE"))) { stopifnot("late analyses use the primary methane fit (real mode, CH4 g/d, all farms)" = mode == "real", !nzchar(tag)); source(file.path(root, "scripts", "lib_late.R"))   # Amendment 44
  late_main(A, feat, tk, md, terms, n, df, cohort=coh, rank=tk$target_rank, prefix=coh, root=root, perm_by=if (src == "difford") "herd" else "farm", m1=TRUE, difford=startsWith(coh, "difford")); quit(save="no") }
suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment)})
orig <- rownames(A); rownames(A) <- sprintf("F%05d", seq_along(orig)); sf <- rownames(A)[match(feat, orig)]
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(A)), colData=md[, terms, drop=FALSE])
fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH", prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
res <- fit$res; stopifnot(nrow(res) > 0); cat("ANCOM-BC2 kept", nrow(res), "of", nrow(A), "features\n"); i <- match(sf, res$taxon)
res$taxon_name <- orig[match(res$taxon, rownames(A))]; write.csv(res, file.path(out, paste0(coh, tag, "_ancombc2.csv")), row.names=FALSE)
tc <- grep("^lfc_CH4_z", names(res), value=TRUE); stopifnot(length(tc) == 1); sfx <- sub("^lfc_", "", tc)
est <- data.frame(id=tk$id, set=tk$set, site=tk$site, cohort=paste0(coh, tag), rank=tk$target_rank, feature=feat, n=n, df=df,
  lfc=res[[tc]][i], se=res[[paste0("se_", sfx)]][i], p_ancom=res[[paste0("p_", sfx)]][i]); est$W <- est$lfc / est$se
est$status <- ifelse(is.na(feat), "no_feature", ifelse(is.na(i), "filtered", ifelse(is.finite(est$W), "ok", "no_estimate")))
if (coh %in% c("difford_arch", "wallace_prokaryotes")) { m1 <- rd(file.path(root, "profiles_16s", coh, paste0(coh, "_m1_counts.tsv"))); m1$cow <- as.character(m1$cow)
  m1 <- m1[match(rownames(md), m1$cow), ]; ok <- !is.na(m1$cow) & (m1$gottschalkii + m1$ruminantium) >= 100
  d1 <- md[ok, ]; d1$lr <- log((m1$gottschalkii[ok] + 0.5) / (m1$ruminantium[ok] + 0.5)); fm <- lm(reformulate(terms, "lr"), d1); s <- summary(fm)$coefficients["CH4_z", ]
  k <- which(est$id == "M1"); est$feature[k] <- "M1_log_ratio"; est$n[k] <- sum(ok); est$df[k] <- fm$df.residual; est$lfc[k] <- s[1]; est$se[k] <- s[2]; est$p_ancom[k] <- s[4]; est$W[k] <- s[3]; est$status[k] <- "ok"
  cat("M1: cows with >= 100 clade-assigned reads", sum(ok), "of", nrow(md), "\n") } else est$status[est$id == "M1"] <- "untestable_here"
if (startsWith(coh, "difford")) est$status[est$id == "M3"] <- "untestable_here"
write.table(est, file.path(out, paste0(coh, tag, "_targets.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
print(est[est$set != "Nelore" | est$status == "ok", c("id", "set", "feature", "status", "n", "lfc", "se", "W", "p_ancom")], row.names=FALSE, digits=3)
