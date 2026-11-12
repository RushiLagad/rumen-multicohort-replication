#!/usr/bin/env Rscript
# Amendment 45 item 3: target-mapping concordance. Per target and cohort: matched feature, prevalence and fit status (from the primary target files);
# for the 30 Charolais animals with both mRNA (GTDB) and 16S (SILVA) data, Spearman correlation of the target's CLR abundance across animals and the
# diet W in each layer. DNA vs RNA is confounded with GTDB vs SILVA here, so agreement is a lower bound. Late analysis. Usage: 109_late_mapping.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "late_mapping"); dir.create(od, showWarnings=FALSE)
stopifnot("commit Amendment 45 first" = any(grepl("^\\*\\*Amendment 45 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
LATE_LABEL <- "late analysis, run after all primary results were seen"; rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
src <- c(nelore="results/replication/nelore_%s_targets.tsv", charolais_mrna="results/replication/charolais_mrna_%s_targets.tsv", alberta="results/alberta/alberta_%s_all_targets.tsv")
T <- data.frame(id=tk$id, set=tk$set, target=tk$target, rank=tk$target_rank, keys=tk$gtdb_keys, key_status=tk$status)
for (c in names(src)) { fs <- sprintf(file.path(root, src[[c]]), c("species", "genus", "family", "phylum")); fs <- fs[file.exists(fs)]
  e <- do.call(rbind, lapply(fs, rd)); e <- e[!duplicated(e$id), ]; j <- match(T$id, e$id); T[[paste0(c, "_feature")]] <- e$feature[j]; T[[paste0(c, "_status")]] <- e$status[j] }
for (c in c("difford_bact", "wallace_prokaryotes")) { e <- rd(file.path(root, "results", "methane", paste0(c, "_targets.tsv"))); j <- match(T$id, e$id); T[[paste0(c, "_feature")]] <- e$feature[j]; T[[paste0(c, "_status")]] <- e$status[j] }
# same-animal Charolais comparison
S <- list()
for (rk in c("genus", "family", "phylum")) { fm <- file.path(od, paste0("charolais_mrna_", rk, "_target_clr.tsv")); fs <- file.path(od, paste0("charolais_16s_", rk, "_target_clr.tsv"))
  if (!file.exists(fm) || !file.exists(fs)) { cat("missing", rk, "\n"); next }
  a <- rd(fm); b <- rd(fs); sm <- rd(sub("_clr", "_status", fm)); s6 <- rd(sub("_clr", "_status", fs))
  for (i in union(sm$id, s6$id)) { x <- a[a$id == i, ]; y <- b[b$id == i, ]; an <- intersect(x$animal, y$animal)
    r <- if (length(an) >= 10) suppressWarnings(cor.test(x$clr[match(an, x$animal)], y$clr[match(an, y$animal)], method="spearman")) else NULL
    S[[paste(rk, i)]] <- data.frame(id=i, rank=rk, mrna_feature=sm$feature[match(i, sm$id)], mrna_status=sm$status[match(i, sm$id)], mrna_prevalence=sm$prevalence[match(i, sm$id)],
      s16_feature=s6$feature[match(i, s6$id)], s16_status=s6$status[match(i, s6$id)], s16_prevalence=s6$prevalence[match(i, s6$id)], animals=length(an),
      spearman_clr=if (is.null(r)) NA else round(unname(r$estimate), 3), p=if (is.null(r)) NA else signif(r$p.value, 3)) } }
S <- do.call(rbind, S)
dw <- function(f, ids) { x <- read.csv(file.path(root, f), check.names=FALSE, stringsAsFactors=FALSE); w <- x$lfc_dietGS / x$se_dietGS; w[match(paste0("TGT_", ids), x$taxon_name)] }
if (!is.null(S)) { S$diet_W_mrna <- NA; S$diet_W_16s <- NA
  for (rk in unique(S$rank)) { k <- S$rank == rk
    fm <- sprintf("results/replication/charolais_mrna_%s_ancombc2.csv", rk); f6 <- sprintf("results/charolais16s/charolais_16s_%s_ancombc2.csv", rk)
    if (file.exists(file.path(root, fm))) S$diet_W_mrna[k] <- round(dw(fm, S$id[k]), 2); if (file.exists(file.path(root, f6))) S$diet_W_16s[k] <- round(dw(f6, S$id[k]), 2) } }
T$charolais_same_animal_spearman <- S$spearman_clr[match(T$id, S$id)]; T$analysis <- LATE_LABEL
write.table(T, file.path(od, "target_mapping_by_cohort.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (!is.null(S)) { S$analysis <- LATE_LABEL; write.table(S, file.path(od, "charolais_mrna_vs_16s_targets.tsv"), sep="\t", quote=FALSE, row.names=FALSE) }
cat("== Amendment 45 item 3: target mapping (", LATE_LABEL, ")\n", sep="")
for (c in c("nelore", "charolais_mrna", "alberta", "difford_bact", "wallace_prokaryotes")) cat(sprintf("%-20s targets with a matched feature %3d | fitted (status ok) %3d\n", c, sum(!is.na(T[[paste0(c, "_feature")]])), sum(T[[paste0(c, "_status")]] %in% "ok")))
if (!is.null(S)) { ok <- is.finite(S$spearman_clr); cat("Charolais same-animal, targets above 10% in both layers:", sum(ok), "| median Spearman of CLR across animals", round(median(S$spearman_clr[ok]), 2),
  "| range", round(min(S$spearman_clr[ok]), 2), "to", round(max(S$spearman_clr[ok]), 2), "| p < 0.05 and positive:", sum(S$spearman_clr[ok] > 0 & S$p[ok] < 0.05), "\n")
  print(S[order(-S$spearman_clr), c("id", "rank", "mrna_feature", "s16_feature", "animals", "spearman_clr", "p", "diet_W_mrna", "diet_W_16s")], row.names=FALSE) }
