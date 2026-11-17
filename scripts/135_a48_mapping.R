#!/usr/bin/env Rscript
# Amendment 48.8: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Mapping confidence (no refit).
# Classes per target and cohort, by rule from plan/target_keys.tsv and results/late_mapping/target_mapping_by_cohort.tsv (Amendment 45):
#   MC1 one key, the exact genus name in that cohort's taxonomy; MC2 GTDB-suffixed genus or several keys summed; MC3 family, phylum, clade or other aggregate; MC4 species.
# Shotgun cohorts use the GTDB keys; 16S cohorts the keys present among the cohort's SILVA genus names (results/16s/<cohort>/genus_lineage_silva.tsv), as 65 and 81 map them.
# Labels (47.2 rules, parametric p) and benchmark exclusions recomputed from the MC1 cohort estimates only, re-pooled with mKH where k >= 2 remains.
# Usage: 135_a48_mapping.R. Writes results/late_a48/mapping/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "mapping"); outs <- file.path(od, c("mapping_classes.tsv", "mc1_pooled.tsv", "summary.tsv")); a48_fresh(outs)
cat(A48_LABEL, "| 48.8 mapping confidence\n")
tk <- rd(file.path(root, "plan", "target_keys.tsv")); ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv"))
PI <- rd(file.path(root, "results", "late_a47", "intervals", "pooled_intervals.tsv")); PI <- PI[PI$version == "primary", ]; dirs <- a48_dirs(root)
mv <- rd(file.path(root, "results", "late_mapping", "charolais_mrna_vs_16s_targets.tsv"))
l47 <- rbind(rd(file.path(root, "results", "late_a47", "labels", "labels_rfi.tsv")), rd(file.path(root, "results", "late_a47", "labels", "labels_methane.tsv")))
lin <- list(alberta_all="alberta", difford="difford_bact", wallace="wallace_prokaryotes")
silva <- lapply(lin, function(c) { f <- file.path(root, "results", "16s", c, "genus_lineage_silva.tsv"); if (file.exists(f)) unique(rd(f)$genus) else NULL })
cls <- function(id, cohort) { x <- tk[tk$id == id, ]; rk <- x$target_rank; keys <- trimws(strsplit(ifelse(is.na(x$gtdb_keys), "", x$gtdb_keys), ";")[[1]]); keys <- keys[nzchar(keys)]
  if (rk == "species") return(c("MC4", "species level"))
  if (rk != "genus" || grepl("^unclassified", x$target)) return(c("MC3", paste0("aggregate (", rk, ")")))
  if (cohort %in% names(lin)) { pres <- intersect(keys, silva[[cohort]]); if (is.null(silva[[cohort]])) return(c(NA, "SILVA genus list missing"))
    if (length(pres) == 1 && pres == x$target) c("MC1", "one exact SILVA genus") else c("MC2", paste0(length(pres), " SILVA names: ", paste(pres, collapse=";"))) } else {
    if (length(keys) == 1 && keys == x$target) c("MC1", "one exact GTDB genus") else c("MC2", paste0(length(keys), " GTDB keys: ", paste(keys, collapse=";"))) } }
MC <- do.call(rbind, lapply(seq_len(nrow(PI)), function(i) { p <- PI[i, ]; do.call(rbind, lapply(strsplit(p$cohorts, ";")[[1]], function(co) { k <- cls(p$id, co)
  data.frame(id=p$id, target=p$target, arm=p$arm, cohort=co, rank=tk$target_rank[tk$id == p$id], class=k[1], basis=k[2],
    charolais_mrna_vs_16s_spearman=if (co == "charolais_mrna" && p$id %in% mv$id) mv$spearman_clr[mv$id == p$id][1] else NA) })) }))
# re-tabulation on MC1 cohort estimates only
key <- paste(MC$id, MC$cohort); est <- ce[paste(ce$id, ce$cohort) %in% key[MC$class %in% "MC1"], ]; est$W <- est$t
lab <- a48_labels(root, est); po <- a48_pool_targets(est, dirs)
lab1 <- lab[!duplicated(lab$id) & lab$label != "untestable", c("id", "label")]; po <- merge(po, lab1, by="id", all.x=TRUE)
full <- PI[, c("id", "k", "excl_15_mkh", "excl_30_mkh")]; names(full) <- c("id", "k_full", "excl_15_full", "excl_30_full"); po <- merge(full, po, by="id", all.x=TRUE)
po$k[is.na(po$k)] <- 0
SU <- data.frame(set=c("all targets (Amendment 47)", "MC1 cohort estimates only"),
  targets_testable=c(nrow(PI), sum(po$k > 0)), rfi_excl_30=c(sum(PI$arm == "RFI" & PI$excl_30_mkh), sum(po$id %in% PI$id[PI$arm == "RFI"] & po$excl_30 %in% TRUE)),
  methane_excl_15=c(sum(PI$arm == "Methane" & PI$excl_15_mkh), sum(po$id %in% PI$id[PI$arm == "Methane"] & po$excl_15 %in% TRUE)),
  replicated=c(sum(l47$label_a47 == "replicated"), sum(po$label %in% "replicated")), cohort_tests=c(nrow(MC), sum(MC$class %in% "MC1")),
  classes=c(paste(names(table(MC$class)), table(MC$class), sep=" ", collapse="; "), ""))
a48_write(MC, outs[1]); a48_write(po, outs[2]); a48_write(SU, outs[3])
cat("== classes\n"); print(table(MC$class, MC$cohort)); cat("\n== summary\n"); print(SU, row.names=FALSE)
