#!/usr/bin/env Rscript
# Amendment 42: pool CLR estimates as 89 (RFI, Alberta all-animal) and 83 (methane), compare with ANCOM-BC2. Usage: 96_clr_compare.R
suppressMessages(library(metafor)); root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")
sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote=""); cd <- file.path(root, "results", "clr"); od <- file.path(root, "results", "clr_pooled"); dir.create(od, showWarnings=FALSE)
stopifnot("commit Amendment 42 first" = any(grepl("^\\*\\*Amendment 42 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
tk0 <- rd(file.path(root, "plan", "target_keys.tsv")); cols <- c("id", "cohort", "W", "df", "n")
# ---- RFI
pR <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
tk <- tk0[!tk0$set %in% c("Methane", "MethaneSecondary"), ]
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available", ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
rf <- c(file.path(cd, paste0(c("nelore", "charolais_mrna"), "_", rep(c("species", "genus", "family", "phylum"), each=2), "_clr_targets.tsv")), file.path(cd, paste0("alberta_", c("genus", "family"), "_all_clr_targets.tsv"))); stopifnot(all(file.exists(rf)))
E <- do.call(rbind, lapply(rf, rd)); E <- E[E$status == "ok", ]; stopifnot(!anyDuplicated(E[, c("id", "cohort")]))
R <- do.call(rbind, lapply(unique(tk$set), function(st) { d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pR(d$direction), testable=d$stratum %in% c("replication", "cross-site"))
  cbind(set=st, target=d$target, classify_set(tg, E[E$id %in% d$id, cols])) }))
old <- rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")); j <- match(R$id, old$id)
write.table(R, file.path(od, "classification_rfi_clr.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
x <- data.frame(id=R$id, target=R$target, k=R$k, r_ancom=round(old$r[j], 3), r_clr=round(R$r, 3), ci_clr=paste0(round(R$ci_lo, 2), " to ", round(R$ci_hi, 2)), q_ancom=round(old$q_bh[j], 2), q_clr=round(R$q_bh, 2), label_ancom=old$label[j], label_clr=R$label)[!startsWith(R$label, "untestable") | !startsWith(old$label[j], "untestable"), ]
cat("== RFI, CLR vs calibrated ANCOM-BC2 (Alberta all-animal)\n"); print(x[order(x$q_clr), ], row.names=FALSE)
cat("label changes:", sum(x$label_ancom != x$label_clr, na.rm=TRUE), "| same sign of pooled r:", sum(sign(x$r_ancom) == sign(x$r_clr), na.rm=TRUE), "of", sum(is.finite(x$r_ancom) & is.finite(x$r_clr)), "| Spearman pooled r:", round(cor(x$r_ancom, x$r_clr, method="spearman", use="complete.obs"), 2), "\n")
# cohort-level W concordance, RFI
P <- do.call(rbind, lapply(c(file.path(root, "results", "replication", paste0(c("nelore", "charolais_mrna"), "_", rep(c("species", "genus", "family", "phylum"), each=2), "_targets.tsv")), file.path(root, "results", "alberta", paste0("alberta_", c("genus", "family"), "_all_targets.tsv"))), rd)); P <- P[P$status == "ok", ]
m <- merge(P[, c("id", "cohort", "W")], E[, c("id", "cohort", "W")], by=c("id", "cohort"), suffixes=c("_ancom", "_clr"))
cat("RFI target x cohort pairs:", nrow(m), "| same sign:", sum(sign(m$W_ancom) == sign(m$W_clr)), "| Spearman W:", round(cor(m$W_ancom, m$W_clr, method="spearman"), 2), "\n")
# ---- methane
pM <- function(x) ifelse(grepl("lower with higher", x), -1, ifelse(grepl("higher with higher", x), 1, NA))
tm <- tk0[tk0$set %in% c("Methane", "MethaneSecondary") | (tk0$set == "Nelore" & tk0$phenotype == "RME_GC"), ]
Em <- rbind(rd(file.path(cd, "difford_bact_clr_targets.tsv")), rd(file.path(cd, "wallace_prokaryotes_clr_targets.tsv"))); Em <- Em[Em$status == "ok", ]
m1 <- rbind(rd(file.path(root, "results", "methane", "difford_arch_targets.tsv")), rd(file.path(root, "results", "methane", "wallace_prokaryotes_targets.tsv"))); m1 <- m1[m1$id == "M1" & m1$status == "ok", ]   # M1 is a linear model in both versions
Em <- rbind(Em[, cols], m1[, cols]); Em$cohort <- ifelse(startsWith(Em$cohort, "difford"), "difford", "wallace"); Em <- Em[!(Em$id == "M2" & Em$cohort == "wallace"), ]   # Amendment 32
M <- do.call(rbind, lapply(unique(tm$set), function(st) { d <- tm[tm$set == st, ]; tg <- data.frame(id=d$id, dir=pM(d$direction), testable=d$status %in% c("mapped", "special"))
  cbind(set=st, target=d$target, classify_set(tg, Em[Em$id %in% d$id, ])) }))
oldm <- rd(file.path(root, "results", "methane_pooled", "classification_primary.tsv")); j <- match(M$id, oldm$id); write.table(M, file.path(od, "classification_methane_clr.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
y <- data.frame(id=M$id, target=M$target, k=M$k, r_ancom=round(oldm$r[j], 3), r_clr=round(M$r, 3), q_ancom=round(oldm$q_bh[j], 2), q_clr=round(M$q_bh, 2), label_ancom=oldm$label[j], label_clr=M$label)[!startsWith(M$label, "untestable") | !startsWith(oldm$label[j], "untestable"), ]
cat("\n== METHANE, CLR vs calibrated ANCOM-BC2\n"); print(y, row.names=FALSE)
cat("label changes:", sum(y$label_ancom != y$label_clr, na.rm=TRUE), "| same sign of pooled r:", sum(sign(y$r_ancom) == sign(y$r_clr), na.rm=TRUE), "of", sum(is.finite(y$r_ancom) & is.finite(y$r_clr)), "\n")
# ---- Wallace farm-robust genera (Amendment 40) under CLR, and transfer to Difford
f <- rd(file.path(root, "results", "wallace_farms", "farm_consistency.tsv")); g <- f$genus[f$farm_robust & f$genus != "TGT_M2"]
wa <- rd(file.path(root, "results", "wallace_farms", "full.tsv")); wc <- rd(file.path(cd, "wallace_prokaryotes_clr_full.tsv")); dc <- rd(file.path(cd, "difford_bact_clr_full.tsv"))
G <- data.frame(genus=g, W_ancom_wallace=round(wa$W[match(g, wa$genus)], 2), W_clr_wallace=round(wc$W[match(g, wc$genus)], 2), W_clr_difford=round(dc$W[match(g, dc$genus)], 2))
cat("\n== WALLACE FARM-ROBUST GENERA UNDER CLR\n"); print(G, row.names=FALSE)
cat("same sign ANCOM vs CLR in Wallace:", sum(sign(G$W_ancom_wallace) == sign(G$W_clr_wallace), na.rm=TRUE), "of", nrow(G), "| CLR |W| > 2 in Wallace:", sum(abs(G$W_clr_wallace) > 2, na.rm=TRUE),
  "| same sign Wallace vs Difford (CLR):", sum(sign(G$W_clr_wallace) == sign(G$W_clr_difford), na.rm=TRUE), "of", sum(is.finite(G$W_clr_difford)), "| CLR |W| > 2 in Difford:", sum(abs(G$W_clr_difford) > 2, na.rm=TRUE), "\n")
s <- intersect(wa$genus, wc$genus); s <- s[!grepl("^TGT_", s)]; cat("all Wallace genera, Spearman ANCOM vs CLR W:", round(cor(wa$W[match(s, wa$genus)], wc$W[match(s, wc$genus)], method="spearman"), 2), "\n")
s <- intersect(wc$genus, dc$genus); s <- s[!grepl("^TGT_", s)]; cat("Difford vs Wallace, CLR W over", length(s), "shared genera: Spearman", round(cor(wc$W[match(s, wc$genus)], dc$W[match(s, dc$genus)], method="spearman"), 2), "\n")
