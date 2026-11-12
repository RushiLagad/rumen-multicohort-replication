#!/usr/bin/env Rscript
# Amendment 45 item 1: positive controls for transfer, from fits already on file (no new model fits). Late analysis, run after all primary results were seen.
# (a) Across cohorts, same platform: Difford vs Wallace genus effects of days in milk and parity (covariates of the primary methane fits), beside the methane term.
# (b) Same animals, different platform and taxonomy: Charolais mRNA (GTDB) vs 16S (SILVA) genus effects of diet and year, beside the RFI class term.
# Raw ANCOM-BC2 W = lfc/se throughout (covariate terms were not permuted). Usage: 108_late_positive_controls.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "late_positive_controls"); dir.create(od, showWarnings=FALSE)
stopifnot("commit Amendment 45 first" = any(grepl("^\\*\\*Amendment 45 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
LATE_LABEL <- "late analysis, run after all primary results were seen"
rdc <- function(f) { x <- read.csv(file.path(root, f), check.names=FALSE, stringsAsFactors=FALSE); x <- x[!grepl("^TGT_|^unassigned", x$taxon_name), ]; x }
W <- function(x, term) { l <- x[[paste0("lfc_", term)]]; s <- x[[paste0("se_", term)]]; stopifnot("term not in table" = !is.null(l)); setNames(l / s, x$taxon_name) }
conc <- function(a, b, ta, tb, comparison, label) {
  wa <- W(a, ta); wb <- W(b, tb); g <- intersect(names(wa), names(wb)); g <- g[is.finite(wa[g]) & is.finite(wb[g])]; x <- wa[g]; y <- wb[g]; N <- length(g)
  ct <- suppressWarnings(cor.test(x, y, method="spearman")); ha <- abs(x) > 1.96; hb <- abs(y) > 1.96; both <- ha & hb; same <- sum(both & sign(x) == sign(y))
  expd <- sum(ha) * sum(hb) / N; p_ov <- phyper(sum(both) - 1, sum(ha), N - sum(ha), sum(hb), lower.tail=FALSE)
  p_sign <- if (sum(both)) binom.test(same, sum(both), 0.5, alternative="greater")$p.value else NA
  pass <- ct$estimate > 0 && ct$p.value < 0.05 && same > expd / 2
  data.frame(comparison=comparison, term=label, n_shared_genera=N, spearman_W=round(unname(ct$estimate), 3), p_spearman=signif(ct$p.value, 3),
    nominal_hits_1=sum(ha), nominal_hits_2=sum(hb), hits_in_both=sum(both), expected_in_both=round(expd, 2), p_overlap=signif(p_ov, 3),
    same_sign_in_both=same, expected_same_sign=round(expd / 2, 2), p_sign=signif(p_sign, 3), pass=pass, analysis=LATE_LABEL) }
dif <- rdc("results/methane/difford_bact_ancombc2.csv"); wal <- rdc("results/methane/wallace_prokaryotes_ancombc2.csv")
cm <- rdc("results/replication/charolais_mrna_genus_ancombc2.csv"); c16 <- rdc("results/charolais16s/charolais_16s_genus_ancombc2.csv")
R <- rbind(
  conc(dif, wal, "dim", "dim", "Difford vs Wallace (16S, SILVA; different herds)", "days in milk"),
  conc(dif, wal, "parity2", "parity2", "Difford vs Wallace (16S, SILVA; different herds)", "parity 2 vs 1"),
  conc(dif, wal, "parity3+", "parity3+", "Difford vs Wallace (16S, SILVA; different herds)", "parity 3+ vs 1"),
  conc(dif, wal, "CH4_z", "CH4_z", "Difford vs Wallace (16S, SILVA; different herds)", "methane (reference: the tested trait)"),
  conc(cm, c16, "dietGS", "dietGS", "Charolais mRNA (GTDB) vs 16S (SILVA), same 30 animals", "diet, grass vs corn silage"),
  conc(cm, c16, "year2020", "year2020", "Charolais mRNA (GTDB) vs 16S (SILVA), same 30 animals", "year 2020 vs 2019"),
  conc(cm, c16, "RFI_classhigh_RFI", "RFI_classhigh_RFI", "Charolais mRNA (GTDB) vs 16S (SILVA), same 30 animals", "RFI class (reference: the tested trait)"))
write.table(R, file.path(od, "transfer_positive_controls.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
for (p in list(list(dif, wal, "difford_vs_wallace", c("dim", "parity2", "parity3+", "CH4_z")), list(cm, c16, "charolais_mrna_vs_16s", c("dietGS", "year2020", "RFI_classhigh_RFI")))) {
  g <- intersect(p[[1]]$taxon_name, p[[2]]$taxon_name); M <- data.frame(genus=g)
  for (t in p[[4]]) { M[[paste0("W1_", t)]] <- W(p[[1]], t)[g]; M[[paste0("W2_", t)]] <- W(p[[2]], t)[g] }
  write.table(M, file.path(od, paste0(p[[3]], "_genus_W.tsv")), sep="\t", quote=FALSE, row.names=FALSE) }
cat("== Amendment 45 item 1: positive controls for transfer (", LATE_LABEL, "; raw ANCOM-BC2 W)\n", sep=""); print(R[, -ncol(R)], row.names=FALSE)
