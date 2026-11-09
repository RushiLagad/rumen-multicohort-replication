#!/usr/bin/env Rscript
# Amendment 38 item 3: de novo discovery scan, exploratory, from the real fits' full ANCOM-BC2 tables plus permutation nulls. No new fits. Usage: 91_denovo_scan.R
# Statistic throughout: W = lfc/se (the W_ column), as in scripts 21/65/81/83/89, referenced to the permutation null. ANCOM-BC2 2.4.0 p-values are not a normal
# transform of W and are reported for the nominal counts and the null false-positive rate only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "denovo_scan"); dir.create(od, showWarnings=FALSE)
stopifnot("commit Amendment 38 first" = any(grepl("^\\*\\*Amendment 38 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
stopifnot("scan nulls not finished" = file.exists(file.path(root, "logs", "scan_null", "DONE")))
rc <- function(f) read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
tabs <- list(
  nelore_genus=list(real="results/replication/nelore_genus_ancombc2.csv", null="results/scan_null", pat="^nelore_genus_seed2[0-9]{3}_ancombc2\\.csv$", trait="RFI"),
  charolais_mrna_genus=list(real="results/replication/charolais_mrna_genus_ancombc2.csv", null="results/scan_null", pat="^charolais_mrna_genus_seed2[0-9]{3}_ancombc2\\.csv$", trait="RFI"),
  alberta_genus=list(real="results/alberta/alberta_genus_ancombc2.csv", null="results/scan_null", pat="^alberta_genus_seed2[0-9]{3}_ancombc2\\.csv$", trait="RFI"),
  difford_bact=list(real="results/methane/difford_bact_ancombc2.csv", null="results/methane_dryrun", pat="^difford_bact_seed1[0-9]{3}_ancombc2\\.csv$", trait="CH4"),
  wallace_prokaryotes=list(real="results/methane/wallace_prokaryotes_ancombc2.csv", null="results/methane_dryrun", pat="^wallace_prokaryotes_seed1[0-9]{3}_ancombc2\\.csv$", trait="CH4"))
pull <- function(f, trait) { r <- rc(f); nm <- if ("taxon_name" %in% names(r)) r$taxon_name else r$taxon
  lc <- grep(if (trait == "RFI") "^lfc_RFI" else "^lfc_CH4_z", names(r), value=TRUE); stopifnot(length(lc) == 1); sfx <- sub("^lfc_", "", lc)
  d <- data.frame(genus=nm, lfc=r[[lc]], se=r[[paste0("se_", sfx)]], p=r[[paste0("p_", sfx)]], stringsAsFactors=FALSE); d$W <- d$lfc / d$se
  d[!grepl("^TGT_", d$genus) & is.finite(d$W), ] }
summ <- list(); scans <- list()
for (tn in names(tabs)) { t <- tabs[[tn]]; d <- pull(file.path(root, t$real), t$trait); d$q <- p.adjust(d$p, "BH")
  fs <- list.files(file.path(root, t$null), t$pat, full.names=TRUE); stopifnot(length(fs) >= 20)
  nul <- lapply(fs, function(f) pull(f, t$trait)); fpr <- sapply(nul, function(x) mean(x$p < 0.05)); nq <- sapply(nul, function(x) sum(p.adjust(x$p, "BH") < 0.05))
  NW <- do.call(rbind, lapply(nul, function(x) x[, c("genus", "W")])); Wn <- NW$W
  d$p_emp <- sapply(d$W, function(w) (sum(abs(Wn) >= abs(w)) + 1) / (length(Wn) + 1)); d$q_emp <- p.adjust(d$p_emp, "BH")   # pooled empirical null
  if (length(fs) >= 100) {   # per-genus calibration, as in 83/89, needs the 200-permutation tables
    sp <- split(NW$W, NW$genus); d$sd_null <- sapply(sp[d$genus], sd); d$k_null <- lengths(sp[d$genus]); d$W_cal <- d$W / pmax(1, d$sd_null)
    d$p_perm <- mapply(function(w, g) { v <- sp[[g]]; (sum(abs(v) >= abs(w)) + 1) / (length(v) + 1) }, d$W, d$genus); d$q_perm <- p.adjust(d$p_perm, "BH")
  } else { d$sd_null <- NA; d$k_null <- length(fs); d$W_cal <- NA; d$p_perm <- NA; d$q_perm <- NA }
  d$table <- tn; d$trait <- t$trait; d <- d[order(d$p_emp, -abs(d$W)), ]; scans[[tn]] <- d; write.table(d, file.path(od, paste0(tn, "_scan.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  summ[[tn]] <- data.frame(table=tn, trait=t$trait, n_genera=nrow(d), n_p05=sum(d$p < 0.05), n_q05=sum(d$q < 0.05), n_perm=length(fs),
    null_fpr_p05=round(mean(fpr), 3), null_max_n_q05=max(nq), null_W_sd=round(sd(Wn), 2), null_W_q99=round(quantile(abs(Wn), 0.99), 1), null_frac_W_gt5=signif(mean(abs(Wn) > 5), 2),
    n_q05_emp=sum(d$q_emp < 0.05), n_q05_perm=if (length(fs) >= 100) sum(d$q_perm < 0.05) else NA, top_genus=d$genus[1], top_W=round(d$W[1], 2), top_q_emp=signif(d$q_emp[1], 2)) }
S <- do.call(rbind, summ); write.table(S, file.path(od, "scan_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== DE NOVO SCAN (exploratory; statistic W = lfc/se; null = within-design permutations)\n"); print(S, row.names=FALSE)
cat("   n_p05/n_q05: ANCOM-BC2 p-values, nominal and BH across genera. null_fpr_p05: fraction of genera at p<0.05 in permuted data (0.05 if calibrated); null_max_n_q05: most q<0.05 genera in one permutation.\n",
    "  null_W_sd / q99 / frac>5: spread and tail of W under the null (1 / 2.6 / 1e-6 if standard normal). n_q05_emp: genera at BH q<0.05 against the pooled null W; n_q05_perm: same with per-genus null (200-permutation tables only).\n")
pairs <- list(c("nelore_genus", "charolais_mrna_genus"), c("difford_bact", "wallace_prokaryotes"), c("alberta_genus", "difford_bact"), c("alberta_genus", "wallace_prokaryotes"))
ov <- do.call(rbind, lapply(pairs, function(p) { a <- scans[[p[1]]]; b <- scans[[p[2]]]; g <- intersect(a$genus, b$genus); a <- a[match(g, a$genus), ]; b <- b[match(g, b$genus), ]
  ha <- a$p_emp < 0.05; hb <- b$p_emp < 0.05; both <- ha & hb; same <- both & sign(a$W) == sign(b$W); N <- length(g)
  exp_both <- sum(ha) * sum(hb) / N; p_hyper <- phyper(sum(both) - 1, sum(ha), N - sum(ha), sum(hb), lower.tail=FALSE)
  rho <- if (N >= 5) cor(a$W, b$W, method="spearman") else NA; rp <- if (N >= 5) cor.test(a$W, b$W, method="spearman", exact=FALSE)$p.value else NA
  data.frame(table_a=p[1], table_b=p[2], comparison=if (a$trait[1] == b$trait[1]) paste("same trait:", a$trait[1]) else "cross-trait (RFI vs CH4)", shared_genera=N, hits_a=sum(ha), hits_b=sum(hb),
    both=sum(both), both_same_sign=sum(same), expected_both=round(exp_both, 2), p_hyper=signif(p_hyper, 2), spearman_W=round(rho, 3), p_spearman=signif(rp, 2),
    same_sign_genera=paste(g[same], collapse=";"), stringsAsFactors=FALSE) }))
write.table(ov, file.path(od, "overlap.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== CROSS-COHORT OVERLAP of empirical-null hits (p_emp<0.05), shared genus names only; Spearman over all shared genera\n"); print(ov[, setdiff(names(ov), "same_sign_genera")], row.names=FALSE)
for (i in seq_len(nrow(ov))) if (nzchar(ov$same_sign_genera[i])) cat("  ", ov$table_a[i], "x", ov$table_b[i], "same-sign:", ov$same_sign_genera[i], "\n")
cat("== top 6 genera per table by empirical p\n"); for (tn in names(scans)) { d <- scans[[tn]]; cat("--", tn, "\n"); print(head(d[, c("genus", "lfc", "W", "p", "q", "p_emp", "q_emp", "sd_null", "q_perm")], 6), digits=3, row.names=FALSE) }
