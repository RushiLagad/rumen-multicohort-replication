#!/usr/bin/env Rscript
# Amendment 38 item 3: de novo discovery scan, exploratory, from the real fits' full ANCOM-BC2 tables plus permutation nulls. No new fits. Usage: 91_denovo_scan.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "denovo_scan"); dir.create(od, showWarnings=FALSE)
stopifnot("commit Amendment 38 first" = any(grepl("^\\*\\*Amendment 38 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
stopifnot("scan nulls not finished" = file.exists(file.path(root, "logs", "scan_null", "DONE")))
rc <- function(f) read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
tabs <- list(
  nelore_genus=list(real="results/replication/nelore_genus_ancombc2.csv", null="results/scan_null", pat="^nelore_genus_seed2[0-9]{3}_ancombc2\\.csv$", trait="RFI", ns="GTDB"),
  charolais_mrna_genus=list(real="results/replication/charolais_mrna_genus_ancombc2.csv", null="results/scan_null", pat="^charolais_mrna_genus_seed2[0-9]{3}_ancombc2\\.csv$", trait="RFI", ns="GTDB"),
  alberta_genus=list(real="results/alberta/alberta_genus_ancombc2.csv", null="results/scan_null", pat="^alberta_genus_seed2[0-9]{3}_ancombc2\\.csv$", trait="RFI", ns="SILVA"),
  difford_bact=list(real="results/methane/difford_bact_ancombc2.csv", null="results/methane_dryrun", pat="^difford_bact_seed1[0-9]{3}_ancombc2\\.csv$", trait="CH4", ns="SILVA"),
  wallace_prokaryotes=list(real="results/methane/wallace_prokaryotes_ancombc2.csv", null="results/methane_dryrun", pat="^wallace_prokaryotes_seed1[0-9]{3}_ancombc2\\.csv$", trait="CH4", ns="SILVA"))
# pull the phenotype coefficient for every genus (aggregates TGT_* and M1 excluded)
pull <- function(f, trait) { r <- rc(f); nm <- if ("taxon_name" %in% names(r)) r$taxon_name else r$taxon
  lc <- grep(if (trait == "RFI") "^lfc_RFI" else "^lfc_CH4_z", names(r), value=TRUE); stopifnot(length(lc) == 1); sfx <- sub("^lfc_", "", lc)
  d <- data.frame(genus=nm, lfc=r[[lc]], se=r[[paste0("se_", sfx)]], W=r[[paste0("W_", sfx)]], p=r[[paste0("p_", sfx)]], stringsAsFactors=FALSE)
  d <- d[!grepl("^TGT_", d$genus) & is.finite(d$W), ]; d }
summ <- list(); scans <- list()
for (tn in names(tabs)) { t <- tabs[[tn]]; d <- pull(file.path(root, t$real), t$trait); d$q <- p.adjust(d$p, "BH")
  fs <- list.files(file.path(root, t$null), t$pat, full.names=TRUE); stopifnot(length(fs) >= 20)
  nul <- lapply(fs, function(f) pull(f, t$trait)); fpr <- sapply(nul, function(x) mean(x$p < 0.05)); Wn <- unlist(lapply(nul, `[[`, "W"))
  lambda <- sd(Wn); d$W_cal <- d$W / max(1, lambda); d$p_cal <- 2 * pnorm(-abs(d$W_cal)); d$q_cal <- p.adjust(d$p_cal, "BH")
  d$p_emp <- sapply(d$W, function(w) (sum(abs(Wn) >= abs(w)) + 1) / (length(Wn) + 1)); d$q_emp <- p.adjust(d$p_emp, "BH")
  d$table <- tn; d$trait <- t$trait; d <- d[order(d$p), ]; scans[[tn]] <- d; write.table(d, file.path(od, paste0(tn, "_scan.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  summ[[tn]] <- data.frame(table=tn, trait=t$trait, n_genera=nrow(d), n_p05=sum(d$p < 0.05), n_q05=sum(d$q < 0.05), min_q=signif(min(d$q), 2),
    n_perm=length(fs), null_fpr_p05=round(mean(fpr), 3), null_fpr_sd=round(sd(fpr), 3), null_W_sd=round(lambda, 2), null_max_n_q05=max(sapply(nul, function(x) sum(p.adjust(x$p, "BH") < 0.05))),
    n_q05_cal=sum(d$q_cal < 0.05), n_q05_emp=sum(d$q_emp < 0.05), top_genus=d$genus[1], top_W=round(d$W[1], 2), top_q_cal=signif(d$q_cal[1], 2)) }
S <- do.call(rbind, summ); write.table(S, file.path(od, "scan_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== DE NOVO SCAN (exploratory; BH across genera within table; null = within-design permutations)\n"); print(S, row.names=FALSE)
cat("   null_fpr_p05: fraction of genera at nominal p<0.05 in permuted data (0.05 if calibrated). null_W_sd: SD of W under the null (1 if calibrated).\n   null_max_n_q05: largest number of q<0.05 genera seen in any single permutation. n_q05_cal / n_q05_emp: real counts after table-level calibration / against the pooled empirical null.\n")
# cross-cohort overlap within namespace; same-trait pairs first, cross-trait labelled
pairs <- list(c("nelore_genus", "charolais_mrna_genus"), c("difford_bact", "wallace_prokaryotes"), c("alberta_genus", "difford_bact"), c("alberta_genus", "wallace_prokaryotes"))
ov <- do.call(rbind, lapply(pairs, function(p) { a <- scans[[p[1]]]; b <- scans[[p[2]]]; g <- intersect(a$genus, b$genus); a <- a[match(g, a$genus), ]; b <- b[match(g, b$genus), ]
  ha <- a$p < 0.05; hb <- b$p < 0.05; both <- ha & hb; same <- both & sign(a$W) == sign(b$W); N <- length(g)
  exp_both <- sum(ha) * sum(hb) / N; p_hyper <- phyper(sum(both) - 1, sum(ha), N - sum(ha), sum(hb), lower.tail=FALSE)
  rho <- if (N >= 5) cor(a$W, b$W, method="spearman") else NA; rp <- if (N >= 5) cor.test(a$W, b$W, method="spearman", exact=FALSE)$p.value else NA
  data.frame(table_a=p[1], table_b=p[2], comparison=if (a$trait[1] == b$trait[1]) paste("same trait:", a$trait[1]) else "cross-trait (RFI vs CH4)", shared_genera=N, hits_a=sum(ha), hits_b=sum(hb),
    both=sum(both), both_same_sign=sum(same), expected_both=round(exp_both, 2), p_hyper=signif(p_hyper, 2), spearman_W=round(rho, 3), p_spearman=signif(rp, 2),
    same_sign_genera=paste(g[same], collapse=";"), stringsAsFactors=FALSE) }))
write.table(ov, file.path(od, "overlap.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== CROSS-COHORT OVERLAP of nominal hits (p<0.05), shared genus names only\n"); print(ov[, setdiff(names(ov), "same_sign_genera")], row.names=FALSE)
for (i in seq_len(nrow(ov))) if (nzchar(ov$same_sign_genera[i])) cat("  ", ov$table_a[i], "x", ov$table_b[i], "same-sign overlap:", ov$same_sign_genera[i], "\n")
cat("== top 5 genera per table (real, calibrated q)\n"); for (tn in names(scans)) { d <- scans[[tn]]; cat("--", tn, "\n"); print(head(d[, c("genus", "lfc", "W", "p", "q", "q_cal", "q_emp")], 5), digits=3, row.names=FALSE) }
