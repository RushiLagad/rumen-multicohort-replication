#!/usr/bin/env Rscript
# Amendment 48.2: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. ANCOM-BC2 configuration and filtered calibration.
# (a) audit of every ancombc2() call (code text, with unset arguments resolved from formals() of the installed ANCOMBC) and a check that the labels use
#     parametric p of calibrated W (k = 1 rows of the Amendment 37 classifications against 2 pt(-|W|, df)); (b) the 200 saved null fits per table (Fig. 3a files
#     and genus filter) and the real and positive-control fits, recounted with the sensitivity filter (passed_ss); (c) Charolais genus refit with
#     group = RFI_class, struc_zero = TRUE, after reproducing the primary fit. Reading (d) as 48.2. Usage: 129_a48_ancombc2.R. Writes results/late_a48/ancombc2/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "ancombc2"); outs <- file.path(od, c("config.tsv", "label_p_check.tsv", "null_per_shuffle.tsv", "null_summary.tsv", "real_fits.tsv", "positive_controls.tsv",
  "charolais_struc_zero.tsv", "charolais_struc_zero_taxa.tsv", "reading.tsv")); a48_fresh(outs)
NSEED <- as.integer(Sys.getenv("A48_TEST_SEEDS", "200")); if (NSEED != 200) cat("NOTE: A48_TEST_SEEDS =", NSEED, "(test run; the plan uses all 200)\n")
cat(A48_LABEL, "| 48.2 ANCOM-BC2 configuration and filtered calibration\n")
rc <- function(f) read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)

# ---------------------------------------------------------------- (a) audit
has_pkg <- requireNamespace("ANCOMBC", quietly=TRUE)
fm <- if (has_pkg) formals(ANCOMBC::ancombc2) else NULL; ver <- if (has_pkg) as.character(packageVersion("ANCOMBC")) else "not installed in this R"
calls <- do.call(rbind, lapply(c("20_positive_controls.R", "21_estimate_targets.R", "65_estimate_alberta.R", "81_estimate_methane.R", "90_positive_controls_16s.R"), function(s) {
  L <- readLines(file.path(root, "scripts", s)); i <- grep("ancombc2\\(", L)[1]; txt <- L[i]; j <- i
  while (sum(gregexpr("\\(", txt)[[1]] > 0) > sum(gregexpr("\\)", txt)[[1]] > 0) && j < length(L)) { j <- j + 1; txt <- paste(txt, trimws(L[j])) }
  txt <- sub(".*ancombc2\\(", "", txt); args <- trimws(strsplit(sub("\\)[^)]*$", "", txt), ",(?![^(]*\\))", perl=TRUE)[[1]]); kv <- do.call(rbind, strsplit(sub("=", "\u0001", args), "\u0001"))
  set <- setNames(trimws(kv[, 2]), trimws(kv[, 1])); nm <- union(names(set), names(fm))
  data.frame(script=s, line=i, argument=nm, value=ifelse(nm %in% names(set), set[nm], "(default)"),
    resolved=sapply(nm, function(a) if (a %in% names(set)) set[[a]] else if (!is.null(fm) && a %in% names(fm)) paste(deparse(fm[[a]]), collapse="") else NA),
    ancombc_version=ver) }))
# labels use parametric p of calibrated W: for k = 1 targets the classification p must equal 2 pt(-|W|, df) of that cohort
er <- rd(file.path(root, "results", "replication_rfi3_calibrated", "estimates_calibrated_all.tsv")); em <- rd(file.path(root, "results", "methane_pooled", "estimates_calibrated.tsv"))
c37r <- rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")); c37m <- rd(file.path(root, "results", "methane_pooled", "classification_primary.tsv"))
lp <- function(cl, es, arm) { x <- cl[cl$k == 1, ]; do.call(rbind, lapply(seq_len(nrow(x)), function(i) { e <- es[es$id == x$id[i] & is.finite(es$W), ]
  e <- e[if (nrow(e) > 1) which.min(abs(e$r - x$r[i])) else seq_len(nrow(e)), ]; if (!nrow(e)) return(NULL)
  data.frame(arm=arm, set=x$set[i], id=x$id[i], cohort=e$cohort, p_classification=x$p[i], p_parametric_calW=2 * pt(-abs(e$W), e$df), p_perm_200=e$p_perm) })) }
LP <- rbind(lp(c37r, er, "RFI"), lp(c37m, em, "Methane")); LP$abs_diff <- abs(LP$p_classification - LP$p_parametric_calW)
LP$verdict <- ifelse(LP$abs_diff < 1e-8, "parametric p of calibrated W", "differs")
cat("== (a) labels: k = 1 classification p equals 2 pt(-|W_cal|, df) in", sum(LP$verdict == "parametric p of calibrated W"), "of", nrow(LP), "rows (max |diff|", signif(max(LP$abs_diff), 3), ")\n")

# ---------------------------------------------------------------- (b) per fit: shares and discoveries with and without the filter
fitstats <- function(f, trait) { x <- rc(f); nm <- if ("taxon_name" %in% names(x)) x$taxon_name else x$taxon
  pc <- grep(paste0("^p_(", trait, ")"), names(x), value=TRUE); stopifnot(length(pc) == 1); t <- sub("^p_", "", pc)
  need <- paste0(c("lfc_", "se_", "q_", "diff_", "passed_ss_"), t); if (!all(need %in% names(x))) return(NULL)
  keep <- !startsWith(as.character(nm), "TGT_") & is.finite(x[[paste0("lfc_", t)]] / x[[paste0("se_", t)]])
  p <- x[[pc]][keep]; q <- x[[paste0("q_", t)]][keep]; d <- as.logical(x[[paste0("diff_", t)]][keep]); ss <- as.logical(x[[paste0("passed_ss_", t)]][keep])
  data.frame(n_genera=sum(keep), share_p05=mean(p < 0.05), share_p05_passed_all=mean(p < 0.05 & ss), share_p05_among_passed=mean(p[ss] < 0.05), share_passed=mean(ss),
    n_bh=sum(d), n_robust=sum(d & ss), any_bh=any(d), any_robust=any(d & ss), diff_vs_q_mismatch=sum(d != (q < 0.05), na.rm=TRUE)) }
TABS <- list(nelore_genus=list(null="results/scan_null", pat="^nelore_genus_seed2[0-9]{3}_ancombc2\\.csv$", s0=2001, real="results/replication/nelore_genus_ancombc2.csv", trait="RFI"),
  charolais_mrna_genus=list(null="results/scan_null", pat="^charolais_mrna_genus_seed2[0-9]{3}_ancombc2\\.csv$", s0=2001, real="results/replication/charolais_mrna_genus_ancombc2.csv", trait="RFI"),
  alberta_genus=list(null="results/scan_null", pat="^alberta_genus_seed2[0-9]{3}_ancombc2\\.csv$", s0=2001, real="results/alberta/alberta_genus_ancombc2.csv", trait="RFI"),
  difford_bact=list(null="results/methane_dryrun", pat="^difford_bact_seed1[0-9]{3}_ancombc2\\.csv$", s0=1001, real="results/methane/difford_bact_ancombc2.csv", trait="CH4"),
  wallace_prokaryotes=list(null="results/methane_dryrun", pat="^wallace_prokaryotes_seed1[0-9]{3}_ancombc2\\.csv$", s0=1001, real="results/methane/wallace_prokaryotes_ancombc2.csv", trait="CH4"))
NS <- do.call(rbind, lapply(names(TABS), function(tn) { t <- TABS[[tn]]; seeds <- t$s0 + 0:199
  fs <- file.path(root, t$null, paste0(tn, "_seed", seeds, "_ancombc2.csv")); fs <- fs[file.exists(fs)][seq_len(min(NSEED, sum(file.exists(fs))))]
  if (NSEED == 200) stopifnot("expected 200 saved null fits" = length(fs) == 200)
  do.call(rbind, lapply(fs, function(f) { s <- fitstats(f, t$trait); if (is.null(s)) stop("passed_ss or diff missing in ", f); cbind(table=tn, seed=as.integer(sub(".*_seed([0-9]+)_.*", "\\1", f)), s) })) }))
cp <- function(x) { n <- length(x); s <- sum(x); ci <- binom.test(s, n)$conf.int; c(s, n, s / n, ci) }
NSUM <- do.call(rbind, lapply(split(NS, NS$table), function(x) { a <- cp(x$any_bh); b <- cp(x$any_robust)
  data.frame(table=x$table[1], shuffles=nrow(x), n_genera=median(x$n_genera),
    median_share_p05=median(x$share_p05), min_share_p05=min(x$share_p05), max_share_p05=max(x$share_p05), mean_share_p05=mean(x$share_p05),
    median_share_p05_passed=median(x$share_p05_passed_all), min_share_p05_passed=min(x$share_p05_passed_all), max_share_p05_passed=max(x$share_p05_passed_all), mean_share_p05_passed=mean(x$share_p05_passed_all),
    median_share_p05_among_passed=median(x$share_p05_among_passed, na.rm=TRUE), median_share_passed=median(x$share_passed),
    median_n_bh=median(x$n_bh), max_n_bh=max(x$n_bh), median_n_robust=median(x$n_robust), max_n_robust=max(x$n_robust),
    P_any_bh=a[3], P_any_bh_lo=a[4], P_any_bh_hi=a[5], P_any_robust=b[3], P_any_robust_lo=b[4], P_any_robust_hi=b[5], diff_vs_q_mismatches=sum(x$diff_vs_q_mismatch)) }))
NSUM <- NSUM[match(names(TABS), NSUM$table), ]
RF <- do.call(rbind, lapply(names(TABS), function(tn) { f <- file.path(root, TABS[[tn]]$real); s <- fitstats(f, TABS[[tn]]$trait)
  if (is.null(s)) data.frame(table=tn, file=TABS[[tn]]$real, note="passed_ss or diff not kept") else cbind(table=tn, file=TABS[[tn]]$real, s, note="") }))
# positive controls (Table 1): the saved fits under results/positive_controls/{profiles_gtdb_full,16s}/; the control terms of positive_controls_summary.tsv and 16s/*_summary.tsv
pcs <- list.files(file.path(root, "results", "positive_controls"), "_ancombc2\\.csv$", recursive=TRUE, full.names=TRUE)
pcs <- pcs[grepl("/(profiles_gtdb_full|16s)/", pcs)]   # the Table 1 fits (profiles_gtdb is the superseded profile)
PC <- do.call(rbind, lapply(pcs, function(f) { x <- rc(f); terms <- unique(sub("^q_", "", grep("^q_", names(x), value=TRUE))); terms <- terms[terms != "(Intercept)"]
  do.call(rbind, lapply(terms, function(t) { has <- all(paste0(c("diff_", "passed_ss_"), t) %in% names(x))
    q <- x[[paste0("q_", t)]]; if (!has) return(data.frame(file=sub(paste0(root, "/"), "", f), term_column=t, n_taxa=nrow(x), n_q05=sum(q < 0.05, na.rm=TRUE), n_robust=NA, note="passed_ss not kept"))
    ss <- as.logical(x[[paste0("passed_ss_", t)]]); data.frame(file=sub(paste0(root, "/"), "", f), term_column=t, n_taxa=nrow(x), n_q05=sum(q < 0.05, na.rm=TRUE), n_robust=sum(q < 0.05 & ss, na.rm=TRUE), note="") })) }))
if (is.null(PC)) PC <- data.frame(file=NA, term_column=NA, n_taxa=NA, n_q05=NA, n_robust=NA, note="no saved positive-control fit found under results/positive_controls/")

# ---------------------------------------------------------------- (c) Charolais genus: reproduce the primary call, then group = RFI_class, struc_zero = TRUE
SZ <- data.frame(step=NA, note="not run"); SZT <- data.frame(taxon=NA, note="not run")
if (has_pkg && Sys.getenv("A48_SKIP_SZ") != "1") { suppressMessages({ library(ANCOMBC); library(TreeSummarizedExperiment) })
  d <- a46_load(root, "charolais_mrna_genus"); A <- as.matrix(d$A); orig <- rownames(A); rownames(A) <- sprintf("F%05d", seq_along(orig)); md <- d$md
  tse <- TreeSummarizedExperiment(assays=list(counts=A), colData=md[, d$terms, drop=FALSE])
  set.seed(20261102); f0 <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(d$terms, collapse=" + "), p_adj_method="BH", prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
  r0 <- f0$res; r0$taxon_name <- orig[match(r0$taxon, rownames(A))]; prim <- rc(file.path(root, "results", "replication", "charolais_mrna_genus_ancombc2.csv"))
  wt <- grep("^W_RFI", names(prim), value=TRUE); stopifnot(length(wt) == 1); j <- match(prim$taxon_name, r0$taxon_name)
  dW <- max(abs(prim[[wt]] - r0[[wt]][j]), na.rm=TRUE); cat("== (c) primary Charolais genus fit reproduced: max |dW| =", signif(dW, 3), "over", sum(!is.na(j)), "taxa\n")
  if (anyNA(j) || nrow(prim) != nrow(r0) || dW > 1e-6) stop("48.2(c): the primary Charolais genus fit is not reproduced (max |dW| ", signif(dW, 3), "); stopping as the amendment requires")
  set.seed(20261102); f1 <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(d$terms, collapse=" + "), p_adj_method="BH", prv_cut=0.10, group="RFI_class", struc_zero=TRUE, verbose=FALSE)
  r1 <- f1$res; r1$taxon_name <- orig[match(r1$taxon, rownames(A))]; z <- f1$zero_ind; z$taxon_name <- orig[match(z$taxon, rownames(A))]
  zc <- grep("^structural_zero", names(z), value=TRUE); zf <- z[rowSums(as.matrix(z[, zc, drop=FALSE])) > 0, c("taxon_name", zc), drop=FALSE]
  tg <- unique(d$feat[!is.na(d$feat)]); k0 <- match(tg, r0$taxon_name); k1 <- match(tg, r1$taxon_name); pq <- function(r, k, p) if (p %in% names(r)) r[[p]][k] else NA
  SZ <- data.frame(target_feature=tg, W_primary=r0[[wt]][k0], W_struc_zero=pq(r1, k1, wt), p_primary=pq(r0, k0, sub("^W_", "p_", wt)), p_struc_zero=pq(r1, k1, sub("^W_", "p_", wt)),
    q_primary=pq(r0, k0, sub("^W_", "q_", wt)), q_struc_zero=pq(r1, k1, sub("^W_", "q_", wt)), flagged_structural_zero=tg %in% zf$taxon_name, primary_reproduced_max_dW=dW,
    taxa_flagged=nrow(zf), taxa_primary=nrow(r0), taxa_struc_zero_fit=nrow(r1), note="sensitivity only; sets nothing")
  SZT <- if (nrow(zf)) cbind(zf, note="flagged as structural zero in a level of RFI_class") else data.frame(taxon_name=NA, note="no taxon flagged") }

# ---------------------------------------------------------------- (d) reading, fixed in 48.2
ok <- NSUM$median_share_p05_passed <= 0.07
RD <- data.frame(rule="median share (p < 0.05 and passed_ss) <= 0.07 in all five tables", tables_meeting=sum(ok), tables=nrow(NSUM),
  tables_failing=paste(NSUM$table[!ok], collapse=";"),
  reading=if (all(ok)) "the inflation concerns unfiltered nominal p values: report filtered rates beside Table 2 and replace 'poorly calibrated'/'miscalibrated' (Background, Results heading, Fig. 3 title, Conclusions) by what was inflated under which settings" else
    "the inflation statement stands, with the configuration stated and the filtered rates reported beside it",
  note="the permutation calibration of W used for effect sizes and labels is unchanged; the filter is a decision rule for calling taxa, not a correction of W")
a48_write(calls, outs[1]); a48_write(LP, outs[2]); a48_write(NS, outs[3]); a48_write(NSUM, outs[4]); a48_write(RF, outs[5]); a48_write(PC, outs[6]); a48_write(SZ, outs[7]); a48_write(SZT, outs[8]); a48_write(RD, outs[9])
cat("== (b) null fits, per table\n"); print(NSUM[, c("table", "shuffles", "median_share_p05", "median_share_p05_passed", "max_share_p05_passed", "median_n_bh", "max_n_bh", "median_n_robust", "max_n_robust", "P_any_bh", "P_any_robust", "P_any_robust_hi", "diff_vs_q_mismatches")], digits=3, row.names=FALSE)
cat("\n== real fits\n"); print(RF, digits=3, row.names=FALSE); cat("\n== positive controls\n"); print(PC, row.names=FALSE)
cat("\n== (c) Charolais structural zeros\n"); print(SZ, digits=3, row.names=FALSE); cat("\n== (d) reading\n"); print(RD[, c("tables_meeting", "tables", "tables_failing", "reading")], row.names=FALSE)
