#!/usr/bin/env Rscript
# Amendment 37 item 4: RFI arm pooled as in 66_pool_rfi3.R, but each cohort W is first calibrated by its 200-permutation null
# (W_cal = W / max(1, SD_null), as in 83_pool_methane.R). Writes results/replication_rfi3_calibrated/; 66's outputs are kept as the uncalibrated version.
suppressMessages(library(metafor))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "replication_rfi3_calibrated"); dir.create(od, showWarnings=FALSE)
sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("commit Amendment 37 first" = any(grepl("^\\*\\*Amendment 37 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
nulldir <- file.path(root, "results", "rfi_null"); stopifnot("RFI nulls not finished" = file.exists(file.path(root, "logs", "rfi_null", "DONE")))
tk <- rd(file.path(root, "plan", "target_keys.tsv")); tk <- tk[!tk$set %in% c("Methane", "MethaneSecondary"), ]  # RFI sets only
base <- list.files(file.path(root, "results", "replication"), "_targets\\.tsv$", full.names=TRUE); stopifnot(length(base) == 8)
pdir <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available",
  ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
# ---- null SD of W per table and target (seed 1001-1200 permutations)
nullsd <- function(tab) { fs <- list.files(nulldir, paste0("^", tab, "_seed1[0-9]{3}_targets\\.tsv$"), full.names=TRUE)
  if (!length(fs)) return(NULL); x <- do.call(rbind, lapply(fs, function(f) { d <- rd(f); d <- d[d$status == "ok", c("id", "W")]; d }))
  s <- aggregate(W ~ id, x, function(w) c(sd=sd(w), k=length(w))); data.frame(id=s$id, sd_null=s$W[, "sd"], k_null=s$W[, "k"], W_null=I(split(x$W, x$id)[s$id])) }
calibrate <- function(f) { tab <- sub("_targets\\.tsv$", "", basename(f)); e <- rd(f); e <- e[e$status == "ok", ]; if (!nrow(e)) return(e)
  ns <- nullsd(tab); stopifnot("no permutation null for table" = !is.null(ns)); j <- match(e$id, ns$id)
  stopifnot("target missing from null" = !anyNA(j)); stopifnot("fewer than 190 null fits" = all(ns$k_null[j] >= 190))
  e$table <- tab; e$sd_null <- ns$sd_null[j]; e$k_null <- ns$k_null[j]; e$W_raw <- e$W
  e$p_perm <- mapply(function(w, i) mean(abs(ns$W_null[[i]]) >= abs(w)), e$W, j)
  e$W <- e$W_raw / pmax(1, e$sd_null); e }
run <- function(tag) {
  af <- file.path(root, "results", "alberta", paste0("alberta_", c("genus", "family"), tag, "_targets.tsv")); stopifnot(all(file.exists(af)))
  est <- do.call(rbind, lapply(c(base, af), calibrate)); stopifnot("a target appears twice for one cohort" = !anyDuplicated(est[, c("id", "cohort")]))
  est$r <- partial_r(est$W, est$df); z <- atanh(est$r); s <- sqrt(1/(est$df - 1)); est$lo <- tanh(z - 1.96*s); est$hi <- tanh(z + 1.96*s)
  write.table(est, file.path(od, paste0("estimates_calibrated", tag, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  coh <- unique(est$cohort); all <- list()
  for (st in unique(tk$set)) {
    d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pdir(d$direction), testable=d$stratum %in% c("replication", "cross-site"))
    e <- est[est$id %in% d$id, ]; res <- classify_set(tg, e[, c("id", "cohort", "W", "df", "n")])
    res <- cbind(set=st, target=d$target, rank=d$target_rank, stratum=d$stratum, published_dir=tg$dir, res)
    for (cc in coh) { x <- e[e$cohort == cc, ]; j <- match(d$id, x$id)
      res[[paste0("r_", cc)]] <- x$r[j]; res[[paste0("lo_", cc)]] <- x$lo[j]; res[[paste0("hi_", cc)]] <- x$hi[j]
      res[[paste0("sd_null_", cc)]] <- x$sd_null[j]; res[[paste0("p_perm_", cc)]] <- x$p_perm[j] }
    ref <- ifelse(is.na(res$published_dir), NA, 0.3 * res$published_dir)
    res$large_effect_not_excluded <- startsWith(res$label, "not replicated") & res$ci_lo <= 0 & res$ci_hi >= 0 &
      ifelse(is.na(ref), res$ci_hi >= 0.3 | res$ci_lo <= -0.3, ref >= res$ci_lo & ref <= res$ci_hi)
    all[[st]] <- res }
  all <- do.call(rbind, all); write.table(all, file.path(od, paste0("classification", if (nzchar(tag)) tag else "_primary", ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  list(cls=all, est=est) }
P <- run(""); S <- run("_all")
cat("== null SD of W by table (median over targets; >1 means ANCOM-BC2 was anti-conservative there)\n")
print(aggregate(sd_null ~ table, P$est, function(x) round(c(median=median(x), max=max(x), n=length(x)), 2)), row.names=FALSE)
cat("== calibration factors applied (sd_null > 1):\n"); x <- P$est[P$est$sd_null > 1, c("id", "table", "sd_null", "W_raw", "W", "p_perm")]; print(x[order(-x$sd_null), ], digits=2, row.names=FALSE)
cat("== PRIMARY, calibrated\n"); print(table(P$cls$set, P$cls$label))
x <- P$cls[!startsWith(P$cls$label, "untestable"), ]
print(x[order(x$set, x$q_bh), c("id", "target", "set", "k", "r", "ci_lo", "ci_hi", "q_bh", "p_perm_nelore", "p_perm_charolais_mrna", "p_perm_alberta", "label")], digits=2, row.names=FALSE)
old <- rd(file.path(root, "results", "replication_rfi3", "classification_primary.tsv")); j <- match(P$cls$id, old$id)
old <- rd(file.path(root, "results", "replication_rfi3", "classification_primary.tsv")); j <- match(P$cls$id, old$id); ok <- !is.na(j)
chg <- ok & old$label[j] != P$cls$label
cat("== labels that change, uncalibrated (66) -> calibrated (89), primary:", sum(chg, na.rm=TRUE), "of", sum(ok), "\n")
if (any(chg, na.rm=TRUE)) print(data.frame(id=P$cls$id, target=P$cls$target, uncalibrated=old$label[j], calibrated=P$cls$label)[which(chg), ], row.names=FALSE)
cat("== pooled r and q, uncalibrated -> calibrated, testable targets:\n")
dr <- data.frame(id=P$cls$id, target=P$cls$target, k=P$cls$k, r_old=old$r[j], r_new=P$cls$r, lo_new=P$cls$ci_lo, hi_new=P$cls$ci_hi, q_old=old$q_bh[j], q_new=P$cls$q_bh)[ok & is.finite(P$cls$r), ]
print(dr[order(dr$q_new), ], digits=2, row.names=FALSE)
cat("== primary vs sensitivity (all 567) label differences:", sum(P$cls$label != S$cls$label), "\n")
cat("== not-replicated targets whose CI still includes the published r = 0.3 (calibrated primary):", sum(P$cls$large_effect_not_excluded, na.rm=TRUE), "of", sum(startsWith(P$cls$label, "not replicated")), "\n")
cat("   same count, uncalibrated (66):", sum(old$large_effect_not_excluded, na.rm=TRUE), "of", sum(startsWith(old$label, "not replicated")), "\n")
