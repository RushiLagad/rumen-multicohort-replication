#!/usr/bin/env Rscript
# Amendments 19, 24, 25: RFI arm pooled over Nelore, Charolais mRNA and Alberta; primary and sensitivity (alberta_all). Usage: 66_pool_rfi3.R
suppressMessages(library(metafor))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "replication_rfi3"); dir.create(od, showWarnings=FALSE)
sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("commit Amendment 25 first" = any(grepl("^\\*\\*Amendment 25 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
base <- list.files(file.path(root, "results", "replication"), "_targets\\.tsv$", full.names=TRUE); stopifnot(length(base) == 8)
pdir <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available",
  ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
run <- function(tag) {
  af <- file.path(root, "results", "alberta", paste0("alberta_", c("genus", "family"), tag, "_targets.tsv")); stopifnot(all(file.exists(af)))
  est <- do.call(rbind, lapply(c(base, af), rd)); est <- est[est$status == "ok", ]; stopifnot("a target appears twice for one cohort" = !anyDuplicated(est[, c("id", "cohort")]))
  est$r <- partial_r(est$W, est$df); z <- atanh(est$r); s <- sqrt(1/(est$df - 1)); est$lo <- tanh(z - 1.96*s); est$hi <- tanh(z + 1.96*s)
  coh <- unique(est$cohort); all <- list()
  for (st in unique(tk$set)) {
    d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pdir(d$direction), testable=d$stratum %in% c("replication", "cross-site"))
    e <- est[est$id %in% d$id, ]; res <- classify_set(tg, e[, c("id", "cohort", "W", "df", "n")])
    res <- cbind(set=st, target=d$target, rank=d$target_rank, stratum=d$stratum, published_dir=tg$dir, res)
    for (cc in coh) { x <- e[e$cohort == cc, ]; j <- match(d$id, x$id); res[[paste0("r_", cc)]] <- x$r[j]; res[[paste0("lo_", cc)]] <- x$lo[j]; res[[paste0("hi_", cc)]] <- x$hi[j] }
    ref <- ifelse(is.na(res$published_dir), NA, 0.3 * res$published_dir)
    res$large_effect_not_excluded <- startsWith(res$label, "not replicated") & res$ci_lo <= 0 & res$ci_hi >= 0 &
      ifelse(is.na(ref), res$ci_hi >= 0.3 | res$ci_lo <= -0.3, ref >= res$ci_lo & ref <= res$ci_hi)
    all[[st]] <- res }
  all <- do.call(rbind, all); write.table(all, file.path(od, paste0("classification", if (nzchar(tag)) tag else "_primary", ".tsv")), sep="\t", quote=FALSE, row.names=FALSE); all }
P <- run(""); S <- run("_all")
cat("== PRIMARY (Alberta without the 35 PRJNA448333 ID numbers)\n"); print(table(P$set, P$label))
kt <- table(P$k[P$k > 0]); cat("== testable targets by number of cohorts (primary):", paste0("k=", names(kt), ": ", kt, collapse="; "), "\n")
x <- P[!startsWith(P$label, "untestable"), ]
print(x[order(x$set, x$q_bh), c("id", "target", "set", "k", "r", "ci_lo", "ci_hi", "q_bh", "loo_all_sig", "label")], digits=2, row.names=FALSE)
cat("== targets whose label differs between primary and sensitivity (all 567):", sum(P$label != S$label), "\n")
if (any(P$label != S$label)) print(data.frame(id=P$id, primary=P$label, sensitivity=S$label)[P$label != S$label, ], row.names=FALSE)
cat("== not-replicated targets whose CI still includes the published r = 0.3 (primary):", sum(P$large_effect_not_excluded, na.rm=TRUE), "of", sum(startsWith(P$label, "not replicated")), "\n")
