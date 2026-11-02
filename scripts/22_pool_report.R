#!/usr/bin/env Rscript
# Pool, classify, report replication targets (Amendments 14, 14b, 14c, 14d). Usage: 22_pool_report.R
suppressMessages(library(metafor))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "replication")
sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
stopifnot("expected 8 target files (4 ranks x 2 cohorts); a fit is still running" = length(list.files(od, "_targets[.]tsv$")) == 8)
est <- do.call(rbind, lapply(list.files(od, "_targets\\.tsv$", full.names=TRUE), rd)); est <- est[est$status == "ok", ]
est$r <- partial_r(est$W, est$df); z <- atanh(est$r); s <- sqrt(1/(est$n - 3 - Q))
est$lo <- tanh(z - 1.96*s); est$hi <- tanh(z + 1.96*s); est$p_t <- 2*pt(-abs(est$W), est$df)
est$p_flag <- abs(log10(pmax(est$p_t, 1e-300)) - log10(pmax(est$p_ancom, 1e-300))) > 1 & pmax(est$p_t, est$p_ancom) > 1e-6
pdir <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available",
  ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart",
  ifelse(tk$site == "feces", "cross-site", "replication")))
cv <- do.call(rbind, lapply(list.files(file.path(root, "logs", "replication"), "_dry\\.log$", full.names=TRUE), function(f) {
  l <- readLines(f); m <- regmatches(l, regexec("^ *(\\S+) .* (Li|Canonical|Nelore) +([0-9]+) +([0-9]+) ", l)); m <- m[lengths(m) == 5]
  if (length(m)) data.frame(id=sapply(m, `[`, 2), frac=as.integer(sapply(m, `[`, 5)) / as.integer(sapply(m, `[`, 4))) }))
cov <- tapply(cv$frac, cv$id, min)
all <- list(); lines <- list()
for (st in unique(tk$set)) {
  d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pdir(d$direction), testable=d$stratum %in% c("replication", "cross-site"))
  e <- est[est$id %in% d$id, ]; res <- classify_set(tg, e[, c("id", "cohort", "W", "df", "n")])
  res$q_testable <- NA; res$q_testable[tg$testable] <- p.adjust(res$p[tg$testable], "BH")
  res <- cbind(set=st, target=d$target, rank=d$target_rank, site=d$site, stratum=d$stratum, published_dir=tg$dir, res[, -1, drop=FALSE])
  res <- cbind(id=d$id, res); res$cross_site <- d$site == "feces"
  res$coverage <- round(as.numeric(cov[d$id]), 2); res$partial_coverage <- res$coverage < 1
  for (cc in c("nelore", "charolais_mrna")) { x <- e[e$cohort == cc, ]; j <- match(d$id, x$id)
    res[[paste0("r_", cc)]] <- x$r[j]; res[[paste0("lo_", cc)]] <- x$lo[j]; res[[paste0("hi_", cc)]] <- x$hi[j] }
  jc <- match(d$id, e$id[e$cohort == "charolais_mrna"]); xc <- e[e$cohort == "charolais_mrna", ]
  res$g_charolais <- xc$W[jc] * sqrt(1/14 + 1/16) * (1 - 3/(4*xc$df[jc] - 1))
  res$p_ancom_disagree <- d$id %in% e$id[e$p_flag]
  ref <- ifelse(is.na(res$published_dir), NA, 0.3 * res$published_dir)
  res$large_effect_not_excluded <- startsWith(res$label, "not replicated") & res$ci_lo <= 0 & res$ci_hi >= 0 &
    ifelse(is.na(ref), res$ci_hi >= 0.3 | res$ci_lo <= -0.3, ref >= res$ci_lo & ref <= res$ci_hi)
  res$reason <- ifelse(res$label == "untestable", ifelse(tg$testable, "below detection, filtered or no estimate", d$stratum), "")
  all[[st]] <- res
}
all <- do.call(rbind, all); write.table(all, file.path(od, "classification_all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(table(all$set, all$label)); print(table(all$stratum, all$label))
print(all[all$label != "untestable", c("id", "target", "stratum", "k", "r", "q_bh", "q_testable", "label")], digits=3, row.names=FALSE)
pdf(file.path(od, "forest_replication.pdf"), width=8, height=10)
for (st in unique(all$set)) for (sg in c("replication", "cross-site")) {
  d <- all[all$set == st & all$stratum == sg & all$k > 0, ]; if (!nrow(d)) next
  L <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) { x <- d[i, ]; e <- est[est$id == x$id, ]
    rbind(data.frame(lab=paste(x$id, x$target, if (x$k > 1) "pooled" else e$cohort[1]), r=x$r, lo=x$ci_lo, hi=x$ci_hi, pooled=TRUE, lb=x$label),
          data.frame(lab=paste0("   ", e$cohort), r=e$r, lo=e$lo, hi=e$hi, pooled=FALSE, lb=""))[if (x$k > 1) 1:(1 + nrow(e)) else 1, ] }))
  y <- rev(seq_len(nrow(L))); par(mar=c(4, 17, 3, 2), cex.axis=0.6)
  plot(NA, xlim=c(-1, 1), ylim=c(0.5, nrow(L) + 0.5), yaxt="n", xlab="partial r (positive = more abundant with higher RFI)", ylab="", main=paste(st, "set,", sg, "stratum"))
  abline(v=0, lty=2); segments(L$lo, y, L$hi, y, col=ifelse(L$pooled, "black", "grey50"))
  points(L$r, y, pch=ifelse(L$pooled, 18, 16), cex=ifelse(L$pooled, 1.4, 0.8), col=ifelse(L$pooled, "black", "grey40"))
  axis(2, at=y, labels=L$lab, las=1, tick=FALSE); text(1, y[L$pooled], L$lb[L$pooled], pos=2, cex=0.45, col="firebrick")
}
invisible(dev.off()); cat("wrote", file.path(od, "classification_all.tsv"), "and forest_replication.pdf\n")
