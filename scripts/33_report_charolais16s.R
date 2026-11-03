#!/usr/bin/env Rscript
# Amendment 17 report: Charolais 16S against Charolais mRNA (same animals). Descriptive; no label is assigned.
setwd("/data/rumen_multicohort"); rd <- function(f) read.delim(f, stringsAsFactors=FALSE, quote="")
g <- function(f) { x <- rd(f); x$r <- x$W / sqrt(x$W^2 + x$df); x$p <- 2 * pt(-abs(x$W), x$df); x }
out <- list()
for (rk in c("genus", "family")) {
  a <- g(sprintf("results/charolais16s/charolais_16s_%s_targets.tsv", rk)); b <- g(sprintf("results/replication/charolais_mrna_%s_targets.tsv", rk))
  m <- merge(a[, c("id", "set", "rank", "feature", "status", "r", "p")], b[, c("id", "rank", "status", "r", "p")], by=c("id", "rank"), suffixes=c("_16s", "_mrna"))
  m$both <- m$status_16s == "ok" & m$status_mrna == "ok"; out[[rk]] <- m; b2 <- m[m$both, ]
  cat(rk, ": targets", nrow(m), "| matched to a 16S feature", sum(!is.na(m$feature)), "| estimable 16S", sum(m$status_16s == "ok"), "| estimable mRNA", sum(m$status_mrna == "ok"), "| both", nrow(b2), "\n")
  if (nrow(b2) > 2) cat("  same sign", sum(sign(b2$r_16s) == sign(b2$r_mrna)), "of", nrow(b2), "| Spearman r", round(cor(b2$r_16s, b2$r_mrna, method="spearman"), 2), "| 16S p<0.05:", sum(b2$p_16s < 0.05), "\n")
}
write.table(do.call(rbind, out), "results/charolais16s/comparison_with_mrna.tsv", sep="\t", quote=FALSE, row.names=FALSE)
