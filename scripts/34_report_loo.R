#!/usr/bin/env Rscript
# Amendment 18 (exploratory): influence of single animals on the Charolais mRNA target estimates.
setwd("/data/rumen_multicohort"); rd <- function(f) read.delim(f, stringsAsFactors=FALSE, quote="")
rp <- function(x) { x$r <- x$W / sqrt(x$W^2 + x$df); x$p <- 2 * pt(-abs(x$W), x$df); x }
res <- list()
for (rk in c("genus", "family", "phylum")) {
  fs <- Sys.glob(sprintf("results/loo/drop_*/charolais_mrna_%s_targets.tsv", rk)); cat(rk, ": fits found", length(fs), "of 30\n"); if (!length(fs)) next
  b <- rp(rd(sprintf("results/replication/charolais_mrna_%s_targets.tsv", rk))); b <- b[b$status == "ok", ]
  d <- do.call(rbind, lapply(fs, function(f) { x <- rp(rd(f)); x$drop <- sub("^drop_", "", basename(dirname(f))); x[x$status == "ok", c("id", "drop", "r", "p")] }))
  for (i in seq_len(nrow(b))) {
    z <- d[d$id == b$id[i], ]; if (!nrow(z)) next; j <- which.max(abs(z$r - b$r[i]))
    res[[length(res) + 1]] <- data.frame(rank=rk, id=b$id[i], set=b$set[i], n_fits=nrow(z), r_full=b$r[i], r_min=min(z$r), r_max=max(z$r), max_abs_change=abs(z$r[j] - b$r[i]),
      most_influential=z$drop[j], sign_flips=sum(sign(z$r) != sign(b$r[i])), p_full=b$p[i], share_p05=if (b$p[i] < 0.05) mean(z$p < 0.05) else NA)
  }
}
stopifnot(length(res) > 0); o <- do.call(rbind, res); write.table(o, "results/loo/influence_charolais_mrna.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("targets", nrow(o), "| any sign flip", sum(o$sign_flips > 0), "| median max change in r", round(median(o$max_abs_change), 3), "| max", round(max(o$max_abs_change), 3), "\n")
print(o[o$sign_flips > 0 | (!is.na(o$share_p05) & o$share_p05 < 1), c("rank", "id", "r_full", "sign_flips", "max_abs_change", "most_influential", "share_p05")], row.names=FALSE)
