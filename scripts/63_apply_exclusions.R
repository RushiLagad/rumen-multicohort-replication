#!/usr/bin/env Rscript
# Amendment 22: apply the 16S sample exclusion rules. Usage: 63_apply_exclusions.R <BioProject> <cohort> [wrong-marker run IDs, comma separated]. Read data only.
a <- commandArgs(TRUE); P <- a[1]; coh <- a[2]; wm <- if (length(a) > 2) strsplit(a[3], ",")[[1]] else character(0)
setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")); R <- file.path("results/16s", P); d <- file.path("profiles_16s", coh)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
trim <- read.delim(file.path(R, "trim_summary.tsv"), header=FALSE, col.names=c("run", "pairs_in", "pairs_kept", "frac_kept"), stringsAsFactors=FALSE)
trk <- rd(file.path(R, "read_tracking.tsv")); gf <- file.path(d, paste0(coh, "_genus_counts.tsv")); uf <- sub("\\.tsv$", "_unfiltered.tsv", gf)
if (!file.exists(uf)) invisible(file.copy(gf, uf)); G <- rd(uf)
stopifnot(all(wm %in% trim$run), all(trim$frac_kept[match(wm, trim$run)] < 0.10))
ex <- rbind(data.frame(run=wm, rule="1 wrong marker", value=round(trim$frac_kept[match(wm, trim$run)], 4)),
            data.frame(run=trk$run[trk$nonchim < 1000], rule="2 under 1000 non-chimeric reads", value=trk$nonchim[trk$nonchim < 1000]))
ex <- ex[!duplicated(ex$run), ]; lost <- setdiff(trim$run, trk$run)
if (length(lost)) ex <- rbind(ex, data.frame(run=lost, rule="no reads after filtering", value=0))
write.table(ex, file.path(R, "excluded_runs.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
keep <- setdiff(names(G)[-1], ex$run); G <- G[, c(names(G)[1], keep)]; G <- G[rowSums(G[, -1, drop=FALSE]) > 0, ]
write.table(G, gf, sep="\t", quote=FALSE, row.names=FALSE)
cat(P, ": excluded", nrow(ex), "of", nrow(trim), "runs |", paste(names(table(ex$rule)), table(ex$rule), sep=": ", collapse="; "), "| kept", length(keep), "samples,", nrow(G), "genera\n")
