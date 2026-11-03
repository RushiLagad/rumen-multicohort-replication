#!/usr/bin/env Rscript
# Amendment 17 step 3: DADA2 on forward reads -> genus counts. No phenotype data are loaded here.
suppressMessages({library(dada2); library(ShortRead)})
setwd("/data/rumen_multicohort"); set.seed(20261102)
O <- "work/charolais16s"; R <- "results/charolais16s"; nt <- 20
runs <- readLines(file.path(O, "runs.txt")); fq <- file.path(O, "trimmed", paste0(runs, "_1.fastq.gz")); stopifnot(length(runs) == 30, all(file.exists(fq)))
W <- 260; pool <- do.call(rbind, lapply(fq, function(f) { s <- FastqSampler(f, 20000); m <- as(quality(yield(s)), "matrix"); close(s)
  m <- m[, seq_len(min(ncol(m), W)), drop=FALSE]; if (ncol(m) < W) m <- cbind(m, matrix(NA_real_, nrow(m), W - ncol(m))); m }))
med <- apply(pool, 2, median, na.rm=TRUE); write.table(data.frame(pos=1:W, medianQ=med), file.path(R, "median_quality_R1.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
tl <- max(which(med[1:220] >= 30)); cat("truncLen", tl, "\n"); stopifnot(is.finite(tl), tl >= 150)
ft <- file.path(O, "filtered", paste0(runs, "_1.fastq.gz")); dir.create(dirname(ft[1]), showWarnings=FALSE)
out <- filterAndTrim(fq, ft, truncLen=tl, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=nt)
stopifnot(all(out[, "reads.out"] > 0), all(file.exists(ft)))
err <- learnErrors(ft, multithread=nt, randomize=TRUE); dd <- dada(ft, err=err, multithread=nt, verbose=0); names(dd) <- runs
st <- makeSequenceTable(dd); st <- removeBimeraDenovo(st, method="consensus", multithread=nt); rownames(st) <- runs
tx <- assignTaxonomy(colnames(st), file.path(O, "silva138_train.fa"), minBoot=50, multithread=nt, taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
g <- unname(tx[, "Genus"]); g <- gsub("[][]", "", g); g <- sub(" group$", "", g); g <- sub(" [0-9]+$", "", g)
ok <- !is.na(g) & nzchar(g) & !grepl("_ge$|^uncultured|^unidentified", g, ignore.case=TRUE)
G <- rowsum(t(st)[ok, , drop=FALSE], g[ok]); stopifnot(nrow(G) > 50, ncol(G) == 30)
d <- "profiles_16s/charolais_16s"; dir.create(d, recursive=TRUE, showWarnings=FALSE)
write.table(data.frame(genus=rownames(G), G, check.names=FALSE), file.path(d, "charolais_16s_genus_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
trk <- data.frame(run=runs, input=out[, 1], filtered=out[, 2], nonchim=rowSums(st), classified_genus=colSums(G), prop_classified=round(colSums(G) / rowSums(st), 3))
write.table(trk, file.path(R, "read_tracking.tsv"), sep="\t", quote=FALSE, row.names=FALSE); saveRDS(list(st=st, tx=tx, truncLen=tl), file.path(R, "dada2.rds"))
cat("genera", nrow(G), " ASVs", ncol(st), " median classified share", median(trk$prop_classified), "\n")
