#!/usr/bin/env Rscript
# Amendments 17, 19, 21: DADA2 for a paired 16S cohort -> genus counts. Usage: 61_dada2_16s_pe.R <BioProject> <cohort> <expected insert bp>. No phenotype data loaded.
suppressMessages({library(dada2); library(ShortRead)})
a <- commandArgs(TRUE); P <- a[1]; coh <- a[2]; INSERT <- as.numeric(a[3]); setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")); set.seed(20261103); nt <- 20
O <- file.path("work/16s", P); R <- file.path("results/16s", P); dir.create(R, recursive=TRUE, showWarnings=FALSE)
runs <- readLines(file.path(O, "runs.txt")); f1 <- file.path(O, "trimmed", paste0(runs, "_1.fastq.gz")); f2 <- sub("_1\\.fastq\\.gz$", "_2.fastq.gz", f1)
ok <- file.exists(f1) & file.exists(f2) & file.size(f1) > 0; runs <- runs[ok]; f1 <- f1[ok]; f2 <- f2[ok]; cat("runs with trimmed pairs:", length(runs), "\n")
medq <- function(fs, W = 320, k = 5000) apply(do.call(rbind, lapply(fs, function(f) { s <- FastqSampler(f, k); q <- as(quality(yield(s)), "matrix"); close(s)
  q <- q[, seq_len(min(ncol(q), W)), drop=FALSE]; if (ncol(q) < W) q <- cbind(q, matrix(NA_real_, nrow(q), W - ncol(q))); q })), 2, median, na.rm=TRUE)
last <- function(q, cap, thr) { i <- which(q[seq_len(cap)] >= thr); if (length(i)) max(i) else NA }
q1 <- medq(f1); q2 <- medq(f2); t1 <- last(q1, 280, 30); t2 <- last(q2, 260, 25)
while (!is.na(t1) && !is.na(t2) && t1 + t2 - INSERT < 40 && t2 < 260) t2 <- t2 + 1
write.table(data.frame(pos=seq_along(q1), R1=q1, R2=q2), file.path(R, "median_quality.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("truncLen R1", t1, "R2", t2, "expected overlap", t1 + t2 - INSERT, "\n"); stopifnot(!is.na(t1), !is.na(t2), t1 >= 150, t2 >= 100)
dir.create(file.path(O, "filtered"), showWarnings=FALSE); g1 <- file.path(O, "filtered", basename(f1)); g2 <- file.path(O, "filtered", basename(f2))
out <- filterAndTrim(f1, g1, f2, g2, truncLen=c(t1, t2), maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=nt)
keep <- out[, 2] > 0 & file.exists(g1); runs <- runs[keep]; f1 <- f1[keep]; g1 <- g1[keep]; g2 <- g2[keep]; out <- out[keep, , drop=FALSE]
e1 <- learnErrors(g1, multithread=nt, randomize=TRUE); e2 <- learnErrors(g2, multithread=nt, randomize=TRUE)
d1 <- dada(g1, err=e1, multithread=nt, verbose=0); d2 <- dada(g2, err=e2, multithread=nt, verbose=0)
mg <- mergePairs(d1, g1, d2, g2, minOverlap=20, maxMismatch=0); names(mg) <- runs
mr <- vapply(mg, function(x) sum(x$abundance), 0) / out[, 2]; cat("median merge rate", round(median(mr), 3), "\n"); mode <- "merged"
if (median(mr) < 0.5) {
  mode <- "forward"; tf <- last(q1, 220, 30); cat("FALLBACK to forward reads (Amendment 21 item 4); truncLen", tf, "\n"); stopifnot(!is.na(tf), tf >= 150)
  dir.create(file.path(O, "filtered_fwd"), showWarnings=FALSE); h1 <- file.path(O, "filtered_fwd", basename(f1))
  out <- filterAndTrim(f1, h1, truncLen=tf, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=nt)
  k2 <- out[, 2] > 0 & file.exists(h1); runs <- runs[k2]; h1 <- h1[k2]; out <- out[k2, , drop=FALSE]
  ef <- learnErrors(h1, multithread=nt, randomize=TRUE); mg <- dada(h1, err=ef, multithread=nt, verbose=0); names(mg) <- runs }
st <- makeSequenceTable(mg); st <- removeBimeraDenovo(st, method="consensus", multithread=nt); rownames(st) <- runs
tx <- assignTaxonomy(colnames(st), "work/charolais16s/silva138_train.fa", minBoot=50, multithread=nt, taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
g <- unname(tx[, "Genus"]); g <- gsub("[][]", "", g); g <- sub(" group$", "", g); g <- sub(" [0-9]+$", "", g)
okg <- !is.na(g) & nzchar(g) & !grepl("_ge$|^uncultured|^unidentified", g, ignore.case=TRUE)
G <- rowsum(t(st)[okg, , drop=FALSE], g[okg]); stopifnot(nrow(G) > 50)
d <- file.path("profiles_16s", coh); dir.create(d, recursive=TRUE, showWarnings=FALSE)
write.table(data.frame(genus=rownames(G), G, check.names=FALSE), file.path(d, paste0(coh, "_genus_counts.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
trk <- data.frame(run=runs, input=out[, 1], filtered=out[, 2], nonchim=rowSums(st), classified_genus=colSums(G), prop_classified=round(colSums(G) / pmax(rowSums(st), 1), 3))
write.table(trk, file.path(R, "read_tracking.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(st=st, tx=tx, truncLen=c(t1, t2), mode=mode, merge_rate=mr), file.path(R, "dada2.rds"))
cat("mode", mode, "| samples", nrow(st), "| ASVs", ncol(st), "| genera", nrow(G), "| median classified share", median(trk$prop_classified), "\n")
