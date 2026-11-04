#!/usr/bin/env Rscript
# Amendments 21, 22, 28, 29, 30: DADA2 for one Difford marker -> genus counts per cow. Error models per flowcell; runs summed per cow before chimera removal.
# Bacteria: merged table (primary unless median merge rate < 0.5) plus forward-read table (sensitivity). Usage: 73_dada2_difford.R <bact|arch> <insert bp> (env NT). No phenotype data loaded.
suppressMessages({library(dada2); library(ShortRead)})
a <- commandArgs(TRUE); M <- a[1]; INSERT <- as.numeric(a[2]); coh <- paste0("difford_", M)
setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")); set.seed(20261104); nt <- as.integer(Sys.getenv("NT", "16"))
O <- file.path("work/16s", coh); R <- file.path("results/16s", coh); dir.create(R, recursive=TRUE, showWarnings=FALSE)
map <- read.delim("harmonized/difford_runs.tsv", colClasses="character"); map <- map[map$in_s4 == "True" & map$marker == M, ]
runs <- readLines(file.path(O, "runs.txt")); f1 <- file.path(O, "trimmed_sub", paste0(runs, "_1.fastq.gz")); f2 <- sub("_1\\.fastq\\.gz$", "_2.fastq.gz", f1)
ok <- file.exists(f1) & file.exists(f2) & file.size(f1) > 0; runs <- runs[ok]; f1 <- f1[ok]; f2 <- f2[ok]
fc <- map$flowcell[match(runs, map$run)]; cow <- map$cow[match(runs, map$run)]; stopifnot(!anyNA(fc), !anyNA(cow))
cat("runs with trimmed pairs:", length(runs), "| cows:", length(unique(cow)), "| flowcells:", length(unique(fc)), "\n")
medq <- function(fs, W = 320, k = 5000) apply(do.call(rbind, lapply(fs, function(f) { s <- FastqSampler(f, k); q <- as(quality(yield(s)), "matrix"); close(s)
  q <- q[, seq_len(min(ncol(q), W)), drop=FALSE]; if (ncol(q) < W) q <- cbind(q, matrix(NA_real_, nrow(q), W - ncol(q))); q })), 2, median, na.rm=TRUE)
last <- function(q, cap, thr) { i <- which(q[seq_len(cap)] >= thr); if (length(i)) max(i) else NA }
q1 <- medq(f1); q2 <- medq(f2); t1 <- last(q1, 280, 30); t2 <- last(q2, 260, 25)
while (!is.na(t1) && !is.na(t2) && t1 + t2 - INSERT < 40 && t2 < 260) t2 <- t2 + 1
write.table(data.frame(pos=seq_along(q1), R1=q1, R2=q2), file.path(R, "median_quality.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("truncLen R1", t1, "R2", t2, "expected overlap", t1 + t2 - INSERT, "\n"); stopifnot(!is.na(t1), !is.na(t2), t1 >= 150, t2 >= 100)
den <- function(a1, a2, fcs, rn) { res <- vector("list", length(a1))
  for (f in unique(fcs)) { i <- which(fcs == f)
    e1 <- learnErrors(a1[i], multithread=nt, randomize=TRUE, verbose=0); d1 <- dada(a1[i], err=e1, multithread=nt, verbose=0); if (length(i) == 1) d1 <- list(d1)
    if (is.null(a2)) x <- d1 else {
      e2 <- learnErrors(a2[i], multithread=nt, randomize=TRUE, verbose=0); d2 <- dada(a2[i], err=e2, multithread=nt, verbose=0); if (length(i) == 1) d2 <- list(d2)
      x <- mergePairs(d1, a1[i], d2, a2[i], minOverlap=20, maxMismatch=0); if (is.data.frame(x)) x <- list(x) }
    res[i] <- x; cat("flowcell", f, "| runs", length(i), "\n") }
  names(res) <- rn; res }
fwd <- function() { tf <- last(q1, 220, 30); cat("forward reads: truncLen", tf, "\n"); stopifnot(!is.na(tf), tf >= 150)
  dir.create(file.path(O, "filtered_fwd"), showWarnings=FALSE); h1 <- file.path(O, "filtered_fwd", basename(f1))
  o <- filterAndTrim(f1, h1, truncLen=tf, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=nt)
  k <- o[, 2] > 0 & file.exists(h1)
  list(mg=den(h1[k], NULL, fc[k], runs[k]), runs=runs[k], cow=cow[k], fc=fc[k], out=o[k, , drop=FALSE], tl=tf) }
finish <- function(mg, rn, cw, fcs, out, tl, mode, tag, mr) {
  st <- makeSequenceTable(mg); rownames(st) <- rn
  stc <- rowsum(st, cw); stc <- removeBimeraDenovo(stc, method="consensus", multithread=nt)
  tx <- assignTaxonomy(colnames(stc), "work/charolais16s/silva138_train.fa", minBoot=50, multithread=nt, taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
  if (M == "arch") { ka <- tx[, "Kingdom"] %in% "Archaea"; cat("off-target (non-archaeal) read share", round(1 - sum(stc[, ka]) / sum(stc), 3), "\n")
    stc <- stc[, ka, drop=FALSE]; tx <- tx[ka, , drop=FALSE] }
  g <- unname(tx[, "Genus"]); g <- gsub("[][]", "", g); g <- sub(" group$", "", g); g <- sub(" [0-9]+$", "", g)
  okg <- !is.na(g) & nzchar(g) & !grepl("_ge$|^uncultured|^unidentified", g, ignore.case=TRUE)
  G <- rowsum(t(stc)[okg, , drop=FALSE], g[okg]); stopifnot(nrow(G) >= if (M == "arch") 1 else 50)
  d <- file.path("profiles_16s", coh); dir.create(d, recursive=TRUE, showWarnings=FALSE)
  write.table(data.frame(genus=rownames(G), G, check.names=FALSE), file.path(d, paste0(coh, tag, "_genus_counts_unfiltered.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  nc <- rowSums(stc); lo <- names(nc)[nc < 1000]
  write.table(data.frame(cow=lo, nonchim=unname(nc[lo]), reason=rep("fewer than 1000 non-chimeric reads (Amendment 22)", length(lo))),
              file.path(R, paste0("excluded", tag, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  G2 <- G[, !colnames(G) %in% lo, drop=FALSE]
  write.table(data.frame(genus=rownames(G2), G2, check.names=FALSE), file.path(d, paste0(coh, tag, "_genus_counts.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(data.frame(run=rn, cow=cw, flowcell=fcs, input=out[, 1], filtered=out[, 2]), file.path(R, paste0("read_tracking_runs", tag, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  cg <- colSums(G)[rownames(stc)]
  write.table(data.frame(cow=rownames(stc), nonchim=nc, classified_genus=cg, prop_classified=round(cg / pmax(nc, 1), 3)),
              file.path(R, paste0("read_tracking_cows", tag, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  saveRDS(list(st=stc, tx=tx, truncLen=tl, mode=mode, merge_rate=mr), file.path(R, paste0("dada2", tag, ".rds")))
  cat("table", if (tag == "") "primary" else tag, "| mode", mode, "| cows", nrow(stc), "| excluded", length(lo), "| ASVs", ncol(stc), "| genera", nrow(G),
      "| median classified share", round(median(cg / pmax(nc, 1)), 3), "\n") }
dir.create(file.path(O, "filtered"), showWarnings=FALSE); g1 <- file.path(O, "filtered", basename(f1)); g2 <- file.path(O, "filtered", basename(f2))
out <- filterAndTrim(f1, g1, f2, g2, truncLen=c(t1, t2), maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=nt)
keep <- out[, 2] > 0 & file.exists(g1)
mg <- den(g1[keep], g2[keep], fc[keep], runs[keep]); mr <- vapply(mg, function(x) sum(x$abundance), 0) / out[keep, 2]
cat("median merge rate", round(median(mr), 3), "\n")
if (median(mr) >= 0.5) {
  finish(mg, runs[keep], cow[keep], fc[keep], out[keep, , drop=FALSE], c(t1, t2), "merged", "", mr)
  if (M == "bact") { f <- fwd(); finish(f$mg, f$runs, f$cow, f$fc, f$out, f$tl, "forward", "_fwd", NA) }
} else { cat("FALLBACK to forward reads (Amendment 21 item 4)\n"); f <- fwd(); finish(f$mg, f$runs, f$cow, f$fc, f$out, f$tl, "forward", "", mr) }
