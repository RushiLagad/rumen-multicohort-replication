#!/usr/bin/env Rscript
# Amendments 22, 28 (items 6-8): DADA2 for one Wallace marker -> genus counts per cow. Reads are pre-merged and primer-free (Amendment 33: real merged qualities, maxEE 2); length window =
# census median +/- 10 bp; error models per sequencing round; abundance-only error model if qualities are flat; runs denoised one at a time;
# SR1 + SR2 summed per cow before chimera removal. Usage: 79_dada2_wallace.R <Prokaryotes|Bacteria|Archaea> (env NT, LIMIT). No phenotype data.
suppressMessages({library(dada2); library(ShortRead)})
M <- commandArgs(TRUE)[1]; coh <- paste0("wallace_", tolower(M)); setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")); set.seed(20261105)
nt <- as.integer(Sys.getenv("NT", "16")); O <- file.path("work/16s", coh); R <- file.path("results/16s", coh)
dir.create(file.path(O, "filtered"), recursive = TRUE, showWarnings = FALSE); dir.create(R, recursive = TRUE, showWarnings = FALSE)
map <- read.delim("harmonized/wallace_runs.tsv", colClasses = "character"); map <- map[map$marker == M, ]; map <- map[order(map$run), ]
if (nzchar(Sys.getenv("LIMIT"))) map <- head(map, as.integer(Sys.getenv("LIMIT")))
stopifnot(all(file.exists(map$file))); med <- as.numeric(map$marker_median_len[1]); lo <- med - 10; hi <- med + 10
qv <- unique(unlist(lapply(head(map$file, 20), function(x) { s <- FastqSampler(x, 2000); q <- as(quality(yield(s)), "matrix"); close(s); unique(as.vector(q)) })))
qv <- qv[!is.na(qv)]; flat <- length(qv) <= 2; ef <- if (flat) noqualErrfun else loessErrfun
cat("marker", M, "| runs", nrow(map), "| cows", length(unique(map$cow)), "| length window", lo, "-", hi, "| distinct quality values", length(qv),
    if (flat) "-> abundance-only error model (noqualErrfun)" else "-> standard error model (loessErrfun)", "\n")
g <- file.path(O, "filtered", paste0(map$run, ".fastq.gz"))
out <- filterAndTrim(map$file, g, minLen = lo, maxLen = hi, maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = 8)
keep <- out[, 2] > 0 & file.exists(g); map <- map[keep, ]; g <- g[keep]; out <- out[keep, , drop = FALSE]
res <- vector("list", nrow(map))
for (rd in unique(map$round)) { i <- which(map$round == rd)
  e <- learnErrors(g[i], multithread = nt, randomize = TRUE, errorEstimationFunction = ef, verbose = 0)
  for (j in i) res[[j]] <- getUniques(dada(g[j], err = e, multithread = nt, verbose = 0))
  invisible(gc()); cat("round", rd, "| runs", length(i), "|", format(Sys.time(), "%H:%M"), "\n") }
names(res) <- map$run; st <- makeSequenceTable(res); rownames(st) <- map$run
stc <- rowsum(st, map$cow); stc <- removeBimeraDenovo(stc, method = "consensus", multithread = nt)
tx <- assignTaxonomy(colnames(stc), "work/charolais16s/silva138_train.fa", minBoot = 50, multithread = nt, taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
if (M == "Archaea") { ka <- tx[, "Kingdom"] %in% "Archaea"; cat("off-target (non-archaeal) read share", round(1 - sum(stc[, ka]) / sum(stc), 3), "\n"); stc <- stc[, ka, drop = FALSE]; tx <- tx[ka, , drop = FALSE] }
gn <- unname(tx[, "Genus"]); gn <- gsub("[][]", "", gn); gn <- sub(" group$", "", gn); gn <- sub(" [0-9]+$", "", gn)
okg <- !is.na(gn) & nzchar(gn) & !grepl("_ge$|^uncultured|^unidentified", gn, ignore.case = TRUE)
G <- rowsum(t(stc)[okg, , drop = FALSE], gn[okg]); stopifnot(nrow(G) >= if (M == "Archaea") 1 else 50)
d <- file.path("profiles_16s", coh); dir.create(d, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(genus = rownames(G), G, check.names = FALSE), file.path(d, paste0(coh, "_genus_counts_unfiltered.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
nc <- rowSums(stc); lw <- names(nc)[nc < 1000]
write.table(data.frame(cow = lw, nonchim = unname(nc[lw]), reason = rep("fewer than 1000 non-chimeric reads (Amendment 22)", length(lw))), file.path(R, "excluded.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
G2 <- G[, !colnames(G) %in% lw, drop = FALSE]
write.table(data.frame(genus = rownames(G2), G2, check.names = FALSE), file.path(d, paste0(coh, "_genus_counts.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(run = map$run, cow = map$cow, round = map$round, input = out[, 1], filtered = out[, 2]), file.path(R, "read_tracking_runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cg <- colSums(G)[rownames(stc)]
write.table(data.frame(cow = rownames(stc), nonchim = nc, classified_genus = cg, prop_classified = round(cg / pmax(nc, 1), 3)), file.path(R, "read_tracking_cows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(st = stc, tx = tx, length_window = c(lo, hi), error_model = if (flat) "noqualErrfun" else "loessErrfun"), file.path(R, "dada2.rds"))
cat("table | marker", M, "| cows", nrow(stc), "| excluded", length(lw), "| ASVs", ncol(stc), "| genera", nrow(G), "| median classified share", round(median(cg / pmax(nc, 1)), 3), "\n")
