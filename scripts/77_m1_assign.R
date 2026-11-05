# Amendment 32 item 4: assign Methanobrevibacter ASVs to the M. gottschalkii or M. ruminantium clade by mismatches to RefSeq type strains
# over the cohort's amplified region; report the share of Methanobrevibacter reads assigned (M1 testable only if >= 70%). No phenotype data.
# Usage: 77_m1_assign.R <cohort> <dada2.rds> <fwd primer> <rev primer>
suppressMessages(library(Biostrings)); setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"))
a <- commandArgs(TRUE); coh <- a[1]; x <- readRDS(a[2]); st <- x$st; tx <- x$tx
ref <- readDNAStringSet("raw/references/mbb_type_strains_16S.fasta"); names(ref) <- sub("^(\\S+) Methanobrevibacter (\\S+).*", "\\1_\\2", names(ref))
rc <- ifelse(grepl("gottschalkii|millerae|thaueri", names(ref)), "gott", "rum")
ext <- function(s, f, r) { fm <- matchPattern(DNAString(f), s, max.mismatch = 2, fixed = FALSE)
  rm <- matchPattern(reverseComplement(DNAString(r)), s, max.mismatch = 2, fixed = FALSE)
  if (!length(fm) || !length(rm)) return(NA_character_); as.character(subseq(s, end(fm)[1] + 1, start(rm)[1] - 1)) }
dif <- function(a, b) if (nchar(a) == nchar(b)) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]]) else {
  al <- pwalign::pairwiseAlignment(a, b); pwalign::nmismatch(al) + sum(unlist(width(pwalign::indel(pwalign::pattern(al))))) + sum(unlist(width(pwalign::indel(pwalign::subject(al))))) }
ri <- vapply(seq_along(ref), function(i) ext(ref[[i]], a[3], a[4]), ""); ok <- !is.na(ri); ri <- ri[ok]; rc <- rc[ok]
mb <- which(tx[, "Genus"] %in% "Methanobrevibacter"); seqs <- colnames(st)[mb]
D <- vapply(seqs, function(s) vapply(ri, function(r) dif(s, r), 0), numeric(length(ri)))
bg <- apply(D[rc == "gott", , drop = FALSE], 2, min); br <- apply(D[rc == "rum", , drop = FALSE], 2, min)
cl <- ifelse(bg < br, "gottschalkii", ifelse(br < bg, "ruminantium", "unassigned"))
rd <- colSums(st[, mb, drop = FALSE]); share <- sum(rd[cl != "unassigned"]) / sum(rd)
dir.create(file.path("results/16s", coh), showWarnings = FALSE, recursive = TRUE)
write.table(data.frame(asv = seq_along(seqs), clade = cl, best_gott = bg, best_rum = br, reads = rd, sequence = seqs),
            file.path("results/16s", coh, "m1_asv_clades.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cnt <- sapply(c("gottschalkii", "ruminantium", "unassigned"), function(k) rowSums(st[, mb[cl == k], drop = FALSE]))
write.table(data.frame(cow = rownames(st), cnt, mbb_total = rowSums(cnt), check.names = FALSE),
            file.path("profiles_16s", coh, paste0(coh, "_m1_counts.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(coh, "| Methanobrevibacter ASVs:", length(seqs), "| share of archaeal reads:", round(sum(rd) / sum(st), 3), "\n")
cat("  ASVs by clade:", paste(names(table(cl)), table(cl), collapse = "; "), "\n")
cat("  reads by clade:", paste(names(rowsum(rd, cl)[, 1]), round(rowsum(rd, cl)[, 1] / sum(rd), 3), collapse = "; "), "\n")
cat("  best-match distance to nearest reference (quantiles 0/25/50/75/100%):", quantile(pmin(bg, br)), "\n")
cat("  share of Methanobrevibacter reads assigned:", round(share, 3), "->", if (share >= 0.7) "M1 TESTABLE (>= 0.70)" else "M1 UNTESTABLE (< 0.70)", "\n")
