# Amendment 32 item 4: does each cohort's amplified region separate the M. gottschalkii and M. ruminantium clades? Reference sequences only; no study data.
suppressMessages(library(Biostrings)); setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"))
ref <- readDNAStringSet("raw/references/mbb_type_strains_16S.fasta")
names(ref) <- sub("^(\\S+) Methanobrevibacter (\\S+).*", "\\1_\\2", names(ref))
clade <- ifelse(grepl("gottschalkii|millerae|thaueri", names(ref)), "gott", "rum")
regions <- list(difford_arch_V4V6 = c("CAGCMGCCGCGGTAA", "GGCCATGCACCWCCTCTC"),       # published primers (Difford 2018)
                wallace_prok_V4_assumed = c("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT")) # ProkA primers unpublished; standard 515F/806R assumed
ext <- function(s, f, r) { fm <- matchPattern(DNAString(f), s, max.mismatch = 2, fixed = FALSE)
  rm <- matchPattern(reverseComplement(DNAString(r)), s, max.mismatch = 2, fixed = FALSE)
  if (!length(fm) || !length(rm)) return(NA_character_); as.character(subseq(s, end(fm)[1] + 1, start(rm)[1] - 1)) }
dif <- function(a, b) { al <- pairwiseAlignment(a, b); nmismatch(al) + sum(unlist(width(indel(pattern(al))))) + sum(unlist(width(indel(subject(al))))) }
for (rg in names(regions)) {
  x <- vapply(seq_along(ref), function(i) ext(ref[[i]], regions[[rg]][1], regions[[rg]][2]), "")
  cat("\n==", rg, "\n  insert lengths:", paste0(names(ref), "=", nchar(x), collapse = "; "), "\n")
  ok <- !is.na(x); if (!any(ok & clade == "gott") || !any(ok & clade == "rum")) { cat("  NOT COVERED by references of both clades\n"); next }
  idx <- which(ok); d <- outer(idx, idx, Vectorize(function(i, j) if (i == j) 0 else dif(x[i], x[j]))); dimnames(d) <- list(names(ref)[idx], names(ref)[idx]); print(d)
  g <- clade[idx] == "gott"
  cat("  min differences BETWEEN clades:", min(d[g, !g]), "| max WITHIN gottschalkii clade:", max(d[g, g]), "| max WITHIN ruminantium clade:", max(d[!g, !g]), "\n")
  cat("  verdict:", if (min(d[g, !g]) >= 1) "region separates the clades" else "region does NOT separate the clades -> M1 untestable", "\n") }
