# Amendment 35 (pending): genus -> family/phylum/kingdom from each 16S cohort's own SILVA assignments (genus names cleaned exactly as in the
# genus tables). A genus with ASVs in more than one SILVA family takes the family holding most of its reads. No phenotype data.
# Usage: 84_silva_lineage.R <cohort> <dada2.rds> <genus_counts.tsv>
setwd(Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort")); a <- commandArgs(TRUE); x <- readRDS(a[2]); tx <- x$tx; st <- x$st
cl <- function(v) { v <- gsub("[][]", "", v); v <- sub(" group$", "", v); sub(" [0-9]+$", "", v) }
g <- cl(unname(tx[, "Genus"])); ok <- !is.na(g) & nzchar(g) & !grepl("_ge$|^uncultured|^unidentified", g, ignore.case = TRUE)
rd <- if (nrow(st) == ncol(tx) || ncol(st) != nrow(tx)) rowSums(st) else colSums(st)   # reads per ASV (st is samples x ASVs)
rd <- colSums(st)[rownames(tx)]; if (anyNA(rd)) rd <- colSums(st)
d <- data.frame(genus = g, family = gsub("[][]", "", unname(tx[, "Family"])), phylum = unname(tx[, "Phylum"]), kingdom = unname(tx[, "Kingdom"]), reads = unname(rd))[ok, ]
d$family[is.na(d$family) | !nzchar(d$family)] <- "unassigned_family"
agg <- aggregate(reads ~ genus + family + phylum + kingdom, d, sum); agg <- agg[order(agg$genus, -agg$reads), ]
multi <- unique(agg$genus[duplicated(agg$genus)]); m <- agg[!duplicated(agg$genus), c("genus", "family", "phylum", "kingdom")]
gt <- read.delim(a[3], check.names = FALSE)[[1]]; miss <- setdiff(gt, m$genus)
dir.create(file.path("results/16s", a[1]), showWarnings = FALSE, recursive = TRUE)
write.table(m, file.path("results/16s", a[1], "genus_lineage_silva.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(a[1], "| genera in table", length(gt), "| mapped", length(gt) - length(miss), "| genera split across SILVA families (majority used):", length(multi), "\n")
