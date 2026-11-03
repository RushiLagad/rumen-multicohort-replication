#!/usr/bin/env Rscript
# Amendment 7 comparison (Amendment 16): Table 7 features, Standard vs capped GTDB, same 52 Nelore samples.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); set.seed(20261102)
od <- file.path(root, "results", "misassignment"); dir.create(od, recursive=TRUE, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
rdm <- function(f) { x <- rd(f); rownames(x) <- x[[1]]; x <- x[, -1]; colnames(x) <- sub("_clean$", "", colnames(x)); as.matrix(x) }
S <- rdm(file.path(root, "profiles", "nelore", "nelore_genus_counts.tsv"))
G <- rdm(file.path(root, "profiles_gtdb", "nelore", "nelore_genus_counts.tsv"))
stopifnot(setequal(colnames(S), colnames(G)), ncol(S) == 52); G <- G[, colnames(S)]
RS <- sweep(S, 2, colSums(S), "/"); RG <- sweep(G, 2, colSums(G), "/")
tk <- rd(file.path(root, "plan", "target_keys.tsv")); tk <- tk[tk$set == "Nelore", ]
env <- c("Actinomadura", "Actinoplanes", "Anaeromyxobacter", "Bosea", "Caulobacter", "Cellulomonas", "Comamonas", "Gordonia", "Lysobacter",
  "Nonomuraea", "Pseudonocardia", "Rubrobacter", "Sinorhizobium", "Sorangium", "Sphingobium", "Brevibacterium", "Burkholderia",
  "Aliivibrio", "Idiomarina", "Marinomonas", "Photobacterium", "Psychromonas")
tk$group <- ifelse(tk$target %in% env, "environmental", ifelse(tk$site == "rumen" & grepl("^RME", tk$phenotype), "stable set (rumen RME)",
  ifelse(tk$site == "rumen", "rumen RFI, rumen taxa", "faeces, other")))
keys <- lapply(strsplit(tk$gtdb_keys, ";"), function(x) { x <- trimws(x); x[nzchar(x)] })
med <- function(M, rows) if (!length(rows)) NA else median(colSums(M[rows, , drop=FALSE]))
tk$std_median <- vapply(tk$target, function(t) if (t %in% rownames(RS)) median(RS[t, ]) else NA_real_, 0)
tk$gtdb_median <- vapply(keys, function(k) { r <- intersect(k, rownames(RG)); if (!length(r)) 0 else med(RG, r) }, 0)
tk$gtdb_keys_present <- vapply(keys, function(k) length(intersect(k, rownames(RG))), 0L)
tk$included <- tk$n_keys > 0 & !is.na(tk$std_median) & tk$std_median > 0
tk$reduction <- ifelse(tk$included, 1 - tk$gtdb_median / tk$std_median, NA)
bm <- function(x) { x <- x[!is.na(x)]; b <- replicate(2000, median(sample(x, replace=TRUE))); c(median(x), quantile(b, c(.025, .975))) }
g <- do.call(rbind, lapply(split(tk, tk$group), function(d) { b <- bm(d$reduction)
  data.frame(group=d$group[1], n_in_family=nrow(d), n_included=sum(d$included), median_reduction=b[1], ci_lo=b[2], ci_hi=b[3], absent_in_gtdb=sum(d$included & d$gtdb_median == 0)) }))
write.table(tk[, c("id", "target", "site", "phenotype", "group", "n_keys", "gtdb_keys_present", "std_median", "gtdb_median", "included", "reduction")],
  file.path(od, "standard_vs_gtdb.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(g, file.path(od, "standard_vs_gtdb_groups.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(g, digits=3, row.names=FALSE)
E <- g$median_reduction[g$group == "environmental"]; Sb <- g$median_reduction[g$group == "stable set (rumen RME)"]
cat(sprintf("environmental median reduction %.0f%%, stable set %.0f%%, difference %.0f points\n", 100*E, 100*Sb, 100*(E - Sb)))
cat(if (E > 0.5 && Sb < 0.25) "PRE-SPECIFIED CRITERION MET: supports misassignment\n" else "PRE-SPECIFIED CRITERION NOT MET: does not support misassignment\n")
cat("excluded (no Standard counts or no GTDB counterpart):", paste(tk$target[!tk$included], collapse=", "), "\n")
