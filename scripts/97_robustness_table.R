#!/usr/bin/env Rscript
# Confirmatory replication labels across every analysis version, both arms. Reads existing classification files only. Usage: 97_robustness_table.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "robustness"); dir.create(od, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
fs <- c(Sys.glob(file.path(root, "results", "replication_rfi3*", "classification_*.tsv")), Sys.glob(file.path(root, "results", "methane_pooled*", "classification_*.tsv")), Sys.glob(file.path(root, "results", "clr_pooled", "classification_*.tsv")))
fs <- fs[!grepl("dryrun|superseded", fs)]
L <- do.call(rbind, lapply(fs, function(f) { d <- rd(f); v <- paste0(basename(dirname(f)), ":", sub("^classification_?", "", sub("[.]tsv$", "", basename(f))))
  arm <- if (grepl("methane", f)) "Methane" else "RFI"
  d <- if (arm == "RFI") d[!d$set %in% c("Methane", "MethaneSecondary"), ] else d[d$set %in% c("Methane", "MethaneSecondary", "Nelore"), ]
  data.frame(arm=arm, version=v, set=d$set, id=d$id, target=d$target, k=d$k, r=round(d$r, 3), q=round(d$q_bh, 3), label=d$label, stringsAsFactors=FALSE) }))
L <- L[!startsWith(L$label, "untestable"), ]; write.table(L, file.path(od, "labels_long.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== versions found\n"); print(table(L$version, L$arm))
cat("\n== label counts per version (testable targets)\n"); S <- as.data.frame.matrix(table(paste(L$arm, L$version, sep=" | "), sub(" [(]direction not stated[)]", "", L$label))); print(S)
write.table(cbind(version=rownames(S), S), file.path(od, "summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
for (a in c("RFI", "Methane")) { x <- L[L$arm == a, ]; W <- reshape(x[, c("id", "target", "version", "label")], idvar=c("id", "target"), timevar="version", direction="wide"); names(W) <- sub("^label[.]", "", names(W))
  write.table(W, file.path(od, paste0("labels_wide_", a, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  lab <- as.matrix(W[, -(1:2), drop=FALSE]); lab <- gsub(" [(]direction not stated[)]", "", lab); dif <- apply(lab, 1, function(z) length(unique(z[!is.na(z)])) > 1)
  cat("\n==", a, ": targets", nrow(W), "| versions", ncol(lab), "| targets whose label differs in any version:", sum(dif), "\n")
  if (any(dif)) print(W[dif, ], row.names=FALSE)
  rr <- reshape(x[, c("id", "version", "r")], idvar="id", timevar="version", direction="wide"); m <- as.matrix(rr[, -1]); sg <- apply(m, 1, function(z) { z <- z[is.finite(z) & abs(z) > 0]; length(unique(sign(z))) > 1 })
  cat("   targets whose pooled r changes sign across versions:", sum(sg), "of", nrow(rr), if (any(sg)) paste0("(", paste(rr$id[sg], collapse=", "), ")") else "", "\n") }
