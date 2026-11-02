#!/usr/bin/env Rscript
# Level 1 report (Amendment 15): classify the Nelore set in Nelore. Usage: 24_reproduce_report.R
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "reproducibility")
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("run both fits first (RFI and RME)" = length(list.files(od, "_targets[.]tsv$")) == 2)
tk <- rd(file.path(root, "plan", "target_keys.tsv")); tk <- tk[tk$set == "Nelore", ]
est <- do.call(rbind, lapply(list.files(od, "_targets[.]tsv$", full.names=TRUE), rd)); est <- est[est$status == "ok", ]
pdir <- function(x) ifelse(grepl("lower with higher", x), -1, ifelse(grepl("higher with higher", x), 1, NA))
d <- cbind(tk[, c("id", "target", "site", "phenotype", "coef")], pub_dir=pdir(tk$direction))
d$stratum <- paste(ifelse(d$site == "feces", "cross-site", "rumen"), sub("_GC", "", d$phenotype))
j <- match(d$id, est$id); d$lfc <- est$lfc[j]; d$W <- est$W[j]; d$df <- est$df[j]
d$testable <- !is.na(d$W); d$r <- d$W / sqrt(d$W^2 + d$df); d$p <- ifelse(d$testable, 2 * pt(-abs(d$W), d$df), 1)
d$q_bh <- p.adjust(d$p, "BH"); d$q_testable <- NA; d$q_testable[d$testable] <- p.adjust(d$p[d$testable], "BH")
d$same_sign <- sign(d$lfc) == d$pub_dir
d$label <- ifelse(!d$testable, "untestable", ifelse(d$q_bh >= 0.05, "not reproduced", ifelse(d$same_sign, "reproduced", "opposite sign")))
d$reason <- ifelse(d$testable, "", ifelse(d$id %in% tk$id[tk$status != "mapped" | tk$n_keys == 0], "no genome-based counterpart", "below detection or filtered"))
write.table(d, file.path(od, "reproduction_all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("family size", nrow(d), "\n"); print(table(d$stratum, d$label))
t <- d[d$testable, ]; b <- binom.test(sum(t$same_sign), nrow(t))
cat(sprintf("same sign as published: %d of %d (%.0f%%, exact 95%% CI %.0f-%.0f%%)\n", sum(t$same_sign), nrow(t), 100*b$estimate, 100*b$conf.int[1], 100*b$conf.int[2]))
cs <- suppressWarnings(cor.test(t$coef, t$lfc, method="spearman")); cat(sprintf("Spearman published coef vs this estimate: rho %.2f (p %.3g), scales differ\n", cs$estimate, cs$p.value))
print(t[order(t$q_bh), c("id", "target", "stratum", "pub_dir", "lfc", "r", "q_bh", "q_testable", "label")][1:min(15, nrow(t)), ], digits=3, row.names=FALSE)
