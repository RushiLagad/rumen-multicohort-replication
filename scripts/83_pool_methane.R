#!/usr/bin/env Rscript
# Amendments 19, 32, 36, 37: calibrate per-table W by the permutation null, pool the methane targets over independent cohorts, label.
# Usage: 83_pool_methane.R <dryrun|real>. dryrun pools the default-seed shuffled fits (a code test on noise); real pools results/methane/.
suppressMessages(library(metafor)); a <- commandArgs(TRUE); mode <- a[1]; stopifnot(mode %in% c("dryrun", "real"))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
plan <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")); stopifnot("commit Amendment 37 first" = any(grepl("^\\*\\*Amendment 37 \\(", plan)))
ind <- file.path(root, "results", if (mode == "dryrun") "methane_dryrun" else "methane"); od <- file.path(root, "results", if (mode == "dryrun") "methane_pooled_dryrun" else "methane_pooled"); dir.create(od, showWarnings=FALSE)
nulldir <- file.path(root, "results", "methane_dryrun")
# ---- null SD of W per table and target from the seed-1001+ permutations
nullsd <- function(tab) { fs <- list.files(nulldir, paste0("^", tab, "_seed1[0-9]{3}_targets\\.tsv$"), full.names=TRUE)
  if (!length(fs)) return(NULL); x <- do.call(rbind, lapply(fs, function(f) { d <- rd(f); d <- d[d$status == "ok", c("id", "W")]; d$seed <- f; d }))
  s <- aggregate(W ~ id, x, function(w) c(sd=sd(w), k=length(w))); data.frame(id=s$id, sd_null=s$W[, "sd"], k_null=s$W[, "k"], W_null=I(split(x$W, x$id)[s$id])) }
tables <- c(difford_arch="difford", difford_bact="difford", wallace_prokaryotes="wallace")
est <- list()
for (tab in names(tables)) { f <- file.path(ind, paste0(tab, "_targets.tsv")); if (!file.exists(f)) { cat("missing", f, "\n"); next }
  e <- rd(f); e <- e[e$status == "ok", ]; if (tab == "difford_bact") e <- e[e$id != "M1", ]; if (tab == "difford_arch") e <- e[e$id == "M1", ]
  ns <- nullsd(tab); stopifnot("no permutation null for table" = !is.null(ns)); j <- match(e$id, ns$id)
  e$sd_null <- ns$sd_null[j]; e$k_null <- ns$k_null[j]; e$W_raw <- e$W
  e$p_perm <- mapply(function(w, i) if (is.na(i)) NA else mean(abs(ns$W_null[[i]]) >= abs(w)), e$W, j)
  e$W <- e$W_raw / pmax(1, e$sd_null, na.rm=TRUE); e$table <- tab; e$cohort <- tables[tab]; est[[tab]] <- e }
est <- do.call(rbind, est); stopifnot(!anyDuplicated(est[, c("id", "cohort")]))
est$reanalysis <- est$id == "M2" & est$cohort == "wallace"   # Amendment 32 item 1
est$r <- partial_r(est$W, est$df); z <- atanh(est$r); s <- sqrt(1/(est$df - 1)); est$lo <- tanh(z - 1.96*s); est$hi <- tanh(z + 1.96*s)
write.table(est, file.path(od, "estimates_calibrated.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
tk <- rd(file.path(root, "plan", "target_keys.tsv")); tk <- tk[tk$set %in% c("Methane", "MethaneSecondary") | (tk$set == "Nelore" & tk$phenotype == "RME_GC"), ]
pdir <- function(x) ifelse(grepl("lower with higher", x), -1, ifelse(grepl("higher with higher", x), 1, NA))
run <- function(e, label) { all <- list()
  for (st in unique(tk$set)) { d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pdir(d$direction), testable=d$status %in% c("mapped", "special"))
    x <- e[e$id %in% d$id, ]; res <- classify_set(tg, x[, c("id", "cohort", "W", "df", "n")])
    res <- cbind(set=st, target=d$target, published_dir=tg$dir, res)
    for (cc in c("difford", "wallace")) { y <- x[x$cohort == cc, ]; j <- match(d$id, y$id); res[[paste0("r_", cc)]] <- y$r[j]; res[[paste0("sd_null_", cc)]] <- y$sd_null[j]; res[[paste0("p_perm_", cc)]] <- y$p_perm[j] }
    all[[st]] <- res }
  all <- do.call(rbind, all); write.table(all, file.path(od, paste0("classification_", label, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE); all }
P <- run(est[!est$reanalysis, ], "primary"); S <- run(est, "with_reanalysis")
cat("== ", toupper(mode), ": primary pooling (Wallace M2 excluded as re-analysis)\n"); print(table(P$set, P$label))
x <- P[!startsWith(P$label, "untestable"), ]
print(x[order(x$set, x$q_bh), c("id", "target", "set", "k", "r", "ci_lo", "ci_hi", "q_bh", "r_difford", "r_wallace", "p_perm_difford", "p_perm_wallace", "label")], digits=2, row.names=FALSE)
cat("== calibration factors applied (sd_null > 1):\n"); print(est[est$sd_null > 1, c("id", "table", "sd_null", "W_raw", "W")], digits=2, row.names=FALSE)
