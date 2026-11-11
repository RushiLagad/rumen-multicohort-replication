#!/usr/bin/env Rscript
# Amendment 44 items 2 to 4: pool the late LinDA or ALDEx2 estimates exactly as the calibrated primary (89 for RFI, Alberta all-animal; 83 for methane),
# after calibrating each table by its within-design permutations (W_cal = W / max(1, null SD), p_perm = share of |null W| >= |W|; Amendment 37),
# and compare with the calibrated ANCOM-BC2 primary. Late analysis, run after all primary results were seen. Usage: 106_late_pool.R <linda|aldex2>
suppressMessages(library(metafor)); M <- commandArgs(TRUE)[1]; stopifnot(M %in% c("linda", "aldex2"))
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
LATE_LABEL <- "late analysis, run after all primary results were seen"
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("commit Amendment 44 first" = any(grepl("^\\*\\*Amendment 44 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))))
ld <- file.path(root, "results", paste0("late_", M)); od <- file.path(ld, "pooled"); dir.create(od, showWarnings=FALSE)
sink(file.path(od, "summary.txt"), split=TRUE); cat(toupper(M), "|", LATE_LABEL, "| Amendment 44 |", format(Sys.time()), "\n")
RFI_T <- c(paste0(rep(c("nelore", "charolais_mrna"), each=4), "_", c("species", "genus", "family", "phylum")), "alberta_genus_all", "alberta_family_all")
CH4_T <- c("difford_bact", "wallace_prokaryotes")
nulls <- function(tab) { fs <- list.files(file.path(ld, "null", tab), "^seed[0-9]+\\.tsv$", full.names=TRUE); stopifnot("no permutations for table" = length(fs) > 0)
  x <- do.call(rbind, lapply(fs, rd)); x <- x[is.finite(x$W), ]; list(K=length(fs), W=split(x$W, x$id)) }
calibrate <- function(tab) { e <- rd(file.path(ld, paste0(tab, "_targets.tsv"))); e <- e[e$status == "ok", ]; if (!nrow(e)) return(NULL); ns <- nulls(tab)
  w0 <- ns$W[e$id]; k <- lengths(w0); stopifnot("target missing from null" = all(k > 0), "fewer than 95% of permutations usable" = all(k >= 0.95 * ns$K))
  e$table <- tab; e$sd_null <- vapply(w0, sd, 0); e$k_null <- k; e$W_raw <- e$W
  e$p_perm <- mapply(function(w, z) mean(abs(z) >= abs(w)), e$W_raw, w0); e$W <- e$W_raw / pmax(1, e$sd_null)
  e$r <- partial_r(e$W, e$df); z <- atanh(e$r); s <- sqrt(1 / (e$df - 1)); e$lo <- tanh(z - 1.96 * s); e$hi <- tanh(z + 1.96 * s); e$analysis <- LATE_LABEL; e }
tk0 <- rd(file.path(root, "plan", "target_keys.tsv")); cols <- c("id", "cohort", "W", "df", "n")
sgn <- function(a, b) { ok <- is.finite(a) & is.finite(b) & a != 0 & b != 0; c(same=sum(sign(a[ok]) == sign(b[ok])), of=sum(ok)) }
report <- function(arm, cls, old, est, oldest) {
  j <- match(cls$id, old$id); te <- !startsWith(cls$label, "untestable") | !startsWith(old$label[j], "untestable")
  chg <- te & cls$label != old$label[j]
  cat("\n==", arm, "|", M, "calibrated | labels:\n"); print(table(cls$set, cls$label))
  cat("replicated:", sum(cls$label == "replicated"), "| labels differing from the calibrated ANCOM-BC2 primary:", sum(chg), "\n")
  if (any(chg)) print(data.frame(id=cls$id, target=cls$target, primary=old$label[j], late=cls$label)[chg, ], row.names=FALSE)
  ps <- sgn(cls$r[te], old$r[j][te]); m <- merge(est[, c("id", "cohort", "W")], oldest[, c("id", "cohort", "W")], by=c("id", "cohort"), suffixes=c("_late", "_ancom")); cs <- sgn(m$W_late, m$W_ancom)
  cat("sign agreement with ANCOM-BC2: pooled r", ps[1], "of", ps[2], "| cohort-level W", cs[1], "of", cs[2], "| Spearman cohort-level W", round(cor(m$W_late, m$W_ancom, method="spearman"), 2), "\n")
  x <- data.frame(id=cls$id, target=cls$target, k=cls$k, r_ancom=round(old$r[j], 3), r_late=round(cls$r, 3), ci_late=paste(round(cls$ci_lo, 2), "to", round(cls$ci_hi, 2)),
    q_ancom=round(old$q_bh[j], 2), q_late=round(cls$q_bh, 2), label_late=cls$label)[te, ]; print(x[order(x$q_late), ], row.names=FALSE)
  data.frame(method=M, arm=arm, testable=sum(!startsWith(cls$label, "untestable")), replicated=sum(cls$label == "replicated"), labels_changed=sum(chg),
    pooled_sign_same=ps[1], pooled_sign_of=ps[2], cohort_sign_same=cs[1], cohort_sign_of=cs[2]) }
# ---- RFI (as 89, Alberta all-animal)
pR <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
tk <- tk0[!tk0$set %in% c("Methane", "MethaneSecondary"), ]
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available", ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
E <- do.call(rbind, lapply(RFI_T, calibrate)); stopifnot("a target appears twice for one cohort" = !anyDuplicated(E[, c("id", "cohort")]))
R <- do.call(rbind, lapply(unique(tk$set), function(st) { d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pR(d$direction), testable=d$stratum %in% c("replication", "cross-site"))
  cbind(set=st, target=d$target, rank=d$target_rank, stratum=d$stratum, published_dir=tg$dir, classify_set(tg, E[E$id %in% d$id, cols])) }))
R$analysis <- LATE_LABEL; write.table(E, file.path(od, "estimates_calibrated_rfi.tsv"), sep="\t", quote=FALSE, row.names=FALSE); write.table(R, file.path(od, "classification_rfi.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
S1 <- report("RFI", R, rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")), E, rd(file.path(root, "results", "replication_rfi3_calibrated", "estimates_calibrated_all.tsv")))
# ---- methane (as 83; M1 is the primary calibrated linear model in every version; Wallace M2 is re-analysis, Amendment 32)
pM <- function(x) ifelse(grepl("lower with higher", x), -1, ifelse(grepl("higher with higher", x), 1, NA))
tm <- tk0[tk0$set %in% c("Methane", "MethaneSecondary") | (tk0$set == "Nelore" & tk0$phenotype == "RME_GC"), ]
Em <- do.call(rbind, lapply(CH4_T, calibrate)); Em$cohort <- ifelse(startsWith(Em$cohort, "difford"), "difford", "wallace")
pm <- rd(file.path(root, "results", "methane_pooled", "estimates_calibrated.tsv")); m1 <- pm[pm$id == "M1", ]; m1$method <- "lm (primary M1)"; m1$analysis <- "primary M1 linear model, calibrated (Amendment 37)"
for (cc in setdiff(names(Em), names(m1))) m1[[cc]] <- NA; Em <- rbind(Em, m1[, names(Em)]); stopifnot(!anyDuplicated(Em[, c("id", "cohort")]))
Em$reanalysis <- Em$id == "M2" & Em$cohort == "wallace"
runm <- function(e, lab) { x <- do.call(rbind, lapply(unique(tm$set), function(st) { d <- tm[tm$set == st, ]; tg <- data.frame(id=d$id, dir=pM(d$direction), testable=d$status %in% c("mapped", "special"))
  cbind(set=st, target=d$target, published_dir=tg$dir, classify_set(tg, e[e$id %in% d$id, cols])) })); x$analysis <- LATE_LABEL
  write.table(x, file.path(od, paste0("classification_methane_", lab, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE); x }
Mp <- runm(Em[!Em$reanalysis, ], "primary"); Ms <- runm(Em, "with_reanalysis"); write.table(Em, file.path(od, "estimates_calibrated_methane.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
pe <- pm[pm$id != "M1", ]; S2 <- report("Methane", Mp, rd(file.path(root, "results", "methane_pooled", "classification_primary.tsv")), Em[Em$id != "M1", ], pe)
cat("methane with re-analysis (Wallace M2 included): replicated", sum(Ms$label == "replicated"), "\n")
# ---- null SD ranges and permutation counts
A <- rbind(E[, c("table", "sd_null", "k_null")], Em[Em$method != "lm (primary M1)", c("table", "sd_null", "k_null")])
ns <- do.call(rbind, lapply(split(A, A$table), function(d) data.frame(table=d$table[1], targets=nrow(d), permutations=max(d$k_null), sd_min=round(min(d$sd_null), 2), sd_median=round(median(d$sd_null), 2), sd_max=round(max(d$sd_null), 2))))
cat("\n== null SD of W by table (calibration applies where > 1)\n"); print(ns, row.names=FALSE)
write.table(ns, file.path(od, "null_sd_by_table.tsv"), sep="\t", quote=FALSE, row.names=FALSE); writeLines(as.character(max(ns$permutations)), file.path(od, "n_permutations.txt"))
S <- rbind(S1, S2); write.table(S, file.path(od, "summary_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE); cat("\n"); print(S, row.names=FALSE); sink()
