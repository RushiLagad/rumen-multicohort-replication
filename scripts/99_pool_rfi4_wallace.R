#!/usr/bin/env Rscript
# Amendment 43 item 4: calibrated RFI pooling as 89 (Alberta all-animal), with Wallace RFI NRC added as a fourth cohort. Sensitivity only.
suppressMessages(library(metafor)); root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "replication_rfi4_wallace_sensitivity"); dir.create(od, showWarnings=FALSE)
sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])); source(file.path(sdir, "lib_replication.R"))
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot("commit Amendment 43 first" = any(grepl("^\\*\\*Amendment 43 \\(", readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")))), file.exists(file.path(root, "logs", "wallace_rfi", "DONE")))
tk <- rd(file.path(root, "plan", "target_keys.tsv")); tk <- tk[!tk$set %in% c("Methane", "MethaneSecondary"), ]
pdir <- function(x) ifelse(grepl("higher in L-RFI|lower with higher", x), -1, ifelse(grepl("higher in H-RFI|higher with higher", x), 1, NA))
tk$stratum <- ifelse(!grepl("^RFI", tk$phenotype), "untestable: phenotype not available", ifelse(tk$status != "mapped" | tk$n_keys == 0, "untestable: no genome-based counterpart", ifelse(tk$site == "feces", "cross-site", "replication")))
calibrate <- function(f, nulldir, tab, cohort=NULL) { e <- rd(f); e <- e[e$status == "ok", ]; if (!nrow(e)) return(NULL)
  fs <- list.files(nulldir, paste0("^", tab, "_seed1[0-9]{3}_targets[.]tsv$"), full.names=TRUE); stopifnot("no null" = length(fs) >= 190)
  x <- do.call(rbind, lapply(fs, function(g) { d <- rd(g); d[d$status == "ok", c("id", "W")] })); sp <- split(x$W, x$id)
  stopifnot("target missing from null" = all(e$id %in% names(sp))); e$table <- tab; e$W_raw <- e$W; e$sd_null <- sapply(sp[e$id], sd); e$k_null <- lengths(sp[e$id])
  e$p_perm <- mapply(function(w, i) (sum(abs(sp[[i]]) >= abs(w)) + 1) / (length(sp[[i]]) + 1), e$W, e$id); e$W <- e$W_raw / pmax(1, e$sd_null)
  if (!is.null(cohort)) e$cohort <- cohort; e }
base <- list.files(file.path(root, "results", "replication"), "_targets[.]tsv$", full.names=TRUE); stopifnot(length(base) == 8)
cols <- c("id", "cohort", "table", "W_raw", "sd_null", "p_perm", "W", "df", "n")
E <- do.call(rbind, c(lapply(base, function(f) calibrate(f, file.path(root, "results", "rfi_null"), sub("_targets[.]tsv$", "", basename(f)))[, cols]),
  lapply(c("genus", "family"), function(rk) calibrate(file.path(root, "results", "alberta", paste0("alberta_", rk, "_all_targets.tsv")), file.path(root, "results", "rfi_null"), paste0("alberta_", rk, "_all"))[, cols]),
  list(calibrate(file.path(root, "results", "methane", "wallace_prokaryotes_rfi_nrc_targets.tsv"), file.path(root, "results", "methane_dryrun"), "wallace_prokaryotes_rfi_nrc", cohort="wallace")[, cols])))
stopifnot(!anyDuplicated(E[, c("id", "cohort")])); E$r <- partial_r(E$W, E$df); write.table(E, file.path(od, "estimates_calibrated.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
R <- do.call(rbind, lapply(unique(tk$set), function(st) { d <- tk[tk$set == st, ]; tg <- data.frame(id=d$id, dir=pdir(d$direction), testable=d$stratum %in% c("replication", "cross-site"))
  res <- cbind(set=st, target=d$target, published_dir=tg$dir, classify_set(tg, E[E$id %in% d$id, c("id", "cohort", "W", "df", "n")]))
  for (cc in unique(E$cohort)) { x <- E[E$cohort == cc, ]; j <- match(d$id, x$id); res[[paste0("r_", cc)]] <- x$r[j] }; res }))
write.table(R, file.path(od, "classification.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("== Wallace RFI NRC, calibrated (sensitivity; Amendment 43)\n"); w <- E[E$cohort == "wallace", ]; w$target <- tk$target[match(w$id, tk$id)]
print(w[, c("id", "target", "n", "W_raw", "sd_null", "W", "p_perm", "r")], digits=3, row.names=FALSE)
old <- rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")); j <- match(R$id, old$id)
x <- data.frame(id=R$id, target=R$target, k3=old$k[j], k4=R$k, r3=round(old$r[j], 3), r4=round(R$r, 3), ci4=paste0(round(R$ci_lo, 2), " to ", round(R$ci_hi, 2)), q3=round(old$q_bh[j], 2), q4=round(R$q_bh, 2), label3=old$label[j], label4=R$label)[R$k > old$k[j], ]
cat("\n== targets gaining Wallace as a cohort: 3-cohort primary vs 4-cohort sensitivity\n"); print(x, row.names=FALSE)
cat("label changes (all targets):", sum(old$label[j] != R$label, na.rm=TRUE), "| replicated with Wallace:", sum(R$label == "replicated"), "\n")
