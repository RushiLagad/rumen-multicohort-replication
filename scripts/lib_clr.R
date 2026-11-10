# Amendment 42: CLR + linear model on the exact objects the primary scripts build before ANCOM-BC2. p is stored in p_ancom for column compatibility.
clr_write <- function(A, feat, tk, md, terms, n, df, cohort, rank, prefix, root, m1=FALSE, difford=FALSE) {
  A <- as.matrix(A); if (all(colnames(A) %in% rownames(md))) md <- md[colnames(A), , drop=FALSE] else stopifnot(nrow(md) == ncol(A))
  B <- A[rowMeans(A > 0) >= 0.10, , drop=FALSE]; L <- log(B + 0.5); C <- sweep(L, 2, colMeans(L[!grepl("^TGT_", rownames(L)), , drop=FALSE]))
  mm <- model.matrix(reformulate(terms), md); stopifnot(nrow(mm) == ncol(C)); v <- terms[length(terms)]; j <- grep(paste0("^", v), colnames(mm)); stopifnot(length(j) == 1)
  XtXi <- solve(crossprod(mm)); Bh <- XtXi %*% crossprod(mm, t(C)); R <- t(C) - mm %*% Bh; dfr <- nrow(mm) - ncol(mm)
  se <- sqrt(colSums(R^2) / dfr * XtXi[j, j]); b <- Bh[j, ]; W <- b / se; p <- 2 * pt(-abs(W), dfr)
  full <- data.frame(genus=rownames(C), lfc=b, se=se, W=W, p=p); i <- match(feat, rownames(C))
  est <- data.frame(id=tk$id, set=tk$set, site=tk$site, cohort=cohort, rank=rank, feature=feat, n=n, df=df, lfc=b[i], se=se[i], p_ancom=p[i]); est$W <- est$lfc / est$se
  est$status <- ifelse(is.na(feat), "no_feature", ifelse(is.na(i), "filtered", ifelse(is.finite(est$W), "ok", "no_estimate")))
  if (m1) est$status[est$id == "M1"] <- "lm_primary"; if (difford) est$status[est$id == "M3"] <- "untestable_here"
  od <- file.path(root, "results", "clr"); dir.create(od, showWarnings=FALSE)
  write.table(est, file.path(od, paste0(prefix, "_clr_targets.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(full, file.path(od, paste0(prefix, "_clr_full.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  cat("CLR (Amendment 42):", prefix, "| animals", ncol(C), "| features", nrow(C), "of", nrow(A), "| df", dfr, "| targets ok", sum(est$status == "ok"), "\n") }
