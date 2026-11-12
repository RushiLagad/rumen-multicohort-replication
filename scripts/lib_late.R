# Amendment 44: late analyses, run after all primary results were seen. LinDA, ALDEx2, Hedges' g (Charolais), Li et al. diversity claims.
# Called from 21, 65 and 81 with LATE=linda|aldex2|hedges|diversity on the objects those scripts build before ANCOM-BC2 (as lib_clr.R).
# LATE_SEEDS=real (default) fits the real trait; LATE_SEEDS=1001:1025 fits those within-design permutations. Writes only to results/late_*/;
# an output file that already exists is never overwritten. No phenotype value is printed.
LATE_LABEL <- "late analysis, run after all primary results were seen"

late_main <- function(A, feat, tk, md, terms, n, df, cohort, rank, prefix, root, perm_by, X=NULL, tax=NULL, lin=NULL, m1=FALSE, difford=FALSE) {
  what <- Sys.getenv("LATE"); stopifnot(what %in% c("linda", "aldex2", "hedges", "diversity", "spikein", "mapping"), !nzchar(Sys.getenv("SHUFFLE_SEED")))
  pl <- readLines(file.path(root, "plan", "ANALYSIS_PLAN.md")); stopifnot("commit Amendment 44 first" = any(grepl("^\\*\\*Amendment 44 \\(", pl)))
  if (what %in% c("spikein", "mapping")) stopifnot("commit Amendment 45 first" = any(grepl("^\\*\\*Amendment 45 \\(", pl)))
  if (startsWith(cohort, "alberta") && what != "spikein") stopifnot("Alberta late analyses use all animals (ALBERTA_ALL=1, Amendment 39)" = cohort == "alberta_all")
  A <- as.matrix(A); stopifnot(all(colnames(A) %in% rownames(md))); md <- md[colnames(A), , drop=FALSE]
  od <- file.path(root, "results", paste0("late_", if (what == "hedges") "hedges_g" else what)); dir.create(od, showWarnings=FALSE)
  if (what == "spikein") return(late_spike(if (is.null(X)) A[!grepl("^TGT_", rownames(A)), , drop=FALSE] else X, md, terms, df, prefix, od, perm_by))
  if (what == "mapping") return(late_mapping(A, feat, tk, cohort, rank, prefix, od))
  cat(LATE_LABEL, "|", what, "|", prefix, "| animals", ncol(A), "| features", nrow(A), "\n")
  if (what %in% c("linda", "aldex2")) late_da(what, A, feat, tk, md, terms, cohort, rank, prefix, od, perm_by, m1, difford)
  else if (what == "hedges") late_hedges(A, feat, tk, md, terms, cohort, rank, prefix, od, root)
  else late_div(X, md, terms, cohort, od, tax, lin)
}

wr <- function(x, f) { tmp <- paste0(f, ".tmp", Sys.getpid()); write.table(x, tmp, sep="\t", quote=FALSE, row.names=FALSE); file.rename(tmp, f) }
perm_md <- function(md, v, g, s) { set.seed(s); g <- as.character(g)
  for (lv in unique(g)) { i <- which(g == lv); if (length(i) > 1) md[[v]][i] <- md[[v]][i[sample.int(length(i))]] }; md }

# ---------------------------------------------------------------- LinDA and ALDEx2
fit_linda <- function(B, md, terms, vcol) {
  L <- MicrobiomeStat::linda(feature.dat=B, meta.dat=md[, terms, drop=FALSE], formula=paste("~", paste(terms, collapse=" + ")), feature.dat.type="count", verbose=FALSE)
  o <- L$output[[vcol]]; stopifnot("trait term missing from linda output" = !is.null(o))
  data.frame(lfc=o$log2FoldChange, se=o$lfcSE, W=o$stat, p=o$pvalue, df=o$df, row.names=rownames(o)) }
# Same per-instance OLS as aldex.glm (one gaussian glm per feature and Monte Carlo instance), averaged over instances as aldex.glm does.
fast_glm <- function(cl, mm, j) {
  mc <- ALDEx2::getMonteCarloInstances(cl); k <- ALDEx2::numMCInstances(cl); stopifnot(identical(names(mc), rownames(mm)))
  XtXi <- solve(crossprod(mm)); H <- XtXi %*% t(mm); dfr <- nrow(mm) - ncol(mm); sE <- sS <- sT <- 0
  for (i in seq_len(k)) { Y <- t(vapply(mc, function(x) x[, i], numeric(nrow(mc[[1]])))); Bh <- H %*% Y; R <- Y - mm %*% Bh
    se <- sqrt(colSums(R^2) / dfr * XtXi[j, j]); b <- Bh[j, ]; sE <- sE + b; sS <- sS + se; sT <- sT + b / se }
  data.frame(Est=sE / k, SE=sS / k, t=sT / k, row.names=rownames(mc[[1]])) }
fit_aldex <- function(B, mm, j, vcol, real) {
  cl <- suppressMessages(ALDEx2::aldex.clr(as.data.frame(B), mm, mc.samples=128, denom="all", verbose=FALSE)); f <- fast_glm(cl, mm, j); dfr <- nrow(mm) - ncol(mm)
  if (!real) return(data.frame(lfc=f$Est, se=f$SE, W=f$t, p=2 * pt(-abs(f$t), dfr), df=dfr, row.names=rownames(f)))
  g <- suppressMessages(ALDEx2::aldex.glm(cl, verbose=FALSE)); cn <- colnames(g); has <- grepl(vcol, cn, fixed=TRUE)
  pick <- function(pat, not=NULL) { x <- cn[has & grepl(pat, cn) & (if (is.null(not)) TRUE else !grepl(not, cn))]; stopifnot("aldex.glm column not found" = length(x) == 1); g[[x]] }
  tg <- pick("t\\.val|t value"); i <- match(rownames(f), rownames(g)); stopifnot(!anyNA(i)); d <- max(abs(tg[i] - f$t))
  cat("ALDEx2 check, aldex.glm t vs per-instance OLS t: max |difference|", signif(d, 3), "over", nrow(f), "features\n")
  stopifnot("aldex.glm and the permutation algebra disagree (Amendment 44 item 4)" = d < 1e-6)
  data.frame(lfc=pick("Est")[i], se=pick(":SE|Std\\. Error")[i], W=tg[i], p=pick("pval|Pr\\(", "padj")[i], df=dfr, row.names=rownames(f)) }
late_da <- function(what, A, feat, tk, md, terms, cohort, rank, prefix, od, perm_by, m1, difford) {
  B <- A[rowMeans(A > 0) >= 0.10, , drop=FALSE]; v <- terms[length(terms)]; g <- interaction(md[, perm_by, drop=FALSE], drop=TRUE)
  mm0 <- model.matrix(reformulate(terms), md); j <- grep(paste0("^", v), colnames(mm0)); stopifnot(length(j) == 1, qr(mm0)$rank == ncol(mm0)); vcol <- colnames(mm0)[j]
  dfr <- nrow(mm0) - ncol(mm0); i <- match(feat, rownames(B)); u <- Sys.getenv("LATE_SEEDS", "real")
  if (what == "linda") suppressMessages(library(MicrobiomeStat)) else suppressMessages(library(ALDEx2))
  fit <- function(m, real) if (what == "linda") fit_linda(B, m, terms, vcol) else fit_aldex(B, model.matrix(reformulate(terms), m), j, vcol, real)
  tlog <- function(unit, sec) cat(paste(prefix, unit, round(sec, 1), nrow(B), ncol(B), sep="\t"), "\n", sep="", file=file.path(od, "timing.tsv"), append=TRUE)
  if (u == "real") {
    f <- file.path(od, paste0(prefix, "_targets.tsv")); if (file.exists(f) && file.size(f) > 0) { cat("exists, not overwritten:", f, "\n"); return(invisible()) }
    set.seed(20261111); t0 <- proc.time()[[3]]; r <- fit(md, TRUE); sec <- proc.time()[[3]] - t0
    if (any(r$df != dfr)) cat("note:", what, "df differs from n - p for", sum(r$df != dfr), "features\n")
    k <- match(feat, rownames(r))
    est <- data.frame(id=tk$id, set=tk$set, site=tk$site, cohort=cohort, rank=rank, feature=feat, n=ncol(B), df=ifelse(is.na(k), dfr, r$df[k]),
      lfc=r$lfc[k], se=r$se[k], p_ancom=r$p[k], W=r$W[k], method=what, analysis=LATE_LABEL)
    est$status <- ifelse(is.na(feat), "no_feature", ifelse(is.na(i), "filtered", ifelse(is.finite(est$W), "ok", "no_estimate")))
    if (m1) est$status[est$id == "M1"] <- "lm_primary"; if (difford) est$status[est$id == "M3"] <- "untestable_here"
    full <- data.frame(genus=rownames(r), r[, c("lfc", "se", "W", "p", "df")], analysis=LATE_LABEL)
    wr(full, file.path(od, paste0(prefix, "_full.tsv"))); wr(est, f); tlog("real", sec)
    cat(what, prefix, "| trait term", vcol, "| df", dfr, "| targets ok", sum(est$status == "ok"), "| seconds", round(sec), "\n")
    print(est[est$status == "ok", c("id", "feature", "lfc", "se", "W", "p_ancom")], row.names=FALSE, digits=3)
  } else {
    se <- as.integer(strsplit(u, ":", fixed=TRUE)[[1]]); seeds <- seq(se[1], se[length(se)]); nd <- file.path(od, "null", prefix); dir.create(nd, recursive=TRUE, showWarnings=FALSE)
    ok <- !is.na(i); t0 <- proc.time()[[3]]; done <- 0
    for (s in seeds) { f <- file.path(nd, paste0("seed", s, ".tsv")); if (file.exists(f) && file.size(f) > 0) next
      mp <- perm_md(md, v, g, s); r <- fit(mp, FALSE); k <- match(feat[ok], rownames(r))
      wr(data.frame(id=tk$id[ok], W=r$W[k], seed=s), f); done <- done + 1 }
    sec <- proc.time()[[3]] - t0; if (done) tlog(paste0("perm_x", done), sec)
    cat(what, prefix, "| permutations within", paste(perm_by, collapse=" x "), "| seeds", u, "| new fits", done, "| seconds", round(sec), "\n") }
}

# ---------------------------------------------------------------- Hedges' g, Charolais (descriptive, not pooled)
late_hedges <- function(A, feat, tk, md, terms, cohort, rank, prefix, od, root) {
  stopifnot(cohort == "charolais_mrna"); f <- file.path(od, paste0(prefix, "_hedges_g.tsv")); if (file.exists(f)) { cat("exists, not overwritten:", f, "\n"); return(invisible()) }
  B <- A[rowMeans(A > 0) >= 0.10, , drop=FALSE]; L <- log(B + 0.5); C <- sweep(L, 2, colMeans(L[!grepl("^TGT_", rownames(L)), , drop=FALSE]))   # CLR as Amendment 42
  hi <- md$RFI_class == "high_RFI"; n1 <- sum(hi); n2 <- sum(!hi); J <- 1 - 3 / (4 * (n1 + n2) - 9)
  mm <- model.matrix(reformulate(terms), md); j <- grep("^RFI_class", colnames(mm)); stopifnot(length(j) == 1)
  XtXi <- solve(crossprod(mm)); Bh <- XtXi %*% crossprod(mm, t(C)); R <- t(C) - mm %*% Bh; dfr <- nrow(mm) - ncol(mm)
  tt <- Bh[j, ] / sqrt(colSums(R^2) / dfr * XtXi[j, j]); i <- match(feat, rownames(C)); ok <- !is.na(i)
  m1 <- rowMeans(C[, hi, drop=FALSE]); m2 <- rowMeans(C[, !hi, drop=FALSE])
  sp <- sqrt(((n1 - 1) * apply(C[, hi, drop=FALSE], 1, var) + (n2 - 1) * apply(C[, !hi, drop=FALSE], 1, var)) / (n1 + n2 - 2))
  out <- data.frame(id=tk$id, set=tk$set, target=tk$target, rank=rank, feature=feat, n_high=n1, n_low=n2, J=J,
    g_unadjusted=unname(J * ((m1 - m2) / sp)[i]), t_adjusted=unname(tt[i]), g_adjusted=unname(J * tt[i] * sqrt(1 / n1 + 1 / n2)), df=dfr)
  out$r_clr <- out$t_adjusted / sqrt(out$t_adjusted^2 + dfr)
  ce <- read.delim(file.path(root, "results", "replication_rfi3_calibrated", "estimates_calibrated_all.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
  ce <- ce[ce$cohort == "charolais_mrna", ]; k <- match(out$id, ce$id); out$r_ancombc2_calibrated <- ce$r[k]; out$W_ancombc2_calibrated <- ce$W[k]
  out$status <- ifelse(is.na(feat), "no_feature", ifelse(!ok, "filtered", "ok")); out$analysis <- LATE_LABEL
  out$note <- "descriptive; extreme-group selection (RFI tails only) inflates g relative to an underlying correlation; positive = higher in high RFI"
  wr(out, f); x <- out[out$status == "ok", ]
  cat("Hedges' g,", prefix, "| high", n1, "low", n2, "| J", round(J, 4), "| targets", nrow(x), "| g unadjusted", round(min(x$g_unadjusted), 2), "to", round(max(x$g_unadjusted), 2),
      "| g adjusted", round(min(x$g_adjusted), 2), "to", round(max(x$g_adjusted), 2), "\n")
  print(x[, c("id", "target", "g_unadjusted", "g_adjusted", "r_clr", "r_ancombc2_calibrated")], row.names=FALSE, digits=3)
}

# ---------------------------------------------------------------- Li et al. diversity claims (Section 3.1)
gtdb_domain <- function(g, tax) {
  sp <- function(l) strsplit(l, "\t|\t", fixed=TRUE); col <- function(v, i) vapply(v, `[`, "", i)
  nd <- sp(readLines(file.path(tax, "nodes.dmp"))); id <- col(nd, 1); par <- setNames(col(nd, 2), id); rnk <- setNames(col(nd, 3), id)
  nl <- readLines(file.path(tax, "names.dmp")); nl <- sp(nl[grepl("scientific name", nl, fixed=TRUE)]); nm <- setNames(col(nl, 2), col(nl, 1))
  gid <- names(rnk)[rnk == "genus"]; gn <- unname(nm[gid]); gid <- gid[!duplicated(gn)]; gn <- gn[!duplicated(gn)]; t <- gid[match(g, gn)]
  for (k in 1:80) { go <- !is.na(t) & t != "1" & !(unname(rnk[t]) %in% c("superkingdom", "domain", "kingdom")) & unname(par[t]) != "1"
    if (!any(go)) break; t[go] <- unname(par[t[go]]) }
  sub("^d__", "", unname(nm[t])) }
rarefy_once <- function(M, d) { out <- matrix(0, nrow(M), ncol(M)); left <- rep(d, ncol(M)); tot <- colSums(M)
  for (i in seq_len(nrow(M) - 1)) { x <- M[i, ]; k <- rhyper(ncol(M), x, tot - x, left); out[i, ] <- k; left <- left - k; tot <- tot - x }
  out[nrow(M), ] <- left; out }
late_div <- function(X, md, terms, cohort, od, tax, lin) {
  f <- file.path(od, paste0(cohort, "_diversity.tsv")); if (file.exists(f)) { cat("exists, not overwritten:", f, "\n"); return(invisible()) }
  X <- round(as.matrix(X)); stopifnot(all(rownames(md) %in% colnames(X))); X <- X[, rownames(md), drop=FALSE]
  dom <- if (!is.null(lin)) lin$kingdom[match(rownames(X), lin$genus)] else gtdb_domain(rownames(X), tax)
  cat("domains of", nrow(X), "genera:\n"); print(table(dom, useNA="ifany"))
  stopifnot("unexpected domain name" = all(dom[!is.na(dom)] %in% c("Bacteria", "Archaea")))
  mm <- model.matrix(reformulate(terms), md); v <- terms[length(terms)]; j <- grep(paste0("^", v), colnames(mm)); stopifnot(length(j) == 1); dfr <- nrow(mm) - ncol(mm)
  set.seed(20261111); rows <- list(); vals <- data.frame(animal=rownames(md))
  for (d in c("Bacteria", "Archaea")) {
    Xd <- X[!is.na(dom) & dom == d, , drop=FALSE]; Xd <- Xd[rowSums(Xd) > 0, , drop=FALSE]; ng <- nrow(Xd); depth <- if (ng) min(colSums(Xd)) else 0
    idx <- if (d == "Bacteria") c("bacterial_shannon_evenness", "bacterial_shannon") else c("archaeal_simpson", "archaeal_inverse_simpson")
    if (ng < 3 || depth < 2) { for (x in idx) rows[[x]] <- data.frame(index=x, domain=d, cohort=cohort, genera_present=ng, depth=depth, n=ncol(X), df=NA, W=NA, p=NA, r=NA,
        status=paste0("untestable: ", if (ng < 3) "fewer than 3 genera present" else "minimum library size below 2")); next }
    acc <- matrix(0, ncol(Xd), 2); cnt <- matrix(0, ncol(Xd), 2)
    for (b in 1:100) { R <- rarefy_once(Xd, depth); P <- sweep(R, 2, depth, "/"); S <- colSums(R > 0)
      H <- -colSums(ifelse(P > 0, P * log(P), 0)); D <- colSums(P^2)
      v1 <- if (d == "Bacteria") ifelse(S > 1, H / log(S), NA) else 1 - D; v2 <- if (d == "Bacteria") H else 1 / D
      for (q in 1:2) { x <- if (q == 1) v1 else v2; ok <- is.finite(x); acc[ok, q] <- acc[ok, q] + x[ok]; cnt[ok, q] <- cnt[ok, q] + 1 } }
    I <- acc / cnt
    for (q in 1:2) { y <- I[, q]; vals[[idx[q]]] <- y; okk <- is.finite(y); fm <- lm.fit(mm[okk, , drop=FALSE], y[okk]); dfq <- sum(okk) - ncol(mm)
      s2 <- sum(fm$residuals^2) / dfq; XtXi <- solve(crossprod(mm[okk, , drop=FALSE])); t <- fm$coefficients[j] / sqrt(s2 * XtXi[j, j])
      rows[[idx[q]]] <- data.frame(index=idx[q], domain=d, cohort=cohort, genera_present=ng, depth=depth, n=sum(okk), df=dfq, W=unname(t), p=2 * pt(-abs(unname(t)), dfq),
        r=unname(t) / sqrt(unname(t)^2 + dfq), status="ok") } }
  out <- do.call(rbind, rows); out$trait_term <- colnames(mm)[j]; out$analysis <- LATE_LABEL
  out$note <- "genus counts, all genera, rarefied per domain to the cohort minimum (100 rarefactions, indices averaged); two-sided"
  wr(vals, file.path(od, paste0(cohort, "_diversity_values.tsv"))); wr(out, f)
  print(out[, c("index", "cohort", "genera_present", "depth", "n", "df", "W", "r", "status")], row.names=FALSE, digits=3)
}


# ---------------------------------------------------------------- Amendment 45 item 2: spike-in (planted effects in real tables)
late_spike <- function(G, md, terms, df, prefix, od, perm_by) {
  suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment)})
  G <- round(as.matrix(G)); G <- G[, rownames(md), drop=FALSE]; elig <- rownames(G)[rowMeans(G > 0) >= 0.10]
  v <- terms[length(terms)]; stopifnot(is.numeric(md[[v]])); g <- interaction(md[, perm_by, drop=FALSE], drop=TRUE)
  se <- as.integer(strsplit(Sys.getenv("LATE_SEEDS", "4001:4001"), ":", fixed=TRUE)[[1]]); seeds <- seq(se[1], se[length(se)])
  bs <- as.numeric(strsplit(Sys.getenv("LATE_B", "0.1,0.25,0.5"), ",", fixed=TRUE)[[1]]); nd <- file.path(od, prefix); dir.create(nd, recursive=TRUE, showWarnings=FALSE)
  k <- max(5, round(0.10 * length(elig)))
  for (s in seeds) {
    mp <- perm_md(md, v, g, s); z <- as.numeric(scale(mp[[v]]))
    set.seed(s + 1e6); sp <- sample(elig, k); dir <- rep(c(1, -1), length.out=k)   # same planted genera for every effect size of one seed (paired)
    for (b in bs) { f <- file.path(nd, sprintf("b%s_seed%d.tsv", format(b), s)); if (file.exists(f) && file.size(f) > 0) next
      G2 <- G; for (i in seq_len(k)) G2[sp[i], ] <- round(G[sp[i], ] * exp(dir[i] * b * z))
      orig <- rownames(G2); rownames(G2) <- sprintf("F%05d", seq_along(orig))
      tse <- TreeSummarizedExperiment(assays=list(counts=G2), colData=mp[, terms, drop=FALSE]); t0 <- proc.time()[[3]]
      fit <- ancombc2(data=tse, assay_name="counts", fix_formula=paste(terms, collapse=" + "), p_adj_method="BH", prv_cut=0.10, group=NULL, struc_zero=FALSE, verbose=FALSE)
      res <- fit$res; tc <- grep(paste0("^lfc_", v), names(res), value=TRUE); stopifnot(length(tc) == 1); sfx <- sub("^lfc_", "", tc)
      gen <- orig[match(res$taxon, rownames(G2))]; planted <- dir[match(gen, sp)] * b; planted[is.na(planted)] <- 0
      out <- data.frame(genus=gen, lfc=res[[tc]], se=res[[paste0("se_", sfx)]], W=res[[tc]] / res[[paste0("se_", sfx)]], df=df, planted_lfc_per_sd=planted, b=b, seed=s, k_planted=k, analysis=LATE_LABEL)
      wr(out, f); cat(prefix, "seed", s, "b", b, "| planted", k, "of", length(elig), "eligible genera | fit", round(proc.time()[[3]] - t0), "s\n") } }
}

# ---------------------------------------------------------------- Amendment 45 item 3: per-animal CLR of every target feature (Charolais mRNA vs 16S, same animals)
late_mapping <- function(A, feat, tk, cohort, rank, prefix, od) {
  f <- file.path(od, paste0(prefix, "_target_clr.tsv")); if (file.exists(f)) { cat("exists, not overwritten:", f, "\n"); return(invisible()) }
  A <- as.matrix(A); B <- A[rowMeans(A > 0) >= 0.10, , drop=FALSE]; L <- log(B + 0.5); C <- sweep(L, 2, colMeans(L[!grepl("^TGT_", rownames(L)), , drop=FALSE]))   # CLR as Amendment 42
  i <- match(feat, rownames(C)); ok <- !is.na(i); prev <- unname(rowMeans(A > 0)[match(feat, rownames(A))])
  out <- do.call(rbind, lapply(which(ok), function(j) data.frame(id=tk$id[j], target=tk$target[j], rank=rank, cohort=cohort, feature=feat[j], prevalence=round(unname(prev[j]), 3), animal=colnames(C), clr=unname(C[i[j], ]))))
  st <- data.frame(id=tk$id, target=tk$target, rank=rank, cohort=cohort, feature=feat, prevalence=round(prev, 3), status=ifelse(is.na(feat), "no_feature", ifelse(ok, "ok", "below_10pct")))
  wr(out, f); wr(st, file.path(od, paste0(prefix, "_target_status.tsv")))
  cat("mapping", prefix, "| targets with a feature", sum(!is.na(feat)), "| above 10% prevalence", sum(ok), "| animals", ncol(C), "\n")
}
