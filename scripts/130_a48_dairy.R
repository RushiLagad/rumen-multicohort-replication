#!/usr/bin/env Rscript
# Amendment 48.3: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Wallace vs Difford on a harmonized genus set.
# Genera (cleaned SILVA names of the A46 dumps; summed target features excluded, as 114) present in >= 10% of animals and with mean relative abundance
# >= 1e-4 in each cohort; CLR within each cohort on that shared set only (pseudo-count 0.5, so both cohorts use the same reference); each cohort's primary
# methane model; partial r with Fisher 95% CI; Freedman-Lane p (9,999 within herd or farm; seeds 4830001 onward: 1 Wallace, 2 Difford, +2 depth, +4 days-in-milk control,
# +6 days-in-milk control with depth). Agreement: Spearman rho with a genus bootstrap (2,000, seed 4839001 + variant - 1), Pearson on Fisher z, sign agreement with
# an exact binomial test, per-genus Q and I2. Usage: 130_a48_dairy.R. Env A48_FL_B (test runs only). Writes results/late_a48/dairy/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "dairy"); outs <- file.path(od, c("genera.tsv", "agreement.tsv", "not_harmonized.tsv", "reading.tsv", "fig_dairy_effects.pdf", "fig_dairy_effects.png")); a48_fresh(outs)
B <- as.integer(Sys.getenv("A48_FL_B", "9999")); if (B != 9999) cat("NOTE: A48_FL_B =", B, "(test run; the plan fixes 9,999)\n")
cat(A48_LABEL, "| 48.3 harmonized Wallace vs Difford | permutations", B, "\n")
load_coh <- function(p) { d <- a46_load(root, p); A <- as.matrix(d$A); ra <- sweep(A, 2, colSums(A), "/")
  # genera are the SILVA genus rows; summed target features (TGT_, aggregates added alongside their genera) are excluded, as the 46.4 scan
  nm <- rownames(A); keep <- !startsWith(nm, "TGT_")
  list(d=d, A=A[keep, , drop=FALSE], nm=nm[keep], prev=rowMeans(A[keep, , drop=FALSE] > 0), mra=rowMeans(ra[keep, , drop=FALSE]), depth=log10(colSums(A))) }
W <- load_coh("wallace_prokaryotes"); D <- load_coh("difford_bact")
stopifnot("duplicated genus names" = !anyDuplicated(W$nm), !anyDuplicated(D$nm))
okW <- W$nm[W$prev >= 0.10 & W$mra >= 1e-4]; okD <- D$nm[D$prev >= 0.10 & D$mra >= 1e-4]; G <- sort(intersect(okW, okD))
cat("shared genera:", length(G), "(Wallace eligible", length(okW), ", Difford eligible", length(okD), ")\n")
fit_coh <- function(X, variant, seed, term="trait") { d <- X$d; md <- d$md; terms <- d$terms; tr <- terms[length(terms)]
  if (grepl("depth", variant)) { md$log_depth <- X$depth[rownames(md)]; terms <- c(terms[-length(terms)], "log_depth", tr) }
  if (term == "dim") terms <- c(setdiff(terms, "dim"), "dim")
  Y <- X$A[match(G, X$nm), rownames(md), drop=FALSE]; L <- log(Y + 0.5); C <- sweep(L, 2, colMeans(L)); rownames(C) <- G
  Dz <- a46_design(md, terms); E <- a46_resid(C, Dz); g <- interaction(md[, d$perm_by, drop=FALSE], drop=TRUE); f <- a46_fl(E, Dz$rx, Dz$qz, g, B, seed=seed)
  ci <- a46_ci(f$r, Dz$df); data.frame(genus=G, r=f$r, lo=ci[, "lo"], hi=ci[, "hi"], df=Dz$df, n=Dz$n, p_t=2 * pt(-abs(r2t(f$r, Dz$df)), Dz$df), p_fl=f$p_fl, q_fl=p.adjust(f$p_fl, "BH"), seed=seed) }
fr <- rd(file.path(root, "results", "wallace_farms", "farm_consistency.tsv")); robust <- fr$genus[fr$farm_robust %in% c(TRUE, "TRUE") & fr$genus != "TGT_M2"]
sc <- rd(file.path(root, "results", "late_a46", "scan", "wallace_prokaryotes_scan.tsv")); scanq <- sc$genus[sc$q_bh < 0.05]
V <- data.frame(variant=c("CH4", "CH4, depth-adjusted", "days in milk (positive control)", "days in milk, depth-adjusted"), term=c("trait", "trait", "dim", "dim"), depth=c(FALSE, TRUE, FALSE, TRUE))
GEN <- list(); AG <- list(); boot_n <- 2000
for (v in seq_len(nrow(V))) { vn <- paste(V$variant[v], if (V$depth[v]) "depth" else "")
  a <- fit_coh(W, vn, 4830001 + 2 * (v - 1), V$term[v]); b <- fit_coh(D, vn, 4830002 + 2 * (v - 1), V$term[v])
  za <- atanh(a$r); zb <- atanh(b$r); va <- 1 / (a$df - 1); vb <- 1 / (b$df - 1); wa <- 1 / va; wb <- 1 / vb; zf <- (wa * za + wb * zb) / (wa + wb)
  Q <- wa * (za - zf)^2 + wb * (zb - zf)^2; I2 <- ifelse(Q > 1, (Q - 1) / Q * 100, 0)
  GEN[[v]] <- data.frame(variant=V$variant[v], genus=G, r_wallace=a$r, lo_wallace=a$lo, hi_wallace=a$hi, p_fl_wallace=a$p_fl, q_fl_wallace=a$q_fl,
    r_difford=b$r, lo_difford=b$lo, hi_difford=b$hi, p_fl_difford=b$p_fl, q_fl_difford=b$q_fl, Q=Q, Q_p=pchisq(Q, 1, lower.tail=FALSE), I2=I2,
    same_sign=sign(a$r) == sign(b$r), farm_robust_wallace=G %in% robust, wallace_scan_q05=G %in% scanq, df_wallace=a$df, df_difford=b$df, seed_wallace=a$seed, seed_difford=b$seed)
  set.seed(4839001 + v - 1); agree <- function(sel, label) { x <- GEN[[v]][sel, ]; n <- nrow(x); if (n < 3) return(data.frame(variant=V$variant[v], subset=label, genera=n, spearman_rho=NA, rho_lo=NA, rho_hi=NA, pearson_z=NA, same_sign=sum(x$same_sign), share_same_sign=NA, p_sign_binom=NA, median_I2=NA, share_I2_gt50=NA, wallace_q05=NA, difford_q05=NA, both_q05_same_sign=NA, bootstrap="", note="fewer than 3 genera"))
    rho <- cor(x$r_wallace, x$r_difford, method="spearman"); bs <- replicate(boot_n, { i <- sample.int(n, n, TRUE); suppressWarnings(cor(x$r_wallace[i], x$r_difford[i], method="spearman")) })
    ci <- quantile(bs, c(0.025, 0.975), na.rm=TRUE); s <- sum(x$same_sign)
    data.frame(variant=V$variant[v], subset=label, genera=n, spearman_rho=rho, rho_lo=ci[[1]], rho_hi=ci[[2]], pearson_z=cor(atanh(x$r_wallace), atanh(x$r_difford)),
      same_sign=s, share_same_sign=s / n, p_sign_binom=binom.test(s, n)$p.value, median_I2=median(x$I2), share_I2_gt50=mean(x$I2 > 50),
      wallace_q05=sum(x$q_fl_wallace < 0.05), difford_q05=sum(x$q_fl_difford < 0.05), both_q05_same_sign=sum(x$q_fl_wallace < 0.05 & x$q_fl_difford < 0.05 & x$same_sign),
      bootstrap=paste(boot_n, "genus resamples, seed", 4839001 + v - 1), note="descriptive: genera are not independent") }
  AG[[v]] <- rbind(agree(rep(TRUE, length(G)), "all shared genera"), agree(G %in% robust, "16 farm-robust Wallace genera (those shared)"), agree(G %in% scanq, "Wallace A46 scan q < 0.05 (those shared)"))
  cat("--", V$variant[v], ": rho", round(AG[[v]]$spearman_rho[1], 3), "(", round(AG[[v]]$rho_lo[1], 3), "to", round(AG[[v]]$rho_hi[1], 3), "), same sign", AG[[v]]$same_sign[1], "of", AG[[v]]$genera[1], "\n") }
GEN <- do.call(rbind, GEN); AG <- do.call(rbind, AG)
NH <- data.frame(item=c("primer set and amplified region", "read length and read use", "breed composition", "production system and diet", "methane measurement method", "sampling (rumen fluid collection)"),
  difford=c("from the source paper (script 76 records the archaeal V4-V6 primer pair used for M1)", "forward reads only in this re-analysis (DADA2, script 73)", "from the source paper", "from the source paper", "from the source paper", "from the source paper"),
  wallace=c("ProkA primers unpublished; 515F/806R assumed for the V4 region (script 76)", "as processed in script 79", "from the source paper", "from the source paper", "from the source paper", "from the source paper"),
  modelled="no", note="tabulated for the manuscript from the source papers; entries marked 'from the source paper' are completed in draft 18, not by this script")
a <- AG[AG$variant == "CH4" & AG$subset == "all shared genera", ]; ad <- AG[AG$variant == "CH4, depth-adjusted" & AG$subset == "all shared genera", ]
rdg <- function(a) if (a$rho_lo <= 0 && a$rho_hi >= 0 && a$p_sign_binom > 0.05) "no detectable agreement" else if (a$rho_hi < 0) "a negative relationship, reported descriptively" else if (a$rho_lo > 0) "a positive relationship, reported descriptively" else "rho interval includes 0 but sign agreement differs from 0.5 (p <= 0.05): reported as it stands"
RD <- data.frame(variant=c("CH4", "CH4, depth-adjusted"), rho=c(a$spearman_rho, ad$spearman_rho), rho_lo=c(a$rho_lo, ad$rho_lo), rho_hi=c(a$rho_hi, ad$rho_hi), p_sign=c(a$p_sign_binom, ad$p_sign_binom),
  reading=c(rdg(a), rdg(ad)), note="the difference is not attributed to biology; Fig. 5 is redrawn from this analysis")
suppressMessages(library(ggplot2)); F <- GEN[GEN$variant == "CH4", ]
g <- ggplot(F, aes(x=r_difford, y=r_wallace)) + geom_hline(yintercept=0, linewidth=0.3) + geom_vline(xintercept=0, linewidth=0.3) +
  geom_errorbar(aes(ymin=lo_wallace, ymax=hi_wallace), colour="grey75", width=0, linewidth=0.3) + geom_errorbarh(aes(xmin=lo_difford, xmax=hi_difford), colour="grey75", height=0, linewidth=0.3) +
  geom_point(aes(colour=farm_robust_wallace), size=1.4) + scale_colour_manual(values=c(`FALSE`="#7f7f7f", `TRUE`="#c55a11"), labels=c("other shared genera", "farm-robust in Wallace"), name=NULL) +
  coord_equal() + labs(x="Difford: CLR partial r with CH4 g/d (95% CI)", y="Wallace: CLR partial r with CH4 g/d (95% CI)",
    caption=sprintf("%d shared genera; common CLR reference. Spearman rho %.2f (bootstrap 95%% CI %.2f to %.2f); same sign %d of %d.", nrow(F), a$spearman_rho, a$rho_lo, a$rho_hi, a$same_sign, a$genera)) +
  theme_minimal(base_size=8) + theme(legend.position="bottom", panel.grid.minor=element_blank(), plot.background=element_rect(fill="white", colour=NA))
ggsave(outs[5], g, width=120, height=130, units="mm", device=cairo_pdf, bg="white"); ggsave(outs[6], g, width=120, height=130, units="mm", dpi=300, bg="white")
a48_write(GEN, outs[1]); a48_write(AG, outs[2]); a48_write(NH, outs[3]); a48_write(RD, outs[4])
cat("== agreement\n"); print(AG[, c("variant", "subset", "genera", "spearman_rho", "rho_lo", "rho_hi", "same_sign", "p_sign_binom", "median_I2", "wallace_q05", "difford_q05")], digits=3, row.names=FALSE)
cat("\n== reading\n"); print(RD, digits=3, row.names=FALSE)
