#!/usr/bin/env Rscript
# Amendment 48.1: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Seven-claim evidence table (no refit,
# except S1 against Wallace CH4 per kg DMI from the 48.9 dump, reconciliation item 12). Inputs: results/late_a46/effect_size/cohort_estimates.tsv (CLR),
# the Amendment 37 calibrated ANCOM-BC2 estimates and classifications, results/late_a47/intervals/pooled_intervals.tsv (mKH, fixed effect).
# Patterns, reasons and secondary criteria exactly as 48.1; every count reported with S1 as a matched claim and as a phenotype-transfer test.
# Usage: 128_a48_seven.R (after the A48 dump step). Writes results/late_a48/seven/ only.
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); sdir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
source(file.path(sdir, "lib_a46.R")); source(file.path(sdir, "lib_a47.R")); source(file.path(sdir, "lib_a48.R")); a48_gate(root)
od <- a48_dir(root, "seven"); outs <- file.path(od, c("seven_claims.tsv", "seven_pooled.tsv", "concordance.tsv", "secondary_matched.tsv", "counts_by_classification.tsv",
  "s1_yield.tsv", "fig_seven_claims.pdf", "fig_seven_claims.png")); a48_fresh(outs)
cat(A48_LABEL, "| 48.1 seven-claim evidence\n")
dirs <- a48_dirs(root); SEVEN <- a48_seven
ce <- rd(file.path(root, "results", "late_a46", "effect_size", "cohort_estimates.tsv"))
PI <- rd(file.path(root, "results", "late_a47", "intervals", "pooled_intervals.tsv")); PI <- PI[PI$version == "primary", ]
er <- rd(file.path(root, "results", "replication_rfi3_calibrated", "estimates_calibrated_all.tsv")); em <- rd(file.path(root, "results", "methane_pooled", "estimates_calibrated.tsv"))
ea <- rbind(data.frame(arm="RFI", er[, c("id", "cohort", "n", "df", "W", "p_perm", "r", "lo", "hi")]), data.frame(arm="Methane", em[, c("id", "cohort", "n", "df", "W", "p_perm", "r", "lo", "hi")]))
c37r <- rd(file.path(root, "results", "replication_rfi3_calibrated", "classification_all.tsv")); c37m <- rd(file.path(root, "results", "methane_pooled", "classification_primary.tsv"))
c37 <- rbind(data.frame(arm="RFI", c37r[, c("set", "id", "k", "r", "ci_lo", "ci_hi", "p", "q_bh", "p_fe", "concordant", "loo_all_sig", "label")]),
             data.frame(arm="Methane", c37m[, c("set", "id", "k", "r", "ci_lo", "ci_hi", "p", "q_bh", "p_fe", "concordant", "loo_all_sig", "label")]))
l47 <- rbind(rd(file.path(root, "results", "late_a47", "labels", "labels_rfi.tsv")), rd(file.path(root, "results", "late_a47", "labels", "labels_methane.tsv")))
used <- function(id) strsplit(PI$cohorts[PI$id == id], ";")[[1]]         # independent test cohorts of each target, as pooled (47.1)
mde <- function(df) tanh((qnorm(0.975) + qnorm(0.80)) / sqrt(df - 1))

# ---------------------------------------------------------------- per cohort (CLR beside calibrated ANCOM-BC2)
SC <- do.call(rbind, lapply(SEVEN, function(id) { cc <- used(id); d <- dirs[[id]]
  do.call(rbind, lapply(cc, function(co) { a <- ce[ce$id == id & ce$cohort == co, ]; b <- ea[ea$id == id & ea$cohort == co, ]; stopifnot(nrow(a) == 1, nrow(b) == 1)
    data.frame(id=id, target=PI$target[PI$id == id], arm=a$arm, cohort=co, published_dir=d, n=a$n, df=a$df,
      r_clr=a$r, lo_clr=a$lo, hi_clr=a$hi, p_t_clr=a$p_t, p_fl_clr=a$p_fl, W_ancombc2_cal=b$W, r_ancombc2=b$r, lo_ancombc2=b$lo, hi_ancombc2=b$hi, p_perm_ancombc2=b$p_perm,
      same_sign_clr=sign(a$r) == d, same_sign_ancombc2=sign(b$r) == d, mde_80=mde(a$df)) })) }))
SC$excl0_pub_clr <- ifelse(SC$published_dir > 0, SC$lo_clr > 0, SC$hi_clr < 0); SC$excl0_pub_ancombc2 <- ifelse(SC$published_dir > 0, SC$lo_ancombc2 > 0, SC$hi_ancombc2 < 0)
SC$bound_pub_clr <- a46_bound(SC$lo_clr, SC$hi_clr, SC$published_dir); SC$bound_pub_ancombc2 <- a46_bound(SC$lo_ancombc2, SC$hi_ancombc2, SC$published_dir)
# leave-one-cohort-out: with k = 2 the other cohort's own estimate (k >= 3 does not occur among the seven; checked)
stopifnot("a seven-claim target has k >= 3" = all(table(SC$id) == 2))
SC$loo_other_r_clr <- ave(SC$r_clr, SC$id, FUN=rev); SC$loo_other_lo_clr <- ave(SC$lo_clr, SC$id, FUN=rev); SC$loo_other_hi_clr <- ave(SC$hi_clr, SC$id, FUN=rev)

# ---------------------------------------------------------------- pattern (first rule that applies), on one estimator
pattern <- function(x, est) { r <- x[[paste0("r_", est)]]; ex <- x[[paste0("excl0_pub_", est)]]; b <- x[[paste0("bound_pub_", est)]]; d <- x$published_dir[1]
  if (any(ex) && !all(ex)) "A: support in one cohort only" else if (all(sign(r) == -d)) "B: against the published direction" else
  if (all(b < 0.15)) "C: bounded (every cohort excludes r = 0.15)" else if (length(unique(sign(r))) > 1) "D: discordant (cohort signs differ)" else "E: consistent, inconclusive" }

# ---------------------------------------------------------------- pooled, reasons, secondary criteria
SP <- do.call(rbind, lapply(SEVEN, function(id) { x <- SC[SC$id == id, ]; p <- PI[PI$id == id, ]; a <- c37[c37$id == id, ]; l <- l47[l47$id == id, ]; d <- x$published_dir[1]
  stopifnot(nrow(p) == 1, nrow(a) == 1, nrow(l) == 1)
  z <- atanh(x$r_clr); v <- 1 / (x$df - 1); w <- 1 / v; zf <- sum(w * z) / sum(w); Q <- sum(w * (z - zf)^2)
  reasons <- c(if (a$q_bh >= 0.05) "pooled q >= 0.05", if (sign(a$r) != d) "pooled sign against the published direction", if (isFALSE(a$loo_all_sig) || is.na(a$loo_all_sig)) "leave-one-cohort-out condition fails")
  data.frame(id=id, target=p$target, arm=p$arm, published_dir=d, k=p$k, cohorts=p$cohorts,
    r_mkh=p$r, lo_mkh=p$lo_mkh, hi_mkh=p$hi_mkh, p_mkh=p$p_mkh, r_fe=p$r_fe, lo_fe=p$lo_fe, hi_fe=p$hi_fe, p_fe_clr=p$p_fe, tau_hat=p$tau_hat, tau2=p$tau2,
    I2=if (Q > 0) max(0, (Q - (p$k - 1)) / Q) * 100 else 0, Q=Q, Q_p=pchisq(Q, p$k - 1, lower.tail=FALSE), statement=p$statement,
    r_a37=a$r, lo_a37=a$ci_lo, hi_a37=a$ci_hi, q_a37=a$q_bh, p_fe_a37=a$p_fe, label_a37=a$label, label_a47=l$label_a47,
    reasons_not_replicated=paste(reasons, collapse="; "),
    pattern_clr=pattern(x, "clr"), pattern_ancombc2=pattern(x, "ancombc2"),
    secondary_Sa_clr=any(x$p_fl_clr < 0.05 & sign(x$r_clr) == d) && sign(p$r) == d, secondary_Sa_ancombc2=any(x$p_perm_ancombc2 < 0.05 & sign(x$r_ancombc2) == d) && sign(a$r) == d,
    secondary_Sb_clr=p$p_fe < 0.05 && sign(p$r_fe) == d, secondary_Sb_ancombc2=a$p_fe < 0.05 && sign(a$r) == d,
    mde_fe=mde(sum(x$df - 1) + 1), mkh_t_crit=qt(0.975, p$k - 1)) }))
SP$estimators_agree_pattern <- substr(SP$pattern_clr, 1, 1) == substr(SP$pattern_ancombc2, 1, 1)

# ---------------------------------------------------------------- S1 against Wallace CH4 per kg DMI (CLR), beside the calibrated ANCOM-BC2 value on file
dy <- a46_load(root, "wallace_prokaryotes_ch4_dmi"); Cy <- a46_clr(dy$A); Dy <- a46_design(dy$md, dy$terms); stopifnot(Dy$df == dy$df)
iy <- match(dy$feat[dy$tk$id == "S1"], rownames(Cy)); stopifnot(length(iy) == 1, !is.na(iy))
ry <- a46_r(a46_resid(Cy[iy, , drop=FALSE], Dy), Dy$rx); ciy <- a46_ci(ry, Dy$df)
ey <- rd(file.path(root, "results", "methane_pooled_ch4_dmi", "estimates_calibrated.tsv")); ey <- ey[ey$id == "S1" & ey$cohort == "wallace", ]; stopifnot(nrow(ey) == 1)
SY <- data.frame(id="S1", target="Sharpea", cohort="wallace", trait="CH4 per kg DMI", published_dir=-1, n=Dy$n, df=Dy$df, r_clr=ry, lo_clr=ciy[1, "lo"], hi_clr=ciy[1, "hi"],
  p_t_clr=2 * pt(-abs(r2t(ry, Dy$df)), Dy$df), r_ancombc2_cal=ey$r, lo_ancombc2_cal=ey$lo, hi_ancombc2_cal=ey$hi, p_perm_ancombc2=ey$p_perm,
  note="source trait type (methane yield); host differs (sheep source, cattle test); k = 1 (Difford has no yield trait)")

# ---------------------------------------------------------------- the 21 matched claims (stratum directed rumen): secondary criteria; counts under both S1 classifications
M <- PI[PI$stratum == "directed rumen", ]
SM <- do.call(rbind, lapply(M$id, function(id) { p <- M[M$id == id, ]; d <- p$published_dir; cc <- strsplit(p$cohorts, ";")[[1]]
  a <- ce[ce$id == id & ce$cohort %in% cc, ]; b <- ea[ea$id == id & ea$cohort %in% cc, ]; a37 <- c37[c37$id == id, ]; stopifnot(nrow(a) == length(cc), nrow(b) == length(cc), nrow(a37) == 1)
  data.frame(id=id, target=p$target, arm=p$arm, k=p$k, published_dir=d,
    Sa_clr=any(a$p_fl < 0.05 & sign(a$r) == d) && sign(p$r) == d, Sa_ancombc2=any(b$p_perm < 0.05 & sign(b$r) == d) && sign(a37$r) == d,
    Sb_clr=if (p$k >= 2) p$p_fe < 0.05 && sign(p$r_fe) == d else NA, Sb_ancombc2=if (p$k >= 2) a37$p_fe < 0.05 && sign(a37$r) == d else NA,
    cohorts_same_sign_clr=sum(sign(a$r) == d), cohorts_same_sign_ancombc2=sum(sign(b$r) == d), cohorts=length(cc), label_a37=a37$label) }))
stopifnot("matched claims are not 21" = nrow(SM) == 21)
cls <- function(s1_matched) { x <- if (s1_matched) SM else SM[SM$id != "S1", ]
  data.frame(classification=if (s1_matched) "S1 matched (current)" else "S1 as phenotype-transfer test", matched=nrow(x), k_ge2=sum(x$k >= 2), k1=sum(x$k == 1),
    replicated_a37=sum(x$label_a37 == "replicated"), Sa_clr=sum(x$Sa_clr), Sa_ancombc2=sum(x$Sa_ancombc2), Sb_clr=sum(x$Sb_clr, na.rm=TRUE), Sb_ancombc2=sum(x$Sb_ancombc2, na.rm=TRUE),
    Sa_clr_ids=paste(x$id[x$Sa_clr], collapse=";"), Sa_ancombc2_ids=paste(x$id[x$Sa_ancombc2], collapse=";")) }
CNT <- rbind(cls(TRUE), cls(FALSE))

# ---------------------------------------------------------------- directional concordance: seven claims (both classifications) and every independent directed test (38 targets)
bt <- function(x) { x <- x[!is.na(x)]; c(n=length(x), same=sum(x), share=mean(x), p_binom=binom.test(sum(x), length(x))$p.value) }
allt <- do.call(rbind, lapply(PI$id[!is.na(PI$published_dir)], function(id) { p <- PI[PI$id == id, ]; cc <- strsplit(p$cohorts, ";")[[1]]
  a <- ce[ce$id == id & ce$cohort %in% cc, c("cohort", "r")]; b <- ea[ea$id == id & ea$cohort %in% cc, c("cohort", "r")]; m <- merge(a, b, by="cohort", suffixes=c("_clr", "_ancombc2"))
  data.frame(id=id, cohort=m$cohort, d=p$published_dir, clr=sign(m$r_clr) == p$published_dir, ancombc2=sign(m$r_ancombc2) == p$published_dir) }))
CC <- do.call(rbind, lapply(list(list("seven claims (S1 matched)", SC$id %in% SEVEN), list("six claims (S1 as transfer test)", SC$id %in% setdiff(SEVEN, "S1"))), function(g) {
  x <- SC[g[[2]], ]; rbind(data.frame(set=g[[1]], estimator="CLR", t(bt(x$same_sign_clr))), data.frame(set=g[[1]], estimator="ANCOM-BC2 calibrated", t(bt(x$same_sign_ancombc2)))) }))
CC <- rbind(CC, data.frame(set="every independent test with a published direction (38 targets)", estimator="CLR", t(bt(allt$clr))),
  data.frame(set="every independent test with a published direction (38 targets)", estimator="ANCOM-BC2 calibrated", t(bt(allt$ancombc2))))
ov <- SC$same_sign_clr & SC$same_sign_ancombc2
CC$overlap_note <- c(paste0("seven claims: published sign under both estimators in ", sum(ov), " of ", nrow(SC), " cohort estimates (CLR ", sum(SC$same_sign_clr), ", ANCOM-BC2 ", sum(SC$same_sign_ancombc2), ")"), rep("", nrow(CC) - 1))
CC$note <- "descriptive: estimates are not independent across targets"

# ---------------------------------------------------------------- forest plot (oriented to the published direction)
suppressMessages({ library(ggplot2) })
FP <- rbind(data.frame(id=SC$id, row=paste(SC$cohort, "CLR"), est="CLR", r=SC$published_dir * SC$r_clr, lo=pmin(SC$published_dir * SC$lo_clr, SC$published_dir * SC$hi_clr), hi=pmax(SC$published_dir * SC$lo_clr, SC$published_dir * SC$hi_clr)),
  data.frame(id=SC$id, row=paste(SC$cohort, "ANCOM-BC2"), est="ANCOM-BC2 calibrated", r=SC$published_dir * SC$r_ancombc2, lo=pmin(SC$published_dir * SC$lo_ancombc2, SC$published_dir * SC$hi_ancombc2), hi=pmax(SC$published_dir * SC$lo_ancombc2, SC$published_dir * SC$hi_ancombc2)),
  data.frame(id=SP$id, row="pooled, modified Knapp-Hartung (CLR)", est="pooled mKH", r=SP$published_dir * SP$r_mkh, lo=pmin(SP$published_dir * SP$lo_mkh, SP$published_dir * SP$hi_mkh), hi=pmax(SP$published_dir * SP$lo_mkh, SP$published_dir * SP$hi_mkh)))
FP$ord <- ifelse(FP$est == "pooled mKH", 99, match(sub(" .*", "", FP$row), unique(SC$cohort)) * 2 + (FP$est != "CLR"))
FP$row <- factor(FP$row, levels=rev(unique(FP$row[order(FP$ord)])))
FP$id <- factor(FP$id, levels=SEVEN, labels=sapply(paste0(SEVEN, " ", sub("Methanobrevibacter gottschalkii clade vs ruminantium clade", "M. gottschalkii : M. ruminantium clade ratio", SP$target[match(SEVEN, SP$id)])), function(x) paste(strwrap(x, 26), collapse="\n"))); FP$lo_c <- pmax(FP$lo, -1); FP$hi_c <- pmin(FP$hi, 1)
g <- ggplot(FP, aes(x=r, y=row, colour=est)) + geom_vline(xintercept=0, linewidth=0.3) + geom_vline(xintercept=0.15, linetype="dotted", linewidth=0.3) + geom_vline(xintercept=0.30, linetype="dashed", linewidth=0.3) +
  geom_errorbarh(aes(xmin=lo_c, xmax=hi_c), height=0.2) + geom_point(size=1.6) + facet_wrap(~ id, ncol=1, scales="free_y", strip.position="left") + coord_cartesian(xlim=c(-1, 1)) +
  scale_colour_manual(values=c("CLR"="#1f4e79", "ANCOM-BC2 calibrated"="#9e9e9e", "pooled mKH"="#c55a11"), name=NULL) +
  labs(x="partial r, oriented so positive = the published direction", y=NULL, caption="Dotted line r = 0.15, dashed line r = 0.30. Intervals are 95%; pooled intervals running off the axis are clipped at -1 and 1.") + theme_minimal(base_size=8) +
  theme(strip.placement="outside", strip.text.y.left=element_text(angle=0, hjust=1, face="bold"), legend.position="bottom", panel.grid.minor=element_blank(), plot.background=element_rect(fill="white", colour=NA))
ggsave(outs[7], g, width=170, height=200, units="mm", device=cairo_pdf, bg="white"); ggsave(outs[8], g, width=170, height=200, units="mm", dpi=300, bg="white")
a48_write(SC, outs[1]); a48_write(SP, outs[2]); a48_write(CC, outs[3]); a48_write(SM, outs[4]); a48_write(CNT, outs[5]); a48_write(SY, outs[6])
cat("== 48.1 seven claims\n"); print(SP[, c("id", "target", "r_mkh", "lo_mkh", "hi_mkh", "tau_hat", "I2", "pattern_clr", "pattern_ancombc2", "secondary_Sa_clr", "secondary_Sa_ancombc2", "secondary_Sb_clr", "secondary_Sb_ancombc2", "reasons_not_replicated")], digits=3, row.names=FALSE)
cat("\n== counts under both S1 classifications\n"); print(CNT, row.names=FALSE)
cat("\n== concordance\n"); print(CC, digits=3, row.names=FALSE)
cat("\n== S1 against Wallace CH4 per kg DMI\n"); print(SY, digits=3, row.names=FALSE)
