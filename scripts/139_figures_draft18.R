#!/usr/bin/env Rscript
# 139: draft 18 figures and Table 4 (Amendment 48.11; reporting only, refits nothing). Reads existing outputs and writes a new folder,
# results/figures_draft18/, refusing to overwrite any file; the draft 17 figures in results/figures_manuscript/ are left as they are.
#   fig2_forest_pooled   Fig. 2, as 102 (Amendment 47 CLR modified Knapp-Hartung), re-laid out to fit 170 mm x FIG2_H_MM (default 150) so that
#                        the figure and its legend fit a 170 x 225 mm page (48.11, reconciliation item 13); legend in one row.
#   fig4_clade_ratio     Fig. 4: single-cohort rows as before, the pooled g-per-day row now the final CLR modified Knapp-Hartung interval
#                        (results/late_a47/intervals/pooled_intervals.tsv, version primary); the pooled row combining Difford g per day with Wallace
#                        per kg intake (different outcomes) is removed from the main figure and stays in Additional file 2 (48.11).
#   fig5_dairy_effects   Fig. 5, redrawn from 48.3 (results/late_a48/dairy/genera.tsv): Wallace vs Difford CLR partial r with 95% intervals on the
#                        harmonized genus set; a, CH4; b, days in milk (the positive control on the same pipeline). rho is descriptive.
#   af4_seven_claims     Additional file 4: the 48.1 forest plot (per-cohort CLR and calibrated ANCOM-BC2, mKH pooled), house style.
#   table4_seven_claims.tsv  Table 4 data from 48.1.
# Run: Rscript scripts/139_figures_draft18.R   (rumen-r env). 170 mm wide, PNG at 300 dpi, no in-image title or source line.
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "figures_draft18")
H2 <- as.numeric(Sys.getenv("FIG2_H_MM", "150"))
outs <- file.path(od, c(paste0(rep(c("fig2_forest_pooled", "fig4_clade_ratio", "fig5_dairy_effects", "af4_seven_claims"), each=2), c(".pdf", ".png")), "table4_seven_claims.tsv"))
if (any(file.exists(outs))) stop("output exists, never overwritten: ", paste(outs[file.exists(outs)], collapse=", "))
dir.create(od, showWarnings=FALSE, recursive=TRUE)
rd <- function(...) read.delim(file.path(root, ...), check.names=FALSE, stringsAsFactors=FALSE, quote="")
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; bg <- "#fcfcfb"; blue <- "#104281"; blue2 <- "#6da7ec"; grey <- "#b4b2a9"; orange <- "#c76a1d"
theme_f <- theme_minimal(base_size=8.5) + theme(
  text = element_text(colour=ink), axis.text = element_text(colour=ink2), axis.title = element_text(colour=ink2, size=8),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour=grid, linewidth=0.3),
  legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(colour=ink2, size=7.5),
  plot.title = element_text(face="bold", size=8.5), plot.title.position = "plot", plot.background = element_rect(fill=bg, colour=NA))
W170 <- 170 / 25.4
save2 <- function(p, name, w, h) { ggsave(file.path(od, paste0(name, ".pdf")), p, width=w, height=h, device=cairo_pdf)
  ggsave(file.path(od, paste0(name, ".png")), p, width=w, height=h, dpi=300, bg=bg); cat(sprintf("wrote %s (%.0f x %.0f mm)\n", name, w * 25.4, h * 25.4)) }
mnum <- function(v, f="%.2f") sub("^-", "−", sprintf(f, v))
mlab <- function(b) sub("^-", "−", format(b, trim=TRUE, drop0trailing=TRUE))
q <- function(s) paste0("'", gsub("'", "\\\\'", s), "'")
COH <- c(nelore="Nelore", charolais_mrna="Charolais", alberta_all="Alberta", difford="Difford", wallace="Wallace")
lab_expr <- function(id, target, dag=FALSE, sep="  ") {
  if (dag) return(paste0(lab_expr(id, target, sep=sep), "*", q("†")))
  if (id == "M1") return(paste0(q(paste0("M1", sep)), "*italic(", q("M. gottschalkii"), ")*", q(" : "), "*italic(", q("M. ruminantium"), ")*", q(" clade ratio")))
  if (grepl("^unclassified ", target)) return(paste0(q(paste0(id, sep, "unclassified ")), "*italic(", q(sub("^unclassified ", "", target)), ")"))
  if (target == "total archaea") return(q(paste0(id, sep, "total archaea")))
  paste0(q(paste0(id, sep)), "*italic(", q(target), ")") }

# ================= Fig. 2 (as 102, compact layout)
p <- rd("results", "late_a47", "intervals", "pooled_intervals.tsv"); p <- p[p$version == "primary", ]
cm <- readLines(file.path(root, "results", "late_a47", "intervals", "chosen_method.txt"))
stopifnot(cm[length(cm)] == "mkh", any(grepl("validated: yes", cm, fixed=TRUE)), sum(p$arm == "RFI") == 27, sum(p$arm == "Methane") == 11)
p$lo <- p$lo_mkh; p$hi <- p$hi_mkh; p$excl_r30 <- p$excl_30_mkh; p$excl_r15 <- p$excl_15_mkh
p$excl_any <- ifelse(is.na(p$published_dir), FALSE, ifelse(p$published_dir > 0, p$hi < 0, p$lo > 0))
desc <- p$k >= 2 & p$statement != "validated"; stopifnot(!any(desc & (p$excl_r30 | p$excl_r15 | p$excl_any)))
cat_lv <- c("Excludes any effect in the published direction", "Excludes r ≥ 0.15", "Excludes r ≥ 0.30", "Not resolved: r = 0.30 inside the interval")
vcol <- setNames(c("#104281", "#256abf", "#6da7ec", "#b4b2a9"), cat_lv)
d <- p %>% mutate(trait=ifelse(arm == "RFI", "RFI", "CH4"), directional=!is.na(published_dir), dir=ifelse(directional, published_dir, 1),
  x=r * dir, xl=pmin(lo * dir, hi * dir), xh=pmax(lo * dir, hi * dir),
  verdict=factor(case_when(excl_any ~ cat_lv[1], excl_r15 ~ cat_lv[2], excl_r30 ~ cat_lv[3], TRUE ~ cat_lv[4]), cat_lv),
  lab=id, expr=mapply(lab_expr, id, target, k >= 2 & statement != "validated"),
  pub=ifelse(directional, sprintf("%s %s", ifelse(published_dir > 0, "higher", "lower"), trait), "not stated"))
stopifnot(identical(d$id[d$k >= 2 & d$statement != "validated"], "C1"), sum(d$arm == "RFI" & d$verdict != cat_lv[4]) == 8, sum(d$arm == "Methane" & d$verdict != cat_lv[4]) == 6)
XL <- c(-0.6, 0.6); xt <- 0.63
d <- d %>% mutate(over_hi=xh > XL[2], over_lo=xl < XL[1], xlp=pmax(xl, XL[1]), xhp=pmin(xh, XL[2]))
theme_2 <- theme_f + theme(panel.grid.major.y=element_blank(), plot.title=element_text(size=8.5, face="bold"))
panel2 <- function(a, title, first, last) {
  s <- d %>% filter(arm == a); lv <- s %>% arrange(desc(xh), desc(x)) %>% pull(lab); s$lab <- factor(s$lab, lv); ex <- setNames(s$expr, s$lab)
  ggplot(s, aes(x, lab)) +
    annotate("rect", xmin=0.15, xmax=0.30, ymin=-Inf, ymax=Inf, fill="#f0efec") + annotate("rect", xmin=0.30, xmax=XL[2], ymin=-Inf, ymax=Inf, fill="#e4e2dc") +
    geom_vline(xintercept=0, colour=ink2, linewidth=0.35, linetype="22") +
    geom_errorbarh(aes(xmin=xlp, xmax=xhp, colour=verdict), height=0, linewidth=0.95) +
    geom_segment(data=filter(s, over_hi), aes(x=XL[2] - 0.04, xend=XL[2], y=lab, yend=lab, colour=verdict), linewidth=0.95, arrow=arrow(length=unit(2.8, "pt"), type="closed"), show.legend=FALSE) +
    geom_segment(data=filter(s, over_lo), aes(x=XL[1] + 0.04, xend=XL[1], y=lab, yend=lab, colour=verdict), linewidth=0.95, arrow=arrow(length=unit(2.8, "pt"), type="closed"), show.legend=FALSE) +
    geom_text(data=filter(s, over_hi), aes(x=XL[2] - 0.01, label=mnum(xh)), hjust=1, vjust=-0.45, size=1.95, colour=ink2) +
    geom_text(data=filter(s, over_lo), aes(x=XL[1] + 0.01, label=mnum(xl)), hjust=0, vjust=-0.45, size=1.95, colour=ink2) +
    geom_point(aes(colour=verdict, shape=directional), size=1.6, stroke=0.6, fill=bg) +
    geom_text(aes(x=xt, label=pub), hjust=0, size=2.3, colour=ink2) + geom_text(aes(x=xt + 0.30, label=k), hjust=1, size=2.3, colour=ink2) +
    { if (first) annotate("text", x=c(xt, xt + 0.30), y=Inf, label=c("Published", "k"), hjust=c(0, 1), vjust=-0.6, size=2.3, fontface="bold", colour=ink2) } +
    { if (first) annotate("text", x=c(0.225, 0.45), y=Inf, label=c("0.15–\n0.30", "r ≥ 0.30\nbenchmark"), vjust=-0.25, size=2.2, colour=ink2, lineheight=0.85) } +
    scale_colour_manual(values=vcol, limits=cat_lv, drop=FALSE,
      labels=c("Excludes any effect in\nthe published direction", "Excludes\nr ≥ 0.15", "Excludes\nr ≥ 0.30", "Not resolved: r = 0.30\ninside the interval")) +
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21), guide="none") +
    scale_x_continuous(breaks=seq(-0.6, 0.6, 0.3), labels=function(b) sub("^-", "−", format(b)), expand=expansion(0)) +
    scale_y_discrete(labels=function(b) parse(text=ex[b]), expand=expansion(add=c(0.6, if (first) 1.2 else 0.6))) +
    coord_cartesian(xlim=XL, clip="off") +
    labs(title=title, y=NULL, x=if (last) "CLR partial r (single cohort or pooled), oriented so positive = the published direction" else NULL) +
    theme_2 + theme(plot.margin=margin(if (first) 4 else 6, 88, 3, 4), legend.box.spacing=unit(14, "pt"), axis.text.x=if (last) element_text() else element_blank(),
      axis.text.y=element_text(size=7.4), legend.position=if (first) "top" else "none", legend.key.width=unit(12, "pt"), legend.text=element_text(size=7.2, lineheight=0.9)) +
    guides(colour=guide_legend(nrow=1, byrow=TRUE))
}
f2 <- panel2("RFI", "a  Feed efficiency (RFI): 27 targets, three beef cohorts, 649 animals", TRUE, FALSE) /
      panel2("Methane", "b  Methane (g per day): 11 targets, two dairy cohorts, 1,758 cows", FALSE, TRUE) + plot_layout(heights=c(27, 11)) + plot_annotation(theme=theme_f)
save2(f2, "fig2_forest_pooled", W170, H2 / 25.4)

# ================= Fig. 4 (clade ratio; final pooled interval)
coh <- function(dirn, lab, cohort, trait) { x <- rd("results", dirn, "estimates_calibrated.tsv"); x <- x[x$id == "M1" & x$cohort == cohort, ]; stopifnot(nrow(x) == 1)
  data.frame(row=lab, r=x$r, lo=x$lo, hi=x$hi, p=(round(x$p_perm * x$k_null) + 1) / (x$k_null + 1), kind="Single cohort", trait=trait) }
m1 <- p[p$id == "M1", ]; stopifnot(nrow(m1) == 1, m1$k == 2, m1$statement == "validated")
M <- rbind(coh("methane_pooled", "Difford, g per day", "difford", "g per day"), coh("methane_pooled", "Wallace, g per day", "wallace", "g per day"),
  coh("methane_pooled_noSE1", "Wallace without SE1, g per day", "wallace", "g per day"), coh("methane_pooled_depth", "Wallace, depth-adjusted, g per day", "wallace", "g per day"),
  coh("methane_pooled_ch4_dmi", "Wallace, per kg intake", "wallace", "per kg intake"),
  data.frame(row="Pooled, g per day (modified Knapp-Hartung)", r=m1$r, lo=m1$lo_mkh, hi=m1$hi_mkh, p=NA, kind="Pooled", trait="g per day"))
M$row <- factor(M$row, rev(M$row)); M$plab <- ifelse(is.na(M$p), "", ifelse(M$p < 0.01, sprintf("p = %.3f", M$p), sprintf("p = %.2f", M$p)))
p4 <- ggplot(M, aes(r, row)) + annotate("rect", xmin=0.15, xmax=0.30, ymin=-Inf, ymax=Inf, fill="#f0efec") + annotate("rect", xmin=0.30, xmax=0.4, ymin=-Inf, ymax=Inf, fill="#e4e2dc") +
  geom_vline(xintercept=0, colour=ink2, linetype="22", linewidth=0.35) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, colour=trait), height=0, linewidth=1.0) + geom_point(aes(colour=trait, shape=kind), size=2.2, fill=bg, stroke=0.8) +
  geom_text(aes(label=plab), x=0.415, hjust=0, size=2.3, colour=ink2) +
  scale_colour_manual(values=c("g per day"=blue, "per kg intake"=orange)) + scale_shape_manual(values=c("Single cohort"=16, "Pooled"=23)) +
  scale_x_continuous(limits=c(-0.32, 0.4), breaks=seq(-0.3, 0.3, 0.1), expand=expansion(0), labels=function(b) mnum(b, "%.1f")) +
  coord_cartesian(clip="off") + labs(x="Partial r with methane (positive = the published direction)", y=NULL) +
  theme_f + theme(panel.grid.major.y=element_blank(), plot.margin=margin(4, 40, 3, 4))
save2(p4, "fig4_clade_ratio", W170, 62 / 25.4)

# ================= Fig. 5 (48.3)
G <- rd("results", "late_a48", "dairy", "genera.tsv"); A <- rd("results", "late_a48", "dairy", "agreement.tsv")
G5 <- G[G$variant %in% c("CH4", "days in milk (positive control)"), ]
lim <- range(c(G5$lo_wallace, G5$hi_wallace, G5$lo_difford, G5$hi_difford)); lim <- c(floor(lim[1] * 20) / 20, ceiling(lim[2] * 20) / 20)
gp <- function(v, title, pos_lab) { s <- G[G$variant == v, ]; a <- A[A$variant == v & A$subset == "all shared genera", ]; stopifnot(nrow(a) == 1, nrow(s) == a$genera)
  s$grp <- ifelse(s$farm_robust_wallace, "Farm-robust in Wallace", "Other shared genera"); s <- s[order(s$farm_robust_wallace), ]
  ggplot(s, aes(r_wallace, r_difford)) +
    geom_abline(slope=1, intercept=0, colour=grid, linewidth=0.5) + geom_hline(yintercept=0, colour=ink2, linewidth=0.3) + geom_vline(xintercept=0, colour=ink2, linewidth=0.3) +
    geom_errorbar(aes(ymin=lo_difford, ymax=hi_difford, colour=grp), width=0, linewidth=0.3, alpha=0.7) +
    geom_errorbarh(aes(xmin=lo_wallace, xmax=hi_wallace, colour=grp), height=0, linewidth=0.3, alpha=0.7) +
    geom_point(aes(colour=grp), size=1.3) +
    scale_colour_manual(values=c("Farm-robust in Wallace"=blue, "Other shared genera"=grey)) +
    scale_x_continuous(limits=lim, labels=mlab) + scale_y_continuous(limits=lim, labels=mlab) + theme(aspect.ratio=1) +
    labs(title=title, subtitle=sprintf("%d shared genera; same sign in %d of %d\nSpearman rho %s, descriptive\n(bootstrap 95%% CI %s to %s)", nrow(s), a$same_sign, nrow(s), mnum(a$spearman_rho), mnum(a$rho_lo), mnum(a$rho_hi)),
      x=sprintf("CLR partial r with %s, Wallace (95%% CI)", pos_lab), y=sprintf("CLR partial r with %s, Difford (95%% CI)", pos_lab)) + theme_f + theme(plot.subtitle=element_text(colour=ink2, size=7.3, lineheight=1.05), aspect.ratio=1) }
f5 <- (gp("CH4", "a  Methane (g per day)", "CH4") | gp("days in milk (positive control)", "b  Days in milk (positive control, same pipeline)", "days in milk")) +
  plot_layout(guides="collect") & theme(legend.position="top")
save2(f5, "fig5_dairy_effects", W170, 105 / 25.4)

# ================= Additional file 4 and Table 4 (48.1)
S7 <- rd("results", "late_a48", "seven", "seven_claims.tsv"); P7 <- rd("results", "late_a48", "seven", "seven_pooled.tsv")
ord <- c("T2", "T4", "T5", "T6", "N15", "S1", "M1"); stopifnot(setequal(P7$id, ord), nrow(S7) == 14)
long <- rbind(
  transmute(S7, id, cohort=COH[cohort], est="CLR", r=r_clr * published_dir, lo=lo_clr, hi=hi_clr, dir=published_dir),
  transmute(S7, id, cohort=COH[cohort], est="Calibrated ANCOM-BC2", r=r_ancombc2 * published_dir, lo=lo_ancombc2, hi=hi_ancombc2, dir=published_dir))
long <- long %>% mutate(xl=pmin(lo * dir, hi * dir), xh=pmax(lo * dir, hi * dir))
pool <- P7 %>% transmute(id, cohort="Pooled (modified Knapp-Hartung, CLR)", est="Pooled", r=r_mkh * published_dir, xl=pmin(lo_mkh * published_dir, hi_mkh * published_dir), xh=pmax(lo_mkh * published_dir, hi_mkh * published_dir))
L <- bind_rows(long %>% select(id, cohort, est, r, xl, xh), pool)
pat <- setNames(sprintf("CLR %s; ANCOM-BC2 %s", substr(P7$pattern_clr, 1, 1), substr(P7$pattern_ancombc2, 1, 1)), P7$id)
tgt <- setNames(P7$target, P7$id)
L$claim <- factor(L$id, ord)
L$row <- ifelse(L$est == "Pooled", "Pooled, modified Knapp-Hartung (CLR)", paste0(L$cohort, ", ", ifelse(L$est == "CLR", "CLR", "ANCOM-BC2")))
L$est <- factor(L$est, c("CLR", "Calibrated ANCOM-BC2", "Pooled"))
L <- L %>% mutate(co=match(cohort, c("Nelore", "Charolais", "Alberta", "Difford", "Wallace", "Pooled (modified Knapp-Hartung, CLR)"))) %>% arrange(claim, co, est)
lv <- c(as.vector(t(outer(c("Nelore", "Charolais", "Alberta", "Difford", "Wallace"), c("CLR", "ANCOM-BC2"), paste, sep=", "))), "Pooled, modified Knapp-Hartung (CLR)")
stopifnot(all(L$row %in% lv)); L$row <- factor(L$row, rev(lv))
strip_lab <- function(ids) sapply(ids, function(i) paste0(lab_expr(i, tgt[[i]], sep="  "), "*", q(paste0("      pattern: ", pat[[i]]))))
L$claim_e <- factor(strip_lab(as.character(L$claim)), strip_lab(ord))
stopifnot(all(table(L$claim) == 5))
pa <- ggplot(L, aes(r, row)) +
  annotate("rect", xmin=0.15, xmax=0.30, ymin=-Inf, ymax=Inf, fill="#f0efec") + annotate("rect", xmin=0.30, xmax=1, ymin=-Inf, ymax=Inf, fill="#e4e2dc") +
  geom_vline(xintercept=0, colour=ink2, linetype="22", linewidth=0.35) +
  geom_errorbarh(aes(xmin=xl, xmax=xh, colour=est), height=0, linewidth=0.8) + geom_point(aes(colour=est, shape=est), size=1.7, fill=bg, stroke=0.7) +
  facet_wrap(~claim_e, ncol=1, scales="free_y", labeller=label_parsed) +
  scale_x_continuous(limits=c(-1, 1), breaks=seq(-1, 1, 0.25), labels=function(b) mnum(b, "%.2f"), expand=expansion(0)) +
  scale_colour_manual(values=c("CLR"=blue, "Calibrated ANCOM-BC2"=orange, "Pooled"=ink), labels=c("CLR", "Calibrated ANCOM-BC2", "Pooled, modified Knapp-Hartung")) +
  scale_shape_manual(values=c("CLR"=16, "Calibrated ANCOM-BC2"=17, "Pooled"=23), labels=c("CLR", "Calibrated ANCOM-BC2", "Pooled, modified Knapp-Hartung")) +
  labs(x="Partial r, oriented so positive = the published direction (95% intervals)", y=NULL) +
  theme_f + theme(panel.grid.major.y=element_blank(), strip.text=element_text(hjust=0, size=7.8, colour=ink), panel.spacing=unit(6, "pt"), axis.text.y=element_text(size=7), plot.margin=margin(5.5, 14, 5.5, 5.5))
save2(pa, "af4_seven_claims", W170, 215 / 25.4)

T4 <- P7 %>% transmute(id, target, arm, published_direction=ifelse(published_dir > 0, "higher", "lower"), k, cohorts=gsub(";", ", ", sapply(strsplit(cohorts, ";"), function(v) paste(COH[v], collapse=";"))),
  pooled_r_mkh=round(r_mkh, 3), pooled_ci_mkh=sprintf("%s to %s", mnum(lo_mkh, "%.3f"), mnum(hi_mkh, "%.3f")), tau_hat=round(tau_hat, 3), I2=round(I2, 1),
  fixed_effect_r=round(r_fe, 3), fixed_effect_p=signif(p_fe_clr, 2), label_a37, label_a47, pattern_clr, pattern_ancombc2,
  secondary_Sa_clr, secondary_Sa_ancombc2, secondary_Sb_clr, secondary_Sb_ancombc2, reasons_not_replicated)
pfmt <- function(p) ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p))   # permutation p as (b + 1)/(B + 1), B = 200 (46.5.1); stored p_perm is b/B
per <- S7 %>% group_by(id) %>% summarise(per_cohort_clr=paste(sprintf("%s %s (%s to %s; n %d)", COH[cohort], mnum(r_clr, "%.3f"), mnum(lo_clr, "%.3f"), mnum(hi_clr, "%.3f"), n), collapse="; "),
  per_cohort_ancombc2=paste(sprintf("%s %s (%s to %s; perm p %s)", COH[cohort], mnum(r_ancombc2, "%.3f"), mnum(lo_ancombc2, "%.3f"), mnum(hi_ancombc2, "%.3f"), pfmt((round(p_perm_ancombc2 * 200) + 1) / 201)), collapse="; "),
  mde_80=paste(sprintf("%s %.3f", COH[cohort], mde_80), collapse="; "), .groups="drop")
T4 <- left_join(T4, per, by="id"); T4 <- T4[match(ord, T4$id), ]
T4$analysis <- "late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen (Amendment 48.1)"
write.table(T4, file.path(od, "table4_seven_claims.tsv"), sep="\t", quote=FALSE, row.names=FALSE); cat("wrote table4_seven_claims.tsv\n")
cat("== Fig. 4\n"); print(M[, c("row", "r", "lo", "hi", "p")], digits=3, row.names=FALSE)
