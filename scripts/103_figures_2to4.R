#!/usr/bin/env Rscript
# 103: Manuscript Figures 2 to 4, same style as 102. Reads existing outputs only; refits nothing.
# Run: Rscript scripts/103_figures_2to4.R   (rumen-r env; writes results/figures_manuscript/fig{2,3,4}_*.{pdf,png})
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "figures_manuscript"); dir.create(od, showWarnings=FALSE)
rd <- function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE, quote=""); rc <- function(f) read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; bg <- "#fcfcfb"; blue <- "#104281"; blue2 <- "#6da7ec"; grey <- "#b4b2a9"; orange <- "#c76a1d"
theme_f <- theme_minimal(base_size=8.5) + theme(
  text = element_text(colour=ink), axis.text = element_text(colour=ink2), axis.title = element_text(colour=ink2, size=8),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour=grid, linewidth=0.3),
  legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(colour=ink2, size=7.5),
  plot.title = element_text(face="bold", size=8.5), plot.subtitle = element_text(colour=ink2, size=8), plot.caption = element_text(colour=muted, size=6.5, hjust=0),
  plot.title.position = "plot", plot.caption.position = "plot", plot.background = element_rect(fill=bg, colour=NA))
save2 <- function(p, name, w, h) { ggsave(file.path(od, paste0(name, ".pdf")), p, width=w, height=h, device=cairo_pdf); ggsave(file.path(od, paste0(name, ".png")), p, width=w, height=h, dpi=300, bg=bg); cat("wrote", name, "\n") }
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
repel <- requireNamespace("ggrepel", quietly=TRUE)
lab_layer <- function(data, ...) if (repel) ggrepel::geom_text_repel(data=data, ..., size=2.3, colour=ink2, min.segment.length=0.2, seed=1, max.overlaps=30) else geom_text(data=data, ..., size=2.3, colour=ink2, vjust=-0.7, check_overlap=TRUE)

# ---------------- Figure 2
tabs <- data.frame(table=c("nelore_genus", "charolais_mrna_genus", "alberta_genus", "difford_bact", "wallace_prokaryotes"),
  label=c("Nelore\nshotgun, RFI", "Charolais\nmRNA, RFI", "Alberta\n16S, RFI", "Difford\n16S, CH4", "Wallace\n16S, CH4"),
  dir=c("scan_null", "scan_null", "scan_null", "methane_dryrun", "methane_dryrun"), s0=c(2001, 2001, 2001, 1001, 1001), stringsAsFactors=FALSE)
N <- do.call(rbind, lapply(seq_len(nrow(tabs)), function(i) { t <- tabs[i, ]
  fs <- file.path(root, "results", t$dir, sprintf("%s_seed%d_ancombc2.csv", t$table, t$s0:(t$s0 + 199))); fs <- fs[file.exists(fs)]; stopifnot(length(fs) == 200)
  do.call(rbind, lapply(fs, function(f) { r <- rc(f); nm <- if ("taxon_name" %in% names(r)) r$taxon_name else r$taxon
    pc <- grep("^p_(RFI|CH4)", names(r), value=TRUE); stopifnot(length(pc) == 1); lc <- sub("^p_", "lfc_", pc); sc <- sub("^p_", "se_", pc); keep <- !grepl("^TGT_", nm) & is.finite(r[[lc]] / r[[sc]]); p <- r[[pc]][keep]   # same genus filter as 91
    data.frame(table=t$table, fpr=mean(p < 0.05), nq=sum(p.adjust(p, "BH") < 0.05)) })) }))
N$lab <- factor(tabs$label[match(N$table, tabs$table)], tabs$label)
sN <- N %>% group_by(lab) %>% summarise(med=median(fpr), mean=mean(fpr), top=max(fpr), maxq=max(nq), n=n(), .groups="drop")
pA <- ggplot(N, aes(lab, fpr)) +
  geom_hline(yintercept=0.05, colour=ink2, linetype="22", linewidth=0.35) +
  geom_jitter(width=0.18, height=0, size=0.45, alpha=0.35, colour=blue2) +
  geom_boxplot(width=0.45, outlier.shape=NA, fill=NA, colour=blue, linewidth=0.45) +
  geom_text(data=sN, aes(lab, top + 0.035, label=paste0("max ", format(maxq, big.mark=","), "\nat q < 0.05")), size=2.25, colour=ink2, lineheight=0.9) +
  annotate("text", x=5.45, y=0.05, label="expected 5%", hjust=1, vjust=-0.6, size=2.25, colour=ink2) +
  scale_y_continuous(labels=function(v) paste0(round(100 * v), "%"), expand=expansion(mult=c(0.02, 0.18))) +
  labs(title="a  ANCOM-BC2 false positives in 200 shuffles of the phenotype", x=NULL, y="Genera at nominal p < 0.05") +
  theme_f + theme(panel.grid.major.x=element_blank())
nsd <- function(dir, fmt, tab, method) { fs <- file.path(root, "results", dir, sprintf(fmt, tab, 1001:1200)); fs <- fs[file.exists(fs)]; stopifnot(length(fs) == 200)
  x <- do.call(rbind, lapply(fs, function(f) { d <- rd(f); d[d$status == "ok", c("id", "W")] })); a <- aggregate(W ~ id, x, sd); names(a)[2] <- "sd"; a$table <- tab; a$method <- method; a }
B <- rbind(nsd("methane_dryrun", "%s_seed%d_targets.tsv", "difford_bact", "ANCOM-BC2"), nsd("methane_dryrun", "%s_seed%d_targets.tsv", "wallace_prokaryotes", "ANCOM-BC2"),
           nsd("clr", "%s_seed%d_clr_targets.tsv", "difford_bact", "CLR linear model"), nsd("clr", "%s_seed%d_clr_targets.tsv", "wallace_prokaryotes", "CLR linear model"))
B <- B[B$id != "M1", ]; B$key <- paste(B$table, B$id); B <- B[B$key %in% B$key[B$method == "CLR linear model"] & B$key %in% B$key[B$method == "ANCOM-BC2"], ]
B$lab <- paste0(tk$target[match(B$id, tk$id)], ", ", ifelse(B$table == "difford_bact", "Difford", "Wallace"))
B$lab <- factor(B$lab, B %>% filter(method == "ANCOM-BC2") %>% arrange(sd) %>% pull(lab))
seg <- B %>% group_by(lab) %>% summarise(lo=min(sd), hi=max(sd), .groups="drop")
pB <- ggplot(B, aes(sd, lab)) +
  geom_vline(xintercept=1, colour=ink2, linetype="22", linewidth=0.35) +
  geom_segment(data=seg, aes(x=lo, xend=hi, y=lab, yend=lab), colour=grid, linewidth=1.3, inherit.aes=FALSE) +
  geom_point(aes(colour=method), size=1.9) +
  scale_colour_manual(values=c("ANCOM-BC2"=blue, "CLR linear model"=orange)) +
  labs(title="b  Spread of each methane target's null statistic", x="SD of W across 200 shuffles (1 = calibrated)", y=NULL) +
  theme_f + theme(panel.grid.major.y=element_blank())
f2 <- (pA | pB) + plot_layout(widths=c(1.15, 1)) + plot_annotation(
  title="ANCOM-BC2 is miscalibrated on these tables; a CLR linear model is not",
  subtitle="The phenotype was permuted 200 times within diet, farm, herd or management group and the full model refitted each time.",
  caption=paste0("a: each point is one shuffle; boxes show the median and quartiles. b: per target, the standard deviation of W = estimate / SE under permutation.\n",
                 "Every estimate in this study is divided by its own null SD when above 1. Sources: results/scan_null/, results/methane_dryrun/, results/clr/."), theme=theme_f)
save2(f2, "fig2_calibration", 7.2, 4.4)

# ---------------- Figure 3
w <- rd(file.path(root, "results", "denovo_scan", "wallace_prokaryotes_scan.tsv")); d <- rd(file.path(root, "results", "denovo_scan", "difford_bact_scan.tsv"))
g <- intersect(w$genus, d$genus); S <- data.frame(genus=g, x=w$W_cal[match(g, w$genus)], y=d$W_cal[match(g, d$genus)], xr=w$W[match(g, w$genus)], yr=d$W[match(g, d$genus)])
S <- S[is.finite(S$x) & is.finite(S$y), ]
fc <- rd(file.path(root, "results", "wallace_farms", "farm_consistency.tsv")); rob <- fc$genus[fc$farm_robust & fc$genus != "TGT_M2"]
S$robust <- S$genus %in% rob
e <- rd(file.path(root, "results", "methane_pooled", "estimates_calibrated.tsv"))
m2 <- data.frame(genus="Succinivibrionaceae (M2)", x=e$W[e$id == "M2" & e$cohort == "wallace"], y=e$W[e$id == "M2" & e$cohort == "difford"]); stopifnot(nrow(m2) == 1)
ct <- cor.test(S$x, S$y, method="spearman", exact=FALSE); ctr <- cor.test(S$xr, S$yr, method="spearman", exact=FALSE)
p3 <- ggplot(S, aes(x, y)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-1.96, ymax=1.96, fill="#f0efec") +
  geom_hline(yintercept=0, colour=ink2, linewidth=0.3) + geom_vline(xintercept=0, colour=ink2, linewidth=0.3) +
  geom_point(data=S[!S$robust, ], colour=grey, size=1.4) +
  geom_point(data=S[S$robust, ], colour=blue, size=1.9) +
  geom_point(data=m2, shape=23, fill=orange, colour=ink, size=2.8) +
  lab_layer(data=S[S$robust, ], aes(label=genus)) + lab_layer(data=m2, aes(label=genus)) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=1.5, size=2.5, colour=ink2,
           label=sprintf("%d shared genera\nSpearman rho = %.2f (p = %.3f)", nrow(S), ct$estimate, ct$p.value)) +
  annotate("text", x=Inf, y=1.96, hjust=1.02, vjust=-0.5, size=2.25, colour=ink2, label="Difford |W| < 1.96") +
  labs(title="Genera associated with methane in Wallace show no association in Difford",
       subtitle="Calibrated W per genus in the two dairy cohorts (methane in g per day). Blue: the 16 genera that kept their sign and\nstrength when each Wallace farm was left out. Diamond: Succinivibrionaceae, reported by Wallace et al. (re-analysis there).",
       caption="Sources: results/denovo_scan/*_scan.tsv (per-genus permutation calibration), results/wallace_farms/farm_consistency.tsv, results/methane_pooled/estimates_calibrated.tsv.",
       x="Calibrated W, Wallace (1,009 cows, 7 farms)", y="Calibrated W, Difford (749 cows, 6 herds)") + theme_f
save2(p3, "fig3_wallace_vs_difford", 6.4, 5.6)

# ---------------- Figure 4
coh <- function(dirn, lab, cohort, trait) { x <- rd(file.path(root, "results", dirn, "estimates_calibrated.tsv")); x <- x[x$id == "M1" & x$cohort == cohort, ]; stopifnot(nrow(x) == 1)
  data.frame(row=lab, r=x$r, lo=x$lo, hi=x$hi, n=x$n, p=x$p_perm, kind="Single cohort", trait=trait) }
pool <- function(dirn, lab, trait) { x <- rd(file.path(root, "results", dirn, "classification_primary.tsv")); x <- x[x$id == "M1", ]; stopifnot(nrow(x) == 1)
  data.frame(row=lab, r=x$r, lo=x$ci_lo, hi=x$ci_hi, n=NA, p=x$p, kind="Pooled", trait=trait) }
M <- rbind(coh("methane_pooled", "Difford, g per day", "difford", "g per day"),
           coh("methane_pooled", "Wallace, g per day", "wallace", "g per day"),
           coh("methane_pooled_noSE1", "Wallace without SE1, g per day", "wallace", "g per day"),
           coh("methane_pooled_depth", "Wallace, depth-adjusted, g per day", "wallace", "g per day"),
           coh("methane_pooled_ch4_dmi", "Wallace, per kg intake", "wallace", "per kg intake"),
           pool("methane_pooled", "Pooled, g per day (primary)", "g per day"),
           pool("methane_pooled_ch4_dmi", "Pooled, Difford g per day + Wallace per kg intake", "per kg intake"))
M$row <- factor(M$row, rev(M$row)); M$plab <- ifelse(M$kind == "Pooled", "", ifelse(M$p < 0.005, "p < 0.005", sprintf("p = %.2f", M$p)))
p4 <- ggplot(M, aes(r, row)) +
  geom_vline(xintercept=0, colour=ink2, linetype="22", linewidth=0.35) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, colour=trait), height=0, linewidth=1.1) +
  geom_point(aes(colour=trait, shape=kind), size=2.4, fill=bg, stroke=0.8) +
  geom_text(aes(label=plab), x=max(M$hi) + 0.03, hjust=0, size=2.35, colour=ink2) +
  scale_colour_manual(values=c("g per day"=blue, "per kg intake"=orange)) + scale_shape_manual(values=c("Single cohort"=16, "Pooled"=23)) +
  scale_x_continuous(expand=expansion(mult=c(0.05, 0.25))) +
  labs(title="The archaeal clade ratio tracks methane per kg intake, not daily output, in one cohort",
       subtitle="Partial r and 95% CI for log(M. gottschalkii clade / M. ruminantium clade) against methane, calibrated; published direction = positive.\nPermutation p shown for single cohorts.",
       caption="Per kg intake is a pre-registered Wallace-only sensitivity trait (intake partly marker-estimated). SE1 is the Swedish farm that may overlap with Danielsson et al. (2017).\nSources: results/methane_pooled*/estimates_calibrated.tsv and classification_primary.tsv.",
       x="Partial r with methane", y=NULL) + theme_f + theme(panel.grid.major.y=element_blank())
save2(p4, "fig4_clade_ratio_trait", 7.2, 3.6)

# ---------------- numbers for the manuscript text
cat("\n== Figure 2a: per-table false positive rate across shuffles\n"); print(as.data.frame(sN), digits=3, row.names=FALSE)
cat("== Figure 2b: null SD range, ANCOM-BC2", paste(round(range(B$sd[B$method == "ANCOM-BC2"]), 2), collapse=" to "), "| CLR", paste(round(range(B$sd[B$method == "CLR linear model"]), 2), collapse=" to "), "| targets", nlevels(B$lab), "\n")
cat("== Figure 3:", nrow(S), "genera; Spearman on calibrated W", round(ct$estimate, 3), "p", signif(ct$p.value, 2), "| on raw W", round(ctr$estimate, 3), "p", signif(ctr$p.value, 2), "| robust plotted", sum(S$robust), "| M2 W Wallace", round(m2$x, 2), "Difford", round(m2$y, 2), "\n")
cat("== Figure 4\n"); print(M[, c("row", "r", "lo", "hi", "n", "p")], digits=3, row.names=FALSE)
