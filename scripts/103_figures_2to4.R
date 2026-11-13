#!/usr/bin/env Rscript
# 103: Manuscript Figures 3 to 5, same style as 102. Reads existing outputs only; refits nothing.
# Run: Rscript scripts/103_figures_2to4.R   (rumen-r env; writes results/figures_manuscript/fig3_calibration,
# fig4_clade_ratio_trait and fig5_wallace_vs_difford, each .pdf and .png). Figure numbers follow the order the manuscript cites them.
# Journal export: 170 mm (6.69 in) wide, PNG at 300 dpi, no in-image title, subtitle or source line (these go in the
# manuscript figure legends). Panel headers a and b in Figure 3 are kept.
# Amendment 46.5 (reporting only): taxon names in italics (46.5.7); permutation p as (b + 1)/(B + 1) in Figure 4 (46.5.1); the Figure 4 pooling of
# Difford g per day with Wallace per kg intake labelled an exploratory comparison of different outcomes (46.5.5); the Figure 5 correlation shown
# without a p value, because genera are not independent (46.5.4).
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
W170 <- 170 / 25.4
mlab <- function(b) sub("^-", "\u2212", format(b, trim=TRUE, drop0trailing=TRUE))   # typographic minus on axes
save2 <- function(p, name, w, h) { ggsave(file.path(od, paste0(name, ".pdf")), p, width=w, height=h, device=cairo_pdf); ggsave(file.path(od, paste0(name, ".png")), p, width=w, height=h, dpi=300, bg=bg); cat("wrote", name, "\n") }
tk <- rd(file.path(root, "plan", "target_keys.tsv"))
repel <- requireNamespace("ggrepel", quietly=TRUE)
# sa: no repel segment for blank (repel-only) labels
lab_layer <- function(data, ...) if (repel) ggrepel::geom_text_repel(data=transform(data, sa=ifelse(nzchar(genus), 1, 0)), ..., parse=TRUE, size=2.3, colour=ink2, min.segment.length=0.2, seed=1, max.overlaps=Inf) else geom_text(data=data, ..., parse=TRUE, size=2.3, colour=ink2, vjust=-0.7, check_overlap=TRUE)
# plotmath for a taxon name: the Latin part in italics, SILVA placeholder suffixes ("UCG-001", "group") and codes in roman (46.5.7)
q <- function(x) paste0("'", gsub("'", "\\\\'", x), "'")
tax_expr <- function(g, suffix="") vapply(seq_along(g), function(i) { x <- g[i]; if (!nzchar(x)) return("")   # empty label: no text and no repel segment
  pp <- strsplit(sub(" \\((.*)\\)$", "_(\\1)", x), "[_ ]")[[1]]; lat <- grepl("^[A-Z][a-z]+$", pp[1]) && !grepl("^(UCG|CAG|NK|FD|NED)", pp[1])
  rest <- paste(pp[-1], collapse=" "); sfx <- suffix[min(i, length(suffix))]
  if (lat && !nzchar(rest) && startsWith(sfx, ",")) return(paste0("italic(", q(paste0(pp[1], ",")), ")*", q(substring(sfx, 2))))   # comma inside the italic run: no italic-correction gap
  e <- if (lat) paste0("italic(", q(pp[1]), ")") else q(pp[1])
  if (nzchar(rest)) e <- paste0(e, "*", q(paste0(" ", rest))); if (nzchar(suffix[min(i, length(suffix))])) e <- paste0(e, "*", q(suffix[min(i, length(suffix))])); e }, "")

# ---------------- Figure 3
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
  geom_text(data=sN, aes(lab, top + 0.035, label=paste0("max ", trimws(format(maxq, big.mark=",")), "\nat q < 0.05")), size=2.25, colour=ink2, lineheight=0.9) +
  annotate("text", x=5.45, y=0.05, label="expected 5%", hjust=1, vjust=-0.6, size=2.25, colour=ink2) +
  scale_y_continuous(labels=function(v) paste0(round(100 * v), "%"), expand=expansion(mult=c(0.02, 0.18))) +
  labs(title="a  False positives across 200 shuffles", x=NULL, y="Genera at nominal p < 0.05") +
  theme_f + theme(panel.grid.major.x=element_blank(), axis.text.x=element_text(size=6.2, lineheight=0.9))
nsd <- function(dir, fmt, tab, method) { fs <- file.path(root, "results", dir, sprintf(fmt, tab, 1001:1200)); fs <- fs[file.exists(fs)]; stopifnot(length(fs) == 200)
  x <- do.call(rbind, lapply(fs, function(f) { d <- rd(f); d[d$status == "ok", c("id", "W")] })); a <- aggregate(W ~ id, x, sd); names(a)[2] <- "sd"; a$table <- tab; a$method <- method; a }
B <- rbind(nsd("methane_dryrun", "%s_seed%d_targets.tsv", "difford_bact", "ANCOM-BC2"), nsd("methane_dryrun", "%s_seed%d_targets.tsv", "wallace_prokaryotes", "ANCOM-BC2"),
           nsd("clr", "%s_seed%d_clr_targets.tsv", "difford_bact", "CLR linear model"), nsd("clr", "%s_seed%d_clr_targets.tsv", "wallace_prokaryotes", "CLR linear model"))
B <- B[B$id != "M1", ]; B$key <- paste(B$table, B$id); B <- B[B$key %in% B$key[B$method == "CLR linear model"] & B$key %in% B$key[B$method == "ANCOM-BC2"], ]
B$lab <- paste0(tk$target[match(B$id, tk$id)], ", ", ifelse(B$table == "difford_bact", "Difford", "Wallace"))
B$lab_e <- tax_expr(tk$target[match(B$id, tk$id)], paste0(", ", ifelse(B$table == "difford_bact", "Difford", "Wallace")))
B$lab <- factor(B$lab, B %>% filter(method == "ANCOM-BC2") %>% arrange(sd) %>% pull(lab))
seg <- B %>% group_by(lab) %>% summarise(lo=min(sd), hi=max(sd), .groups="drop")
pB <- ggplot(B, aes(sd, lab)) +
  geom_vline(xintercept=1, colour=ink2, linetype="22", linewidth=0.35) +
  geom_segment(data=seg, aes(x=lo, xend=hi, y=lab, yend=lab), colour=grid, linewidth=1.3, inherit.aes=FALSE) +
  geom_point(aes(colour=method), size=1.9) +
  scale_colour_manual(values=c("ANCOM-BC2"=blue, "CLR linear model"=orange)) +
  scale_y_discrete(labels=function(b) parse(text=B$lab_e[match(b, B$lab)])) +
  labs(title="b  Null spread per methane target", x="SD of W across 200 shuffles\n(1 = calibrated)", y=NULL) +
  theme_f + theme(panel.grid.major.y=element_blank(), legend.justification=c(1, 0), legend.direction="vertical",
                  legend.background=element_rect(fill=NA, colour=NA)) +
  # legend in the panel's lower right so the 170 mm width fits: rows are sorted by ANCOM-BC2 null SD, so the bottom rows sit near
  # x = 1 (all CLR SDs are 0.94 to 1.08) while the axis runs to about 5. legend.position.inside needs ggplot2 >= 3.5.
  if (packageVersion("ggplot2") >= "3.5.0") theme(legend.position="inside", legend.position.inside=c(0.99, 0.02)) else theme(legend.position=c(0.99, 0.02))
f2 <- (pA | pB) + plot_layout(widths=c(1.35, 1)) + plot_annotation(theme=theme_f)   # title, subtitle and sources moved to the manuscript legend
save2(f2, "fig3_calibration", W170, 3.6)

# ---------------- Figure 4 (clade ratio)
coh <- function(dirn, lab, cohort, trait) { x <- rd(file.path(root, "results", dirn, "estimates_calibrated.tsv")); x <- x[x$id == "M1" & x$cohort == cohort, ]; stopifnot(nrow(x) == 1)
  data.frame(row=lab, r=x$r, lo=x$lo, hi=x$hi, n=x$n, p=(round(x$p_perm * x$k_null) + 1) / (x$k_null + 1), kind="Single cohort", trait=trait) }   # (b + 1)/(B + 1)
pool <- function(dirn, lab, trait) { x <- rd(file.path(root, "results", dirn, "classification_primary.tsv")); x <- x[x$id == "M1", ]; stopifnot(nrow(x) == 1)
  data.frame(row=lab, r=x$r, lo=x$ci_lo, hi=x$ci_hi, n=NA, p=x$p, kind="Pooled", trait=trait) }
M <- rbind(coh("methane_pooled", "Difford, g per day", "difford", "g per day"),
           coh("methane_pooled", "Wallace, g per day", "wallace", "g per day"),
           coh("methane_pooled_noSE1", "Wallace without SE1, g per day", "wallace", "g per day"),
           coh("methane_pooled_depth", "Wallace, depth-adjusted, g per day", "wallace", "g per day"),
           coh("methane_pooled_ch4_dmi", "Wallace, per kg intake", "wallace", "per kg intake"),
           pool("methane_pooled", "Pooled, g per day (primary)", "g per day"),
           pool("methane_pooled_ch4_dmi", "Exploratory: pooled Difford g per day\nwith Wallace per kg intake (different outcomes)", "per kg intake"))
M$row <- factor(M$row, rev(M$row)); M$plab <- ifelse(M$kind == "Pooled", "", ifelse(M$p < 0.01, sprintf("p = %.3f", M$p), sprintf("p = %.2f", M$p)))
p_ratio <- ggplot(M, aes(r, row)) +
  geom_vline(xintercept=0, colour=ink2, linetype="22", linewidth=0.35) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, colour=trait), height=0, linewidth=1.1) +
  geom_point(aes(colour=trait, shape=kind), size=2.4, fill=bg, stroke=0.8) +
  geom_text(aes(label=plab), x=max(M$hi) + 0.03, hjust=0, size=2.35, colour=ink2) +
  scale_colour_manual(values=c("g per day"=blue, "per kg intake"=orange)) + scale_shape_manual(values=c("Single cohort"=16, "Pooled"=23)) +
  scale_x_continuous(expand=expansion(mult=c(0.05, 0.25)), labels=function(b) sub("^-", "\u2212", sprintf("%.2f", b))) +
  labs(x="Partial r with methane", y=NULL) + theme_f + theme(panel.grid.major.y=element_blank())
save2(p_ratio, "fig4_clade_ratio_trait", W170, 2.7)

# ---------------- Figure 5 (Wallace vs Difford)
w <- rd(file.path(root, "results", "denovo_scan", "wallace_prokaryotes_scan.tsv")); d <- rd(file.path(root, "results", "denovo_scan", "difford_bact_scan.tsv"))
g <- intersect(w$genus, d$genus); S <- data.frame(genus=g, x=w$W_cal[match(g, w$genus)], y=d$W_cal[match(g, d$genus)], xr=w$W[match(g, w$genus)], yr=d$W[match(g, d$genus)])
S <- S[is.finite(S$x) & is.finite(S$y), ]
fc <- rd(file.path(root, "results", "wallace_farms", "farm_consistency.tsv")); rob <- fc$genus[fc$farm_robust & fc$genus != "TGT_M2"]
S$robust <- S$genus %in% rob
e <- rd(file.path(root, "results", "methane_pooled", "estimates_calibrated.tsv"))
m2 <- data.frame(genus="Succinivibrionaceae (M2)", x=e$W[e$id == "M2" & e$cohort == "wallace"], y=e$W[e$id == "M2" & e$cohort == "difford"]); stopifnot(nrow(m2) == 1)
ct <- cor.test(S$x, S$y, method="spearman", exact=FALSE); ctr <- cor.test(S$xr, S$yr, method="spearman", exact=FALSE)
p_scatter <- ggplot(S, aes(x, y)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-1.96, ymax=1.96, fill="#f0efec") +
  geom_hline(yintercept=0, colour=ink2, linewidth=0.3) + geom_vline(xintercept=0, colour=ink2, linewidth=0.3) +
  geom_point(data=S[!S$robust, ], colour=grey, size=1.4) +
  geom_point(data=S[S$robust, ], colour=blue, size=1.9) +
  geom_point(data=m2, shape=23, fill=orange, colour=ink, size=2.8) +
  # Robust-genus labels repel from every plotted point (blank labels for the grey points and for the M2 diamond and its
  # label area), not only from each other. The M2 label sits at a fixed spot centred just below the diamond: repelled
  # placement had put it on the grey point at about (-4.5, 0).
  # blank labels along y = 0 keep the repelled labels off the axis line
  lab_layer(data=rbind(data.frame(x=S$x, y=S$y, genus=tax_expr(ifelse(S$robust, S$genus, ""))),
                       data.frame(x=m2$x + c(0, -1.2, 0, 1.2), y=m2$y + c(0, -0.25, -0.25, -0.25), genus=""),
                       data.frame(x=seq(min(S$x), max(S$x) * 1.15, length.out=80), y=0, genus="")), aes(label=genus, segment.alpha=sa), box.padding=0.5, point.padding=0.45, force=3, force_pull=0.4, max.iter=50000) +
  geom_text(data=transform(m2, genus=tax_expr(genus)), aes(label=genus), parse=TRUE, vjust=1, nudge_y=-0.2, size=2.3, colour=ink2) +
  annotate("text", x=min(c(S$x, m2$x)), y=2.85, hjust=0, vjust=1, size=2.5, colour=ink2, lineheight=1.05,   # hjust 0: both lines left-aligned
           label=sprintf("%d shared genera\nSpearman rho = %s (descriptive)", nrow(S), sub("^-", "\u2212", sprintf("%.2f", ct$estimate)))) +
  annotate("text", x=Inf, y=1.96, hjust=1.02, vjust=-0.5, size=2.25, colour=ink2, label="Difford |W| < 1.96") +
  labs(x="Calibrated W, Wallace (1,009 cows, 7 farms)", y="Calibrated W, Difford (749 cows, 6 herds)") + scale_y_continuous(limits=c(min(c(S$y, m2$y, -2)) - 0.3, 2.9), labels=mlab) +
  scale_x_continuous(expand=expansion(mult=c(0.04, 0.16)), labels=mlab) + theme_f   # room on the right for the crowded positive-W labels
save2(p_scatter, "fig5_wallace_vs_difford", W170, 5.0)

# ---------------- numbers for the manuscript text
cat("\n== Figure 3a: per-table false positive rate across shuffles\n"); print(as.data.frame(sN), digits=3, row.names=FALSE)
cat("== Figure 3b: null SD range, ANCOM-BC2", paste(round(range(B$sd[B$method == "ANCOM-BC2"]), 2), collapse=" to "), "| CLR", paste(round(range(B$sd[B$method == "CLR linear model"]), 2), collapse=" to "), "| targets", nlevels(B$lab), "\n")
cat("== Figure 5:", nrow(S), "genera; Spearman on calibrated W", round(ct$estimate, 3), "p", signif(ct$p.value, 2), "| on raw W", round(ctr$estimate, 3), "p", signif(ctr$p.value, 2), "| robust plotted", sum(S$robust), "| M2 W Wallace", round(m2$x, 2), "Difford", round(m2$y, 2), "\n")
cat("== Figure 4\n"); print(M[, c("row", "r", "lo", "hi", "n", "p")], digits=3, row.names=FALSE)
