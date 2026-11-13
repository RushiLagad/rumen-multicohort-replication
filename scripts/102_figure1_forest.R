#!/usr/bin/env Rscript
# 102: Manuscript Figure 2. Pooled CLR partial r and 95% interval for the 27 feed efficiency and 11 methane targets, oriented to the published
# direction, with the r = 0.15 to 0.30 and r >= 0.30 benchmark bands shaded and k (independent cohorts) beside each target.
# Amendment 46 (46.1, 46.5.7): the exclusion statements rest on the CLR partial correlation with Hartung-Knapp-Sidik-Jonkman intervals
# (k = 1: the cohort's own interval), whose coverage passed the 46.2 check, so the figure shows that estimand (results/late_a46/effect_size/
# pooled.tsv, from 112); the ANCOM-BC2-derived intervals that earlier versions of this figure showed are in Additional file 2. Refits nothing.
# The x axis is truncated at r = +/-0.6 so that narrow intervals stay readable; an interval reaching past the edge ends in an arrow with its true
# bound printed beside it (full bounds in Additional file 2).
# Run: Rscript scripts/102_figure1_forest.R      (rumen-r env; writes results/figures_manuscript/fig2_forest_pooled.{pdf,png})
# Journal export: 170 mm (6.69 in) wide, PNG at 300 dpi, no in-image title, subtitle or source line (these go in the manuscript figure legend).
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "figures_manuscript"); dir.create(od, showWarnings=FALSE)
p <- read.delim(file.path(root, "results", "late_a46", "effect_size", "pooled.tsv"), check.names=FALSE, stringsAsFactors=FALSE, quote="")
pass <- readLines(file.path(root, "results", "late_a46", "coverage_fix1", "pass.txt"))
stopifnot("46.2 coverage check must pass before exclusion classes are drawn" = any(grepl("^46.2 overall: PASS", pass)),
          sum(p$arm == "RFI") == 27, sum(p$arm == "Methane") == 11)

ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; bg <- "#fcfcfb"
cat_lv <- c("Excludes any effect in the published direction", "Excludes r ≥ 0.15",
            "Excludes r ≥ 0.30", "Not resolved: r = 0.30 inside the interval")
vcol <- setNames(c("#104281", "#256abf", "#6da7ec", "#b4b2a9"), cat_lv)   # ordinal blue ramp = strength of exclusion; grey = no conclusion
# y labels as plotmath: taxon names in italics (46.5.7)
q <- function(s) paste0("'", gsub("'", "\\\\'", s), "'")
lab_expr <- function(id, target) {
  if (id == "M1") return(paste0(q("M1  "), "*italic(", q("M. gottschalkii"), ")*", q(" : "), "*italic(", q("M. ruminantium"), ")*", q(" clade ratio")))
  if (grepl("^unclassified ", target)) return(paste0(q(paste0(id, "  unclassified ")), "*italic(", q(sub("^unclassified ", "", target)), ")"))
  if (target == "total archaea") return(q(paste0(id, "  total archaea")))
  paste0(q(paste0(id, "  ")), "*italic(", q(target), ")") }

d <- p %>% mutate(
  trait = ifelse(arm == "RFI", "RFI", "CH4"),
  directional = !is.na(published_dir), dir = ifelse(directional, published_dir, 1),
  x = r * dir, xl = pmin(lo * dir, hi * dir), xh = pmax(lo * dir, hi * dir),
  verdict = factor(case_when(excl_any ~ cat_lv[1], excl_r15 ~ cat_lv[2], excl_r30 ~ cat_lv[3], TRUE ~ cat_lv[4]), cat_lv),
  lab = id, expr = mapply(lab_expr, id, target),
  pub = ifelse(directional, sprintf("%s %s", ifelse(published_dir > 0, "higher", "lower"), trait), "not stated"))
stopifnot(!anyNA(d$verdict), sum(d$arm == "RFI" & d$verdict != cat_lv[4]) == 9, sum(d$arm == "Methane" & d$verdict != cat_lv[4]) == 7)

theme_f <- theme_minimal(base_size=8.5) + theme(
  text = element_text(colour=ink), axis.text = element_text(colour=ink2), axis.title = element_text(colour=ink2, size=8),
  panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour=grid, linewidth=0.3),
  legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(colour=ink2, size=7.5),
  plot.title = element_text(face="bold", size=10), plot.subtitle = element_text(colour=ink2, size=8), plot.caption = element_text(colour=muted, size=6.5, hjust=0),
  plot.title.position = "plot", plot.caption.position = "plot", plot.background = element_rect(fill=bg, colour=NA))
XL <- c(-0.6, 0.6); xt <- 0.63   # truncated axis; text columns sit just right of the panel
mnum <- function(v) sub("^-", "\u2212", sprintf("%.2f", v))
d <- d %>% mutate(over_hi = xh > XL[2], over_lo = xl < XL[1], xlp = pmax(xl, XL[1]), xhp = pmin(xh, XL[2]))

panel <- function(a, title, first, last) {
  s <- d %>% filter(arm == a); lv <- s %>% arrange(desc(xh), desc(x)) %>% pull(lab); s$lab <- factor(s$lab, lv); ex <- setNames(s$expr, s$lab)
  ggplot(s, aes(x, lab)) +
    annotate("rect", xmin=0.15, xmax=0.30, ymin=-Inf, ymax=Inf, fill="#f0efec") +
    annotate("rect", xmin=0.30, xmax=XL[2], ymin=-Inf, ymax=Inf, fill="#e4e2dc") +
    geom_vline(xintercept=0, colour=ink2, linewidth=0.35, linetype="22") +
    geom_errorbarh(aes(xmin=xlp, xmax=xhp, colour=verdict), height=0, linewidth=1.1) +
    geom_segment(data=filter(s, over_hi), aes(x=XL[2] - 0.04, xend=XL[2], y=lab, yend=lab, colour=verdict), linewidth=1.1,
                 arrow=arrow(length=unit(3.2, "pt"), type="closed"), show.legend=FALSE) +
    geom_segment(data=filter(s, over_lo), aes(x=XL[1] + 0.04, xend=XL[1], y=lab, yend=lab, colour=verdict), linewidth=1.1,
                 arrow=arrow(length=unit(3.2, "pt"), type="closed"), show.legend=FALSE) +
    geom_text(data=filter(s, over_hi), aes(x=XL[2] - 0.01, label=mnum(xh)), hjust=1, vjust=-0.55, size=2.0, colour=ink2) +
    geom_text(data=filter(s, over_lo), aes(x=XL[1] + 0.01, label=mnum(xl)), hjust=0, vjust=-0.55, size=2.0, colour=ink2) +
    geom_point(aes(colour=verdict, shape=directional), size=1.9, stroke=0.7, fill=bg) +
    geom_text(aes(x=xt, label=pub), hjust=0, size=2.35, colour=ink2) +
    geom_text(aes(x=xt + 0.30, label=k), hjust=1, size=2.35, colour=ink2) +
    { if (first) annotate("text", x=c(xt, xt + 0.30), y=Inf, label=c("Published", "k"), hjust=c(0, 1), vjust=-0.6,
                          size=2.35, fontface="bold", colour=ink2) } +
    { if (first) annotate("text", x=c(0.225, 0.45), y=Inf, label=c("0.15–\n0.30", "r ≥ 0.30\nbenchmark"),
                          vjust=-0.25, size=2.3, colour=ink2, lineheight=0.9) } +
    scale_colour_manual(values=vcol, limits=cat_lv, drop=FALSE, labels=sub("in the published", "in\nthe published", sub("inside the", "inside\nthe", cat_lv))) +
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21), guide="none") +
    scale_x_continuous(breaks=seq(-0.6, 0.6, 0.3), labels=function(b) sub("^-", "\u2212", format(b)), expand=expansion(0)) +
    scale_y_discrete(labels=function(b) parse(text=ex[b]), expand=expansion(add=c(0.7, if (first) 1.3 else 0.7))) +
    coord_cartesian(xlim=XL, clip="off") +
    labs(title=title, y=NULL, x=if (last) "Pooled CLR partial r, oriented so positive = the published direction" else NULL) +
    theme_f + theme(plot.margin=margin(if (first) 6 else 10, 92, 4, 4), legend.box.spacing=unit(16, "pt"), axis.text.x=if (last) element_text() else element_blank(),
                    plot.title=element_text(size=8.5, face="bold"), legend.position=if (first) "top" else "none") +
    guides(colour=guide_legend(ncol=2, byrow=FALSE))
}
f1 <- panel("RFI", "a  Feed efficiency (RFI): 27 targets, three beef cohorts, 649 animals", TRUE, FALSE) /
      panel("Methane", "b  Methane (g per day): 11 targets, two dairy cohorts, 1,758 cows", FALSE, TRUE) +
  plot_layout(heights=c(27, 11)) +
  plot_annotation(theme=theme_f)   # title, subtitle and source moved to the manuscript legend
W <- 170 / 25.4; H <- 7.9
ggsave(file.path(od, "fig2_forest_pooled.pdf"), f1, width=W, height=H, device=cairo_pdf)
ggsave(file.path(od, "fig2_forest_pooled.png"), f1, width=W, height=H, dpi=300, bg=bg)
cat("wrote", file.path(od, "fig2_forest_pooled.{pdf,png}"), "\n")
cat("== Figure 2 (CLR, HKSJ): classes by arm\n"); print(table(d$arm, d$verdict))
