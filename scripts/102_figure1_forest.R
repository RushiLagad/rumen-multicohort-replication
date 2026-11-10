#!/usr/bin/env Rscript
# 102: Manuscript Figure 1. Pooled r and 95% CI for the 27 feed efficiency and 11 methane targets, oriented to the published
# direction, with r = 0.15 and r = 0.30 bands shaded. Reads results/power/power_equivalence.tsv (from 93); refits nothing.
# Run: Rscript scripts/102_figure1_forest.R      (rumen-r env; writes results/figures_manuscript/fig1_forest_pooled.{pdf,png})
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "figures_manuscript"); dir.create(od, showWarnings=FALSE)
p <- read.delim(file.path(root, "results", "power", "power_equivalence.tsv"), check.names=FALSE, stringsAsFactors=FALSE, quote="")
stopifnot(sum(p$arm == "RFI") == 27, sum(p$arm == "Methane") == 11, !any(grepl("^replicated", p$label)))

ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; bg <- "#fcfcfb"
cat_lv <- c("Excludes any effect in the published direction", "Excludes r ≥ 0.15",
            "Excludes r ≥ 0.30", "Underpowered: r = 0.30 not excluded")
vcol <- setNames(c("#104281", "#256abf", "#6da7ec", "#b4b2a9"), cat_lv)   # ordinal blue ramp = strength of exclusion; grey = no conclusion
short <- c("Methanobrevibacter gottschalkii clade vs ruminantium clade" = "M. gottschalkii : M. ruminantium clade ratio",
           "unclassified Methanomassiliicoccaceae" = "Unclassified Methanomassiliicoccaceae")

d <- p %>% mutate(
  trait = ifelse(arm == "RFI", "RFI", "CH4"),
  directional = !is.na(published_dir), dir = ifelse(directional, published_dir, 1),
  x = r * dir, xl = pmin(ci_lo * dir, ci_hi * dir), xh = pmax(ci_lo * dir, ci_hi * dir),
  verdict = factor(case_when(grepl("^CI excludes", statement) ~ cat_lv[1], grepl("0[.]15", statement) ~ cat_lv[2],
                             grepl("0[.]30 in", statement) ~ cat_lv[3], TRUE ~ cat_lv[4]), cat_lv),
  name = ifelse(target %in% names(short), short[target], target),
  lab = sprintf("%s  %s", id, name),
  pub = ifelse(directional, sprintf("%s %s", ifelse(published_dir > 0, "higher", "lower"), trait), "not stated"),
  neff = format(n_eff, big.mark=","))
stopifnot(!anyNA(d$verdict), sum(d$arm == "RFI" & d$verdict != cat_lv[4]) == 17, sum(d$arm == "RFI" & d$verdict == cat_lv[4]) == 10,
          all(d$verdict[d$arm == "Methane"] %in% cat_lv[1:2]))

theme_f <- theme_minimal(base_size=8.5) + theme(
  text = element_text(colour=ink), axis.text = element_text(colour=ink2), axis.title = element_text(colour=ink2, size=8),
  panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour=grid, linewidth=0.3),
  legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(colour=ink2, size=7.5),
  plot.title = element_text(face="bold", size=10), plot.subtitle = element_text(colour=ink2, size=8), plot.caption = element_text(colour=muted, size=6.5, hjust=0),
  plot.title.position = "plot", plot.caption.position = "plot", plot.background = element_rect(fill=bg, colour=NA))
XL <- c(-0.7, 0.7); xt <- 0.74   # text columns sit just right of the panel

panel <- function(a, title, first, last) {
  s <- d %>% filter(arm == a); s$lab <- factor(s$lab, s %>% arrange(desc(xh)) %>% pull(lab))
  ggplot(s, aes(x, lab)) +
    annotate("rect", xmin=0.15, xmax=0.30, ymin=-Inf, ymax=Inf, fill="#f0efec") +
    annotate("rect", xmin=0.30, xmax=XL[2], ymin=-Inf, ymax=Inf, fill="#e4e2dc") +
    geom_vline(xintercept=0, colour=ink2, linewidth=0.35, linetype="22") +
    geom_errorbarh(aes(xmin=xl, xmax=xh, colour=verdict), height=0, linewidth=1.1) +
    geom_point(aes(colour=verdict, shape=directional), size=1.9, stroke=0.7, fill=bg) +
    geom_text(aes(x=xt, label=pub), hjust=0, size=2.35, colour=ink2) +
    geom_text(aes(x=xt + 0.40, label=neff), hjust=1, size=2.35, colour=ink2) +
    { if (first) annotate("text", x=c(xt, xt + 0.40), y=Inf, label=c("Published", "Eff. n"), hjust=c(0, 1), vjust=-0.6,
                          size=2.35, fontface="bold", colour=ink2) } +
    { if (first) annotate("text", x=c(0.225, 0.50), y=Inf, label=c("0.15–0.30", "≥ 0.30, typical\npublished size"),
                          vjust=-0.25, size=2.3, colour=ink2, lineheight=0.9) } +
    scale_colour_manual(values=vcol, limits=cat_lv, drop=FALSE) +
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21), guide="none") +
    scale_x_continuous(breaks=seq(-0.6, 0.6, 0.3), expand=expansion(0)) + scale_y_discrete(expand=expansion(add=c(0.7, if (first) 1.3 else 0.7))) +
    coord_cartesian(xlim=XL, clip="off") +
    labs(title=title, y=NULL, x=if (last) "Pooled partial r, oriented so positive = the published direction" else NULL) +
    theme_f + theme(plot.margin=margin(if (first) 22 else 10, 112, 4, 4), legend.box.spacing=unit(16, "pt"), axis.text.x=if (last) element_text() else element_blank(),
                    plot.title=element_text(size=8.5, face="bold"), legend.position=if (first) "top" else "none") +
    guides(colour=guide_legend(ncol=2, byrow=FALSE))
}
f1 <- panel("RFI", "a  Feed efficiency (RFI): 27 targets, three beef cohorts, 649 animals", TRUE, FALSE) /
      panel("Methane", "b  Methane (g per day): 11 targets, two dairy cohorts, 1,758 cows", FALSE, TRUE) +
  plot_layout(heights=c(27, 11)) +
  plot_annotation(
    title="No published association replicated; most published-size effects are excluded",
    subtitle=paste0("Calibrated pooled partial correlation and 95% CI per target (REML random effects; single-cohort targets show that cohort).\n",
                    "Shaded: effects of r = 0.15 to 0.30 and of r ≥ 0.30 in the published direction. Hollow points: no published direction,\n",
                    "shown on the estimated scale. Eff. n: effective sample size implied by the pooled interval."),
    caption=paste0("Feed efficiency: Nelore (52), Charolais (30), Alberta (567, all animals). Methane: Difford (749) and Wallace (1,009); Wallace is\n",
                   "re-analysis for Succinivibrionaceae (M2) and enters it only in a sensitivity pooling. Source: results/power/power_equivalence.tsv."),
    theme=theme_f)
ggsave(file.path(od, "fig1_forest_pooled.pdf"), f1, width=7.2, height=8.6, device=cairo_pdf)
ggsave(file.path(od, "fig1_forest_pooled.png"), f1, width=7.2, height=8.6, dpi=300, bg=bg)
cat("wrote", file.path(od, "fig1_forest_pooled.{pdf,png}"), "\n")
