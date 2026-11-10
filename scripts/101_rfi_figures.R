#!/usr/bin/env Rscript
# 101: Figures for the multi-cohort feed-efficiency (RFI) replication, primary calibrated analysis: Alberta all-animal tables (Amendment 39 item 3), as used by 93.
# Reads only finished result tables; refits nothing.
#   Fig R1  forest plot: per-cohort calibrated partial r and pooled random-effects estimate, all 27 testable RFI targets
#   Fig R2  equivalence: pooled 95% CI oriented to the published direction against r = 0.15 and r = 0.30
#   Fig R3  robustness: replication label for every target under each analysis version
# Run: Rscript scripts/101_rfi_figures.R      (rumen-r env; writes results/figures_rfi/*.pdf and *.png)
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
root <- Sys.getenv("RUMEN_ROOT", "/data/rumen_multicohort"); od <- file.path(root, "results", "figures_rfi"); dir.create(od, showWarnings=FALSE)
rd <- function(f) read.delim(file.path(root, "results", f), check.names=FALSE, stringsAsFactors=FALSE, quote="")

cls <- rd("replication_rfi3_calibrated/classification_all.tsv") %>% filter(k > 0)
est <- rd("replication_rfi3_calibrated/estimates_calibrated_all.tsv") %>% filter(status == "ok", id %in% cls$id)
pow <- rd("power/power_equivalence.tsv") %>% filter(arm == "RFI")
rob <- rd("robustness/labels_wide_RFI.tsv")
stopifnot(nrow(cls) == 27, setequal(cls$id, pow$id), setequal(cls$id, rob$id), all(abs(pow$r[match(cls$id, pow$id)] - cls$r) < 1e-9))

# ---- shared styling (dataviz reference palette; validated light-mode, three slots pass CVD and normal-vision checks) ----
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"
n_of <- function(co) unique(est$n[est$cohort == co])
coh <- c(nelore = sprintf("Nelore, shotgun (n = %d)", n_of("nelore")),
         charolais_mrna = sprintf("Charolais, metatranscriptome (n = %d)", n_of("charolais_mrna")),
         alberta_all = sprintf("Alberta, 16S (n = %d)", n_of("alberta_all")))
col <- setNames(c("#2a78d6", "#eb6834", "#1baf7a"), coh); shp <- setNames(c(16, 15, 17), coh)  # shape = secondary encoding (aqua is < 3:1)
set_lab <- c(Canonical = "Canonical rumen taxa", Li = "Li et al. 2019, Microbiome 7:6", Nelore = "Nelore discovery taxa")
theme_r <- theme_minimal(base_size=8.5) + theme(
  text = element_text(colour=ink), axis.text = element_text(colour=ink2), axis.title = element_text(colour=ink2, size=8),
  panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour=grid, linewidth=0.3),
  strip.text = element_text(face="bold", hjust=0, colour=ink, size=8.5), strip.placement = "outside",
  legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(colour=ink2),
  plot.title = element_text(face="bold", size=10), plot.subtitle = element_text(colour=ink2, size=8), plot.caption = element_text(colour=muted, size=6.5, hjust=0),
  plot.background = element_rect(fill="#fcfcfb", colour=NA))
save2 <- function(p, stem, w, h) { ggsave(file.path(od, paste0(stem, ".pdf")), p, width=w, height=h, device=cairo_pdf)
  ggsave(file.path(od, paste0(stem, ".png")), p, width=w, height=h, dpi=300, bg="#fcfcfb"); cat("wrote", stem, "\n") }

pub_txt <- function(d) ifelse(is.na(d), "not stated", ifelse(d > 0, "higher RFI", "lower RFI"))
cls <- cls %>% mutate(set_f = factor(set_lab[set], set_lab),
  tlab = sprintf("%s  %s", id, target), pub = pub_txt(published_dir),
  ord = rank(-r, ties.method="first"))

# ---- Fig R1: forest ----------------------------------------------------------------------------------------------------
lv <- c("Pooled (random effects)", coh)
fp <- bind_rows(
  cls %>% transmute(id, set_f, tlab, src = lv[1], r, lo = ci_lo, hi = ci_hi),
  est %>% transmute(id, src = coh[cohort], r, lo, hi) %>% left_join(cls %>% select(id, set_f, tlab), by="id")) %>%
  filter(!(src == lv[1] & id %in% cls$id[cls$k == 1])) %>%      # single-cohort targets: the cohort estimate is the result
  mutate(src = factor(src, rev(lv)))
ylev <- cls %>% arrange(set_f, r) %>% pull(tlab); fp$tlab <- factor(fp$tlab, ylev); cls$tlab <- factor(cls$tlab, ylev)
dg <- position_dodge(width=0.75)
panel_r1 <- function(s, last, first) { d <- fp %>% filter(set_f == s) %>% droplevels(); cl <- cls %>% filter(set_f == s) %>% droplevels()
  ggplot(d, aes(r, tlab, group=src)) +
    geom_vline(xintercept=0, colour=ink2, linewidth=0.35, linetype="22") +
    geom_errorbarh(aes(xmin=lo, xmax=hi, colour=src), height=0, linewidth=0.45, position=dg) +
    geom_point(aes(colour=src, shape=src, size=src), position=dg, stroke=0.6) +
    geom_text(data=cl, aes(x=1.06, y=tlab, label=pub), inherit.aes=FALSE, hjust=0, size=2.4, colour=ink2) +
    annotate("text", x=1.06, y=Inf, label="Published direction", hjust=0, vjust=-0.5, size=2.4, fontface="bold", colour=ink2) +
    scale_colour_manual(values=c(setNames(ink, lv[1]), col), breaks=lv, limits=lv, drop=FALSE) +
    scale_shape_manual(values=c(setNames(18, lv[1]), shp), breaks=lv, limits=lv, drop=FALSE) +
    scale_size_manual(values=c(setNames(3.2, lv[1]), setNames(rep(1.5, 3), coh)), breaks=lv, limits=lv, drop=FALSE) +
    scale_x_continuous(breaks=seq(-1, 1, 0.5)) + scale_y_discrete(expand=expansion(add=0.6)) +
    coord_cartesian(xlim=c(-1, 1), clip="off") +
    labs(title=s, x=if (last) "Partial r with RFI (calibrated; positive = more abundant in less efficient animals)" else NULL, y=NULL) +
    theme_r + theme(plot.title=element_text(size=8.5, face="bold", hjust=0), plot.title.position="plot",
                    plot.margin=margin(4, 78, 4, 4), axis.text.x=if (last) element_text() else element_blank(),
                    legend.position=if (first) "top" else "none") +
    guides(colour=guide_legend(nrow=2, byrow=TRUE), shape=guide_legend(nrow=2, byrow=TRUE), size=guide_legend(nrow=2, byrow=TRUE)) }
sets <- levels(droplevels(cls$set_f))
r1 <- wrap_plots(lapply(seq_along(sets), function(i) panel_r1(sets[i], i == length(sets), i == 1)), ncol=1,
                 heights=as.numeric(table(droplevels(cls$set_f))[sets])) +
  plot_annotation(title="No published feed-efficiency association replicated across cohorts",
    subtitle=paste0("Per-cohort ANCOM-BC2 effects calibrated to the within-design permutation null (200 shuffles) and pooled by REML\n",
                    "random effects. Bars: 95% CI. Black diamond = pooled; single-cohort targets show only that cohort."),
    caption=paste0("Charolais contributes a converted binary-contrast effect; its extreme-group selection inflates it relative to a correlation.\n",
                   "Alberta: all phenotyped animals (Amendment 39 item 3). Source: results/replication_rfi3_calibrated/*_all.tsv."),
    theme=theme_r)
save2(r1, "figR1_forest_rfi", 8.2, 10.2)

# ---- Fig R2: equivalence ------------------------------------------------------------------------------------------------
cat_lv <- c("CI excludes any effect in the published direction", "Excludes r ≥ 0.15 in the published direction",
            "Excludes r ≥ 0.30 in the published direction", "Underpowered: r = 0.30 not excluded")
eq <- pow %>% left_join(cls %>% select(id, set_f, tlab), by="id") %>%
  mutate(dir = ifelse(is.na(published_dir), 1, published_dir), directional = !is.na(published_dir),
         x = r * dir, xl = pmin(ci_lo * dir, ci_hi * dir), xh = pmax(ci_lo * dir, ci_hi * dir),
         verdict = factor(case_when(grepl("^CI excludes", statement) ~ cat_lv[1], grepl("0[.]15", statement) ~ cat_lv[2],
                                    grepl("0[.]30 in", statement) ~ cat_lv[3], TRUE ~ cat_lv[4]), cat_lv),
         lab = sprintf("%s  %s%s", id, target, ifelse(directional, "", "  (direction not stated)")))
stopifnot(!anyNA(eq$verdict), sum(eq$verdict == cat_lv[4]) == 10)
eq$lab <- factor(eq$lab, eq %>% arrange(desc(xh)) %>% pull(lab))
vcol <- setNames(c("#104281", "#256abf", "#6da7ec", "#b4b2a9"), cat_lv)  # ordinal blue ramp for strength of exclusion; grey = no conclusion
r2 <- ggplot(eq, aes(x, lab)) +
  annotate("rect", xmin=0.30, xmax=Inf, ymin=-Inf, ymax=Inf, fill="#f0efec") +
  geom_vline(xintercept=0, colour=ink2, linewidth=0.35, linetype="22") +
  geom_vline(xintercept=c(0.15, 0.30), colour=muted, linewidth=0.3) +
  annotate("text", x=c(0.15, 0.30), y=Inf, label=c("r = 0.15", "r = 0.30 (typical published size)"), vjust=-0.5, hjust=c(1.05, -0.03), size=2.3, colour=ink2, lineheight=0.9) +
  geom_errorbarh(aes(xmin=xl, xmax=xh, colour=verdict), height=0, linewidth=1.1) +
  geom_point(aes(colour=verdict, shape=directional), size=1.9, stroke=0.7, fill="#fcfcfb") +
  scale_colour_manual(values=vcol, drop=FALSE) + scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21), guide="none") +
  scale_x_continuous(limits=c(-0.7, 0.7), breaks=seq(-0.6, 0.6, 0.3)) + coord_cartesian(clip="off") +
  labs(x="Pooled partial r, oriented so positive = the published direction", y=NULL,
       title="17 of 27 published RFI associations can be ruled out at their published size",
       subtitle="Pooled 95% CI per target. Shaded: effects as large as typically reported.\nHollow points: no published direction, shown on the estimated scale.",
       caption="Same pooled estimates as Fig R1 (results/power/power_equivalence.tsv).\nThe other 10 targets had too few animals to exclude r = 0.30.") +
  theme_r + theme(plot.margin = margin(22, 12, 6, 6), plot.title.position="plot", plot.caption.position="plot") + guides(colour=guide_legend(ncol=2, byrow=TRUE))
save2(r2, "figR2_equivalence_rfi", 7.6, 7.4)

# ---- Fig R3: robustness -----------------------------------------------------------------------------------------------
ver <- c(`replication_rfi3_calibrated:all` = "Calibrated\nall Alberta\n(primary)", `replication_rfi3_calibrated:primary` = "Calibrated\nAlberta\nsubset",
         `replication_rfi3:all` = "Uncalibrated\nall Alberta", `replication_rfi3:primary` = "Uncalibrated\nAlberta\nsubset", `clr_pooled:rfi_clr` = "CLR + linear\nmodel")
stopifnot(setequal(names(ver), setdiff(names(rob), c("id", "target"))))
rl <- rob %>% tidyr::pivot_longer(-c(id, target), names_to="version", values_to="label") %>%
  left_join(cls %>% select(id, set_f, tlab), by="id") %>% mutate(version = factor(ver[version], ver), tlab = factor(tlab, ylev),
  label = factor(ifelse(grepl("^replicated", label), "Replicated", ifelse(grepl("direction not stated", label), "Not replicated (no published direction)", "Not replicated")),
                 c("Replicated", "Not replicated", "Not replicated (no published direction)")))
r3 <- ggplot(rl, aes(version, tlab, fill=label)) + geom_tile(colour="#fcfcfb", linewidth=0.8) +
  scale_fill_manual(values=c(Replicated="#2a78d6", `Not replicated`="#898781", `Not replicated (no published direction)`="#c9c7be")) +
  scale_x_discrete(position="top") + facet_grid(set_f ~ ., scales="free_y", space="free_y") +
  labs(x=NULL, y=NULL, title="The verdict does not depend on the analysis version",
       subtitle="Replication label for every testable RFI target under each pooled version (Amendments 35-42).\nNo target is labelled replicated under any version.",
       caption="Calibrated = rescaled to the permutation null. All Alberta = every phenotyped animal (primary, Amendment 39);\nsubset = the earlier Alberta subset. CLR = second, independent estimator (Amendment 42).") +
  theme_r + theme(plot.title.position="plot", plot.caption.position="plot", panel.grid = element_blank(), axis.text.x.top = element_text(size=7, lineheight=0.9), strip.text.y = element_text(angle=0, hjust=0))
save2(r3, "figR3_robustness_rfi", 7.4, 8.4)
