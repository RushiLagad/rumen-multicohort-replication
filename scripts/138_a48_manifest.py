#!/usr/bin/env python3
"""Amendment 48: numbers-only manifest of the Amendment 48 results (late analysis, run after all primary, late,
Amendment 46 and Amendment 47 results were seen). Reads results/late_a48/ (no new fit) and writes
results/late_a48/manifest_a48.tsv with one row per number: section, item, file, column, row_key, value, rounding.
Draft 18 is checked against this file. In clade_ratio/grid.tsv the analysis is decoded from the seed
(seed = 4890000 + 100 x analysis + cell; see the Amendment 48 result, item 1).
Usage: 138_a48_manifest.py [root]"""
import os, sys
import pandas as pd

root = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"
d = os.path.join(root, "results", "late_a48"); out = os.path.join(d, "manifest_a48.tsv")
if os.path.exists(out):
    sys.exit(f"output exists, Amendment 48 never overwrites: {out}")
rd = lambda f: pd.read_csv(os.path.join(d, f), sep="\t")
R = []

def fmt(v, dp):
    if isinstance(v, str):
        return v
    if pd.isna(v):
        return "NA"
    if dp is None:
        try:
            return str(int(v)) if float(v).is_integer() and not isinstance(v, bool) and str(v) not in ("True", "False") else str(v)
        except (TypeError, ValueError):
            return str(v)
    return f"{float(v):.{dp}f}"

def add(sec, item, f, col, key, v, dp=3):
    num = not isinstance(v, (bool, str)) and not pd.isna(v) and dp is not None
    try:
        import numpy as np
        num = num and not isinstance(v, np.bool_)
    except ImportError:
        pass
    R.append((sec, item, f, col, key, fmt(v, dp if num else None), f"{dp} dp" if num else "as stored"))

def rows(sec, f, df, keycols, cols, dp=3, label=None):
    for _, x in df.iterrows():
        key = ";".join(f"{k}={x[k]}" for k in keycols)
        for c in cols:
            add(sec, label(x) if label else key, f, c, key, x[c], dp if c not in INTCOLS else None)

INTCOLS = {"n", "df", "k", "cohorts", "shuffles", "n_genera", "n_bh", "n_robust", "median_n_robust", "max_n_robust", "max_n_bh",
           "genera", "same_sign", "wallace_q05", "difford_q05", "features", "estimable_pi", "consistent_across_farms", "cells",
           "refits", "refits_holding", "most_influential_position", "matched", "k_ge2", "k1", "replicated_a37", "Sa_clr",
           "Sa_ancombc2", "Sb_clr", "Sb_ancombc2", "taxa_flagged", "n_taxa", "n_q05", "cows", "removed", "k_farms",
           "farms_same_sign", "asvs", "asvs_assigned", "targets_testable", "rfi_excl_30", "methane_excl_15", "replicated",
           "cohort_tests", "cells_testable", "cells_interval_excludes_0", "influence_refits_excluding_0", "tables_meeting", "tables"}

# ---- 48.1 seven claims
f = "seven/seven_claims.tsv"; S = rd(f)
rows("48.1 per cohort", f, S, ["id", "cohort"], ["n", "df", "r_clr", "lo_clr", "hi_clr", "p_t_clr", "p_fl_clr", "W_ancombc2_cal",
     "r_ancombc2", "lo_ancombc2", "hi_ancombc2", "p_perm_ancombc2", "mde_80", "same_sign_clr", "same_sign_ancombc2",
     "excl0_pub_clr", "excl0_pub_ancombc2"])
f = "seven/seven_pooled.tsv"; P = rd(f)
rows("48.1 pooled", f, P, ["id"], ["k", "r_mkh", "lo_mkh", "hi_mkh", "p_mkh", "r_fe", "lo_fe", "hi_fe", "p_fe_clr", "tau_hat", "I2", "Q_p",
     "r_a37", "lo_a37", "hi_a37", "q_a37", "label_a37", "label_a47", "reasons_not_replicated", "pattern_clr", "pattern_ancombc2",
     "estimators_agree_pattern", "secondary_Sa_clr", "secondary_Sa_ancombc2", "secondary_Sb_clr", "secondary_Sb_ancombc2", "mde_fe"])
f = "seven/concordance.tsv"; rows("48.1 concordance", f, rd(f), ["set", "estimator"], ["n", "same", "share", "p_binom"])
add("48.1 concordance", "seven claims: published sign under both estimators", f, "overlap_note", "set=seven claims (S1 matched)",
    rd(f).overlap_note.dropna().iloc[0], None)
f = "seven/counts_by_classification.tsv"; rows("48.1 counts", f, rd(f), ["classification"],
     ["matched", "k_ge2", "k1", "replicated_a37", "Sa_clr", "Sa_ancombc2", "Sb_clr", "Sb_ancombc2", "Sa_clr_ids", "Sa_ancombc2_ids"])
f = "seven/s1_yield.tsv"; rows("48.1 S1 against methane yield", f, rd(f), ["id", "cohort", "trait"],
     ["n", "df", "r_clr", "lo_clr", "hi_clr", "p_t_clr", "r_ancombc2_cal", "lo_ancombc2_cal", "hi_ancombc2_cal", "p_perm_ancombc2"])

# ---- 48.2 ANCOM-BC2
f = "ancombc2/config.tsv"; C = rd(f)
add("48.2(a) audit", "ANCOMBC version", f, "ancombc_version", "all rows", C.ancombc_version.iloc[0], None)
for (s, a), x in C.groupby(["script", "argument"]):
    add("48.2(a) audit", f"{s} {a}", f, "resolved", f"script={s};argument={a}", x.resolved.iloc[0], None)
f = "ancombc2/label_p_check.tsv"; L = rd(f)
add("48.2(a) audit", "labels use the parametric p of the calibrated W (rows checked; rows agreeing)", f, "verdict", "all rows",
    f"{len(L)}; {(L.verdict == 'parametric p of calibrated W').sum()}", None)
f = "ancombc2/null_summary.tsv"; rows("48.2(b) filtered null", f, rd(f), ["table"], ["shuffles", "n_genera", "median_share_p05",
     "min_share_p05", "max_share_p05", "median_share_p05_passed", "min_share_p05_passed", "max_share_p05_passed", "median_share_passed",
     "median_n_bh", "max_n_bh", "median_n_robust", "max_n_robust", "P_any_bh", "P_any_bh_lo", "P_any_bh_hi", "P_any_robust",
     "P_any_robust_lo", "P_any_robust_hi"])
f = "ancombc2/real_fits.tsv"; rows("48.2(b) real fits", f, rd(f), ["table"], ["n_genera", "share_p05", "share_p05_passed_all",
     "share_passed", "n_bh", "n_robust"])
f = "ancombc2/positive_controls.tsv"; rows("48.2(b) positive controls", f, rd(f), ["file", "term_column"], ["n_taxa", "n_q05", "n_robust"])
f = "ancombc2/charolais_struc_zero.tsv"; Z = rd(f)
add("48.2(c) structural zeros", "taxa flagged as structural zeros", f, "taxa_flagged", "all rows", Z.taxa_flagged.iloc[0], None)
add("48.2(c) structural zeros", "target features flagged", f, "flagged_structural_zero", "sum over targets", int(Z.flagged_structural_zero.sum()), None)
add("48.2(c) structural zeros", "largest |W change| over targets", f, "W_struc_zero - W_primary", "max over targets",
    (Z.W_struc_zero - Z.W_primary).abs().max(), 4)
add("48.2(c) structural zeros", "targets with q < 0.05 after refit", f, "q_struc_zero", "count < 0.05", int((Z.q_struc_zero < 0.05).sum()), None)
add("48.2(c) structural zeros", "primary call reproduced, largest |dW|", f, "primary_reproduced_max_dW", "all rows", Z.primary_reproduced_max_dW.max(), None)
f = "ancombc2/reading.tsv"; rows("48.2(d) reading", f, rd(f), ["rule"], ["tables_meeting", "tables", "tables_failing", "reading"])

# ---- 48.3 dairy
f = "dairy/agreement.tsv"; rows("48.3 agreement", f, rd(f), ["variant", "subset"], ["genera", "spearman_rho", "rho_lo", "rho_hi",
     "pearson_z", "same_sign", "share_same_sign", "p_sign_binom", "median_I2", "share_I2_gt50", "wallace_q05", "difford_q05"])
f = "dairy/reading.tsv"; rows("48.3 reading", f, rd(f), ["variant"], ["rho", "rho_lo", "rho_hi", "p_sign", "reading"])

# ---- 48.4 farms
f = "farms/reading.tsv"; rows("48.4 reading", f, rd(f), ["estimand"], ["features", "estimable_pi", "consistent_across_farms",
     "consistent_ids", "share_pi_includes_0", "reading"])
f = "farms/random_effects.tsv"; rows("48.4 random effects", f, rd(f), ["genus", "estimand"], ["k_farms", "method", "estimate",
     "lo_mkh", "hi_mkh", "p_mkh", "tau2", "I2", "Q_p", "pi_lo", "pi_hi", "pi_excludes_0", "outlying_farms", "W_full", "lofo_sign_kept"])

# ---- 48.5 CLR grid
f = "clr_grid/robustness.tsv"; rows("48.5 robustness", f, rd(f), ["id"], ["k", "statement_reference", "cells", "holds_in_all_cells",
     "failing_cells", "reading", "r_min", "r_max", "bound_min", "bound_max", "labels"])
f = "clr_grid/changes.tsv"; CH = rd(f)
for i, x in CH.groupby("id"):
    add("48.5 changes from pc0.5 R0", f"{i}: cells with any change", f, "cell", f"id={i}", ";".join(x.cell), None)
    for c in ["sign_changed", "excl15_changed", "excl30_changed", "excl0_changed", "label_changed"]:
        y = x[x[c]]
        if len(y):
            add("48.5 changes from pc0.5 R0", f"{i}: {c}", f, c, f"id={i}", ";".join(y.cell), None)
    for _, y in x[x.label_changed].iterrows():
        add("48.5 changes from pc0.5 R0", f"{i} {y.cell}: label", f, "label", f"id={i};cell={y.cell}", y.label, None)
f = "clr_grid/cohort_cells.tsv"; CC = rd(f)
for i in ["N06", "N13", "C3", "T4"]:
    x = CC[(CC.id == i) & (CC.cohort == "charolais_mrna")]
    ex = (x.lo > 0) | (x.hi < 0)
    add("48.5 Charolais intervals excluding 0", f"{i}: cells where the Charolais interval excludes 0", f, "lo;hi",
        f"id={i};cohort=charolais_mrna", f"{int(ex.sum())} of {len(x)}", None)
    add("48.5 Charolais intervals excluding 0", f"{i}: cells where it fails", f, "lo;hi", f"id={i};cohort=charolais_mrna",
        ";".join(x.cell[~ex]) if (~ex).any() else "none", None)

# ---- 48.6 Charolais
f = "charolais/summary.tsv"; rows("48.6 leave-one-animal-out", f, rd(f), ["id"], ["r_full", "lo_full", "hi_full", "bound_full", "statement",
     "r_min", "r_max", "bound_min", "bound_max", "refits", "refits_holding", "robust_to_single_animals", "most_influential_position",
     "max_abs_dfbetas", "grid_holds_in_all_cells", "grid_failing_cells", "at_threshold"])

# ---- 48.7 power
f = "power/mde_cohort.tsv"; rows("48.7 minimum detectable r", f, rd(f), ["cohort", "table"], ["n", "df", "mde_80"])
f = "power/mde_pooled.tsv"; rows("48.7 minimum detectable r, pooled", f, rd(f), ["id"], ["k", "dfs", "mde_fe", "mkh_t_crit",
     "mkh_expected_halfwidth_r_tau0"])
f = "power/published_effects.tsv"; PE = rd(f)
rows("48.7 published effects", f, PE, ["set", "id", "cohort"], ["source_value", "r_published", "r_clr", "lo_clr", "hi_clr", "comparison"])
for s, x in PE.groupby("set"):
    for c, n in x.comparison.value_counts().items():
        add("48.7 published effects, counts", f"{s}: {c}", f, "comparison", f"set={s}", int(n), None)
f = "power/attenuation.tsv"; A = rd(f)
rows("48.7 attenuation", f, A, ["quantity", "what", "trait_reliability", "feature_reliability"], ["r_observed", "r_true"])

# ---- 48.8 mapping
f = "mapping/summary.tsv"; rows("48.8 mapping", f, rd(f), ["set"], ["targets_testable", "rfi_excl_30", "methane_excl_15", "replicated",
     "cohort_tests", "classes"])
f = "mapping/mc1_pooled.tsv"; rows("48.8 MC1 only", f, rd(f), ["arm", "id"], ["k_full", "excl_15_full", "excl_30_full", "k", "cohorts",
     "r", "lo", "hi", "excl_0", "excl_15", "excl_30", "label"])

# ---- 48.9 clade ratio
NAMES = {1: "Wallace CH4 per kg DMI", 2: "Wallace CH4 g/d", 3: "Difford CH4 g/d", 4: "Wallace CH4 per kg DMI without SE1"}
f = "clade_ratio/grid.tsv"; G = rd(f)
G["analysis_no"] = (G.seed - 4890000) // 100
if not ((G.seed - 4890000) % 100 == G.cell).all():
    sys.exit("clade_ratio/grid.tsv: seed does not decode to the cell; stop")
G["analysis_name"] = G.analysis_no.map(NAMES)
rows("48.9 grid", f, G, ["analysis_name", "cell"], ["margin", "pseudo_count", "min_share", "share_assigned", "testable", "cows", "df",
     "t", "r", "lo", "hi", "p_t", "p_fl", "seed"])
for a, x in G.groupby("analysis_name"):
    ex = (x.lo > 0) | (x.hi < 0)
    add("48.9 grid summary", f"{a}: cells with interval excluding 0", f, "lo;hi", f"analysis={a}", f"{int(ex.sum())} of {len(x)}", None)
    add("48.9 grid summary", f"{a}: r range", f, "r", f"analysis={a}", f"{x.r.min():.3f} to {x.r.max():.3f}", None)
    add("48.9 grid summary", f"{a}: lower bound range", f, "lo", f"analysis={a}", f"{x.lo.min():.3f} to {x.lo.max():.3f}", None)
    add("48.9 grid summary", f"{a}: largest Freedman-Lane p", f, "p_fl", f"analysis={a}", x.p_fl.max(), 4)
f = "clade_ratio/influence.tsv"; rows("48.9 influence (Wallace CH4 per kg DMI)", f, rd(f), ["fit"], ["cows", "removed", "estimate", "t",
     "p", "r", "lo", "hi", "excludes_0"])
f = "clade_ratio/assignment.tsv"; rows("48.9 assignment", f, rd(f), ["cohort", "margin"], ["share_assigned", "asvs_assigned", "asvs"])
f = "clade_ratio/reading.tsv"; rows("48.9 reading (Wallace CH4 per kg DMI)", f, rd(f), ["cells"], ["cells_testable",
     "cells_interval_excludes_0", "influence_refits_excluding_0", "robust", "reading"])

M = pd.DataFrame(R, columns=["section", "item", "file", "column", "row_key", "value", "rounding"])
M.insert(0, "analysis", "late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen (Amendment 48)")
tmp = out + f".tmp{os.getpid()}"; M.to_csv(tmp, sep="\t", index=False); os.replace(tmp, out)
print(f"manifest_a48.tsv: {len(M)} rows, {M.section.nunique()} sections")
