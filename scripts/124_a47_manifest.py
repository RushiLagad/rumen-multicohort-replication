#!/usr/bin/env python3
"""Amendment 47: numbers-only manifest of the Amendment 47 results (late analysis, run after all primary, late and Amendment 46 results were seen).
Reads results/late_a47/ (no new fit) and writes results/late_a47/manifest_a47.tsv: section, item, value, source. Draft 13 is checked against this file.
Usage: 124_a47_manifest.py [root]"""
import os, sys
import pandas as pd

root = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"
d = os.path.join(root, "results", "late_a47"); out = os.path.join(d, "manifest_a47.tsv")
if os.path.exists(out):
    sys.exit(f"output exists, Amendment 47 never overwrites: {out}")
rd = lambda f: pd.read_csv(os.path.join(d, f), sep="\t")
R = []
def add(sec, item, val, src):
    R.append((sec, item, val, src))
f3 = lambda x: f"{x:.3f}"
iv = lambda lo, hi: f"{lo:.3f} to {hi:.3f}"

# ---- 47.1 validation
V = rd("intervals/validation.tsv"); cm = open(os.path.join(d, "intervals", "chosen_method.txt")).read().splitlines()
add("47.1 validation", "designated procedure", "modified Knapp-Hartung (REML, metafor test = adhoc)", "intervals/chosen_method.txt")
add("47.1 validation", "validated (tau <= 0.10, every rho, primary configurations and targets, coverage + 2 MCSE >= 0.93)", cm[1].split("validated: ")[1], "intervals/chosen_method.txt")
add("47.1 validation", "cells extended from 1,000 to 5,000 replicates", cm[3].split(", ")[1].split(" cells")[0], "intervals/chosen_method.txt")
for _, v in V.iterrows():
    add("47.1 validation", f"{v['name']}: rule cells passed, Check 1", f"{v.pass_check1} of {v.cells_check1}", "intervals/validation.tsv")
    add("47.1 validation", f"{v['name']}: rule cells passed, Check 2", f"{v.pass_check2} of {v.cells_check2}", "intervals/validation.tsv")
    add("47.1 validation", f"{v['name']}: lowest coverage (Check 1; Check 2)", f"{f3(v.min_coverage_check1)}; {f3(v.min_coverage_check2)}", "intervals/validation.tsv")
    add("47.1 validation", f"{v['name']}: worst Check 2 cell", v.worst_check2, "intervals/validation.tsv")
    add("47.1 validation", f"{v['name']}: passes the rule", "yes" if v.passes_rule else "no", "intervals/validation.tsv")
G = rd("intervals/generator_check.tsv")
add("47.1 generator", "corrected generator, largest |population r - rho| over all tables (tolerance 0.005)", f"{G.diff_corrected.abs().max():.4f}", "intervals/generator_check.tsv")
for c, x in G.groupby("cohort"):
    add("47.1 generator", f"Amendment 46 (46.2) generator, largest |population r - rho|, {c}", f"{x.diff_a46.abs().max():.4f}", "intervals/generator_check.tsv")
CO = rd("intervals/check2_cohort.tsv")
for c, x in CO.groupby("cohort"):
    m = x.groupby("rho").coverage.mean()
    add("47.1 per-cohort coverage (k = 1 intervals, corrected generator)", f"{c}: coverage averaged over targets, range over rho", f"{m.min():.3f} to {m.max():.3f}", "intervals/check2_cohort.tsv")
    add("47.1 per-cohort coverage (k = 1 intervals, corrected generator)", f"{c}: type I error at alpha 0.05 (rho = 0)", f3(x[x.rho == 0].rejection.mean()), "intervals/check2_cohort.tsv")
C2 = rd("intervals/check2_cells.tsv"); P = rd("intervals/pooled_intervals.tsv"); pp = P[P.version == "primary"]
for _, t in pp[pp.k >= 2].sort_values(["arm", "id"]).iterrows():
    x = C2[(C2.role == "primary") & (C2.id == t.id) & (C2.arm == t.arm)]
    for m in ["mkh", "hksj", "wald", "fe"]:
        y = x[(x.method == m) & (x.tau <= 0.10 + 1e-9)]
        add("47.1 per-target coverage (Check 2, lowest over rho, tau <= 0.10)", f"{t.arm} {t.id} {t.target}: {m}", f3(y.coverage.min()), "intervals/check2_cells.tsv")
    for tau in (0.20, 0.30):
        y = x[(x.method == "mkh") & ((x.tau - tau).abs() < 1e-9)]
        add("47.1 per-target coverage (Check 2, lowest over rho)", f"{t.arm} {t.id}: mkh at tau {tau:.2f}", f3(y.coverage.min()), "intervals/check2_cells.tsv")
    add("47.1 per-target statement", f"{t.arm} {t.id} {t.target}: tau-hat; tau_req; statement", f"{t.tau_hat:.3f}; {t.tau_req:.2f}; {t.statement}", "intervals/pooled_intervals.tsv")
CV = rd("intervals/configuration_validation.tsv")
for r, x in CV.groupby("role"):
    add("47.1 configurations (modified Knapp-Hartung)", f"{r}: pools passing at tau <= 0.10; <= 0.20; <= 0.30", f"{int(x.pass_tau010.sum())}; {int(x.pass_tau020.sum())}; {int(x.pass_tau030.sum())} of {len(x)}", "intervals/configuration_validation.tsv")

# ---- 47.1 real data: intervals and exclusion counts
for _, t in pp.sort_values(["arm", "k", "id"]).iterrows():
    add("47.1 intervals (primary)", f"{t.arm} {t.id} {t.target} (k = {t.k}): r; modified KH interval; HKSJ interval", f"{t.r:.3f}; {iv(t.lo_mkh, t.hi_mkh)}; {iv(t.lo_hksj, t.hi_hksj)}", "intervals/pooled_intervals.tsv")
EC = rd("intervals/exclusion_counts.tsv"); e = EC[(EC.version == "primary")]
for arm in ["RFI", "Methane"]:
    for kc in ["all", "single cohort (k = 1)", "pooled (k >= 2)"]:
        for cl in ["all", "directional", "magnitude"]:
            for m in ["mkh", "hksj"]:
                y = e[(e.arm == arm) & (e.k_class == kc) & (e.claim == cl) & (e.method == m)].iloc[0]
                add("47.1 exclusion counts (pointwise; BH does not adjust them)", f"{arm}, {kc}, {cl}, {m}: excludes r = 0.30; excludes r = 0.15; interval excludes 0 (of targets)",
                    f"{y.excl_r30}; {y.excl_r15}; {y.interval_excludes_0} (of {y.targets})", "intervals/exclusion_counts.tsv")
for arm, de in [("RFI", "30"), ("RFI", "15"), ("Methane", "15"), ("Methane", "30")]:
    x = pp[(pp.arm == arm) & pp[f"excl_{de}_mkh"]]
    add("47.1 exclusions by target (modified KH)", f"{arm}, excludes r = 0.{de}", "; ".join(f"{a} {b} (k = {k})" for a, b, k in zip(x.id, x.target, x.k)), "intervals/pooled_intervals.tsv")
x = pp[pp.excl0_mkh]; add("47.1 exclusions by target (modified KH)", "interval excludes 0 (association), both arms", "; ".join(f"{a} {b} {c} (k = {k})" for a, b, c, k in zip(x.arm, x.id, x.target, x.k)), "intervals/pooled_intervals.tsv")
x = pp[(pp.k >= 2) & (pp.excl_30_hksj | pp.excl_15_hksj | pp.excl0_hksj) & ~(pp.excl_30_mkh & pp.excl_15_mkh)]
add("47.1 exclusions by target (modified KH)", "pooled exclusions under HKSJ (Amendment 46) that modified KH loses at r = 0.30 or 0.15, or HKSJ intervals excluding 0",
    "; ".join(f"{a} {b}: HKSJ {iv(lo, hi)}, modified KH {iv(l2, h2)}" for a, b, lo, hi, l2, h2 in zip(x.id, x.target, x.lo_hksj, x.hi_hksj, x.lo_mkh, x.hi_mkh)), "intervals/pooled_intervals.tsv")

# ---- 47.2 labels
for f, arm in [("labels/labels_rfi.tsv", "RFI"), ("labels/labels_methane.tsv", "Methane")]:
    L = rd(f); L = L[~L.label_a47.str.startswith("untestable")]
    add("47.2 labels (CLR estimand, modified KH)", f"{arm}: testable targets", str(len(L)), f)
    for c in ["label_a47", "label_a47_fl_k1", "label_a42_clr", "label_a37_primary"]:
        add("47.2 labels (CLR estimand, modified KH)", f"{arm}: {c}", "; ".join(f"{k} {v}" for k, v in L[c].value_counts().sort_index().items()), f)
    dd = L[L.label_a47 != L.label_a37_primary]
    add("47.2 labels (CLR estimand, modified KH)", f"{arm}: targets whose label differs from Amendment 37", "; ".join(f"{a} {b}: {c} (A37 {e})" for a, b, c, e in zip(dd.id, dd.target, dd.label_a47, dd.label_a37_primary)) or "none", f)
    add("47.2 labels (CLR estimand, modified KH)", f"{arm}: replicated", str(int(L.label_a47.str.startswith("replicated").sum())), f)
    nv = L[L.label_note.fillna("") != ""]
    add("47.2 labels (CLR estimand, modified KH)", f"{arm}: labels resting on a configuration not validated", "; ".join(f"{a} ({n})" for a, n in zip(nv.id, nv.label_note)) or "none", f)
K = rd("labels/check_a42.tsv")
for _, k in K.iterrows():
    add("47.2 check", f"{k.arm}: REML-Wald pooling reproduces the Amendment 42 CLR labels", f"{k.labels_identical} of {k.targets}", "labels/check_a42.tsv")

# ---- 47.3 null scan
S = rd("null_scan/summary.tsv")
for _, s in S.iterrows():
    add("47.3 null scan (complete null)", f"{s.table}: P(any BH rejection at q < 0.05), Clopper-Pearson 95% CI", f"{s.any_rejection} of {s.replicates} = {s.p_any_rejection:.3f} ({s.cp_lo:.3f} to {s.cp_hi:.3f})", "null_scan/summary.tsv")
    add("47.3 null scan (complete null)", f"{s.table}: discovery counts (count:replicates); max", f"{s.discovery_counts}; {s.max_discoveries}", "null_scan/summary.tsv")
    add("47.3 null scan (complete null)", f"{s.table}: genera; inner permutations; mean share of null p < 0.05", f"{s.genera}; {s.inner_permutations}; {s.mean_frac_p05:.3f}", "null_scan/summary.tsv")
pj = open(os.path.join(d, "null_scan", "projection.txt")).read().splitlines()
add("47.3 null scan (complete null)", "extension to 1,000 replicates (CI containing both 0.05 and 0.10)", pj[-1].split(": ", 1)[1], "null_scan/projection.txt")

# ---- 47.4
W = rd("intervals/wallace_rfi_estimates.tsv"); add("47.4 four-cohort RFI", "Wallace RFI NRC (CLR): cows; df; targets", f"{W.n.iloc[0]}; {W.df.iloc[0]}; {len(W)}", "intervals/wallace_rfi_estimates.tsv")
P4 = rd("rfi4/pooled.tsv")
for _, t in P4.sort_values("id").iterrows():
    add("47.4 four-cohort RFI", f"{t.id} {t.target} (k = {t.k}): r; modified KH interval; excludes 0.30; excludes 0.15; statement",
        f"{t.r:.3f}; {iv(t.lo_mkh, t.hi_mkh)}; {'yes' if t.excl_30_mkh else 'no'}; {'yes' if t.excl_15_mkh else 'no'}; {t.statement}", "rfi4/pooled.tsv")
CS = rd("rfi4/canonical.tsv")
for _, c in CS.iterrows():
    add("47.4 four-cohort RFI", f"canonical C1 to C6, {c.method}: excludes r = 0.30; excludes r = 0.15 (of 6)", f"{int(c.excl_r30)}; {int(c.excl_r15) if pd.notna(c.excl_r15) else 'not computed'}", "rfi4/canonical.tsv")
L4 = rd("rfi4/labels.tsv"); add("47.4 four-cohort RFI", "replicated (four-cohort labels, modified KH)", f"{int(L4.label.str.startswith('replicated').sum())} of {len(L4)}", "rfi4/labels.tsv")
DS = rd("rfi4/depth_statement.tsv")
for _, x in DS.iterrows():
    add("47.4 depth adjustment", f"{x.estimator}: Wallace genera at q < 0.05, unadjusted to depth-adjusted", f"{x.wallace_q05_unadjusted} to {x.wallace_q05_depth_adjusted}", x.source)

M = pd.DataFrame(R, columns=["section", "item", "value", "source"]); M["analysis"] = "late analysis, run after all primary, late and Amendment 46 results were seen (Amendment 47)"
M.to_csv(out + ".tmp", sep="\t", index=False); os.replace(out + ".tmp", out); print(f"wrote {out}: {len(M)} rows")
