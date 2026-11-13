#!/usr/bin/env python3
"""119: results manifest (review item A6). One row per headline number used in the manuscript, figures, Additional files, README or cover letter,
computed from the committed result files (no fit, no hand-typed value) with the file and column it comes from and where it is used.
Draft text, figures and spreadsheets are checked against this file before submission. The manifest is an index, regenerated from the results
whenever they change (git keeps earlier versions); it changes no result.
Run: python3 scripts/119_manifest.py   (pandas, openpyxl). Writes results/manifest/headline_numbers.tsv and target_estimates.tsv.
"""
import os, re, subprocess
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort"); P = lambda *a: os.path.join(ROOT, *a)
od = P("results", "manifest"); os.makedirs(od, exist_ok=True)
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
A46 = "results/late_a46"
rows = []
def add(section, key, value, desc, src, used=""):
    if isinstance(value, (bool, np.bool_)): value = bool(value)
    elif isinstance(value, (np.integer,)): value = int(value)
    elif isinstance(value, (float, np.floating)): value = float(value)
    rows.append(dict(section=section, key=key, value=value, description=desc, source=src, used_in=used))
def fmt_ci(r, lo, hi): return f"{r:.3f} ({lo:.3f} to {hi:.3f})"

# ---------------------------------------------------------------- A. cohorts and animals
ce = rd(f"{A46}/effect_size/cohort_estimates.tsv"); ce_main = ce[~ce.table.str.endswith("_M1")]
n = ce_main.groupby("cohort").n.max()
lab = {"nelore": "Nelore", "charolais_mrna": "Charolais (mRNA)", "alberta_all": "Alberta (all animals)", "difford": "Difford", "wallace": "Wallace"}
for c in ["nelore", "charolais_mrna", "alberta_all", "difford", "wallace"]:
    add("A cohorts", f"n_{c}", int(n[c]), f"animals analysed, {lab[c]} (primary target tables)", f"{A46}/effect_size/cohort_estimates.tsv: max n per cohort", "Abstract; Methods; Fig. 1a; Fig. 2")
add("A cohorts", "n_rfi_animals", int(n[["nelore", "charolais_mrna", "alberta_all"]].sum()), "animals in the RFI arm", "sum of the above", "Fig. 2a title")
add("A cohorts", "n_methane_cows", int(n[["difford", "wallace"]].sum()), "cows in the methane arm", "sum of the above", "Fig. 2b title")
add("A cohorts", "n_total", int(n.sum()), "animals in all five cohorts ('about 2,400')", "sum of the above", "Abstract; Fig. 1a; Conclusions")
m1 = ce[ce.table.str.endswith("_M1")].set_index("cohort").n
for c in m1.index: add("A cohorts", f"n_M1_{c}", int(m1[c]), f"cows entering the M1 clade ratio, {c}", f"{A46}/effect_size/cohort_estimates.tsv (table *_M1)", "Methods (M1)")

# ---------------------------------------------------------------- B. targets
tk = rd("plan/target_keys.tsv")
add("B targets", "targets_prespecified", len(tk), "pre-specified targets, all sets", "plan/target_keys.tsv rows", "Fig. 1b; Additional file 1")
for s, g in tk.groupby("set"): add("B targets", f"targets_set_{s}", len(g), f"targets in set {s}", "plan/target_keys.tsv", "Fig. 1b; Methods")
pl = rd(f"{A46}/effect_size/pooled.tsv")
for a in ["RFI", "Methane"]:
    x = pl[pl.arm == a]
    add("B targets", f"testable_{a}", len(x), f"testable {a} targets (a tested feature in at least one independent cohort)", f"{A46}/effect_size/pooled.tsv rows", "Abstract; Results")
    for k, g in x.groupby("k"): add("B targets", f"k{k}_{a}", len(g), f"{a} targets tested in {k} independent cohort(s)", f"{A46}/effect_size/pooled.tsv: k", "Abstract ('14 and 6 ... one cohort'); Results; Fig. 2 k column")
led = rd(f"{A46}/ledger/evidence_ledger.tsv")
for (a, cat), g in led.groupby(["arm", "evidence_category"]):
    add("B targets", f"category_{a}_{cat.split(' (')[0].replace(' ', '_').lower()}", len(g), f"{a}: {cat}", f"{A46}/ledger/evidence_ledger.tsv: evidence_category", "Results; Additional file 1 (Evidence ledger)")
ov = rd(f"{A46}/ledger/overlap_audit.tsv")
for (role, ip), g in ov.groupby(["role", "in_primary_estimate"]):
    add("B targets", f"overlap_{role.replace(' ', '_').replace(',', '').lower()}", len(g), f"target-cohort pairs: {role} (in primary estimate: {ip})", f"{A46}/ledger/overlap_audit.tsv: role", "Methods; Additional file 1 (Overlap audit)")

# ---------------------------------------------------------------- C. replication labels
for a, f in (("RFI", "results/replication_rfi3_calibrated/classification_all.tsv"), ("Methane", "results/methane_pooled/classification_primary.tsv")):
    c = rd(f); add("C labels", f"replicated_{a}_primary", int((c.label == "replicated").sum()), f"{a} targets labelled replicated, primary (Amendment 37)", f"{f}: label", "Abstract; Results; Conclusions")
af2 = pd.read_excel(P("results", "additional_files", "Additional_file_2_estimates.xlsx"), sheet_name=None)
PLv = af2["Pooled results"]
add("C labels", "replicated_any_version", int((PLv.label == "replicated").sum()), "target-version rows labelled replicated, all versions in Additional file 2", "Additional_file_2: Pooled results, label", "Abstract; Discussion")
add("C labels", "versions", PLv.version.nunique(), "analysis versions in Additional file 2", "Additional_file_2: Pooled results, version", "Methods; Additional file 2")
ch = PLv[(PLv.label == "opposite sign")]
add("C labels", "label_changes_any_version", ";".join(f"{i} ({v})" for i, v in zip(ch.id, ch.version)) or "none", "targets labelled opposite sign in any version", "Additional_file_2: Pooled results, label", "Results; Discussion")

# ---------------------------------------------------------------- D. calibration (Figure 3)
lg = open(P("logs", "late_a46", "103_figs.log")).read()
fp = re.findall(r"^\s*(\w+)\\n.*?\s([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+(\d+)\s+200$", lg, re.M)
for nm, med, mean, top, maxq in fp:
    add("D calibration", f"shuffle_fpr_mean_{nm}", float(mean), f"{nm}: mean share of genera at nominal p < 0.05 across 200 shuffles (ANCOM-BC2)", "logs/late_a46/103_figs.log (Figure 3a table)", "Abstract ('6 to 26%'); Results; Fig. 3a")
    add("D calibration", f"shuffle_maxq_{nm}", int(maxq), f"{nm}: most genera at BH q < 0.05 in one shuffle", "logs/late_a46/103_figs.log", "Fig. 3a")
m = re.search(r"null SD range, ANCOM-BC2 ([0-9.]+) to ([0-9.]+) \| CLR ([0-9.]+) to ([0-9.]+)", lg)
add("D calibration", "null_sd_ancombc2_methane", f"{m.group(1)} to {m.group(2)}", "range of null SD of W, methane targets, ANCOM-BC2", "logs/late_a46/103_figs.log (Figure 3b)", "Results; Fig. 3b")
add("D calibration", "null_sd_clr_methane", f"{m.group(3)} to {m.group(4)}", "range of null SD of W, methane targets, CLR", "logs/late_a46/103_figs.log (Figure 3b)", "Results; Fig. 3b")
sp = rd("results/late_spikein/spikein_summary.tsv")
for _, r in sp[sp.b_lfc_per_sd == 0.25].iterrows():
    add("D calibration", f"spikein_fpr_{r.table}", f"{r.fpr_raw:.3f} raw, {r.fpr_calibrated:.3f} calibrated", f"spike-in, b = 0.25: false-positive rate among unplanted genera, {r.table}", "results/late_spikein/spikein_summary.tsv", "Results; Discussion ('about 25% to 5%')")
    add("D calibration", f"spikein_recovery_{r.table}", f"{r.power_calibrated:.3f} (median calibrated r {r.median_r_calibrated:.3f})", f"spike-in, b = 0.25: calibrated recovery of planted genera, {r.table}", "results/late_spikein/spikein_summary.tsv", "Discussion ('nine in ten ... r near 0.17')")
q = rd(f"{A46}/effect_size/a37_perm_p_recomputed.tsv")
add("D calibration", "perm_p_zero_primary_before", int((q.p_perm_a37 == 0).sum()), "primary per-cohort permutation p that were 0 as b/B", f"{A46}/effect_size/a37_perm_p_recomputed.tsv", "Methods (Amendment 46)")
add("D calibration", "perm_p_min_primary", float(q.p_perm_b1.min()), "smallest primary permutation p as (b + 1)/(B + 1)", f"{A46}/effect_size/a37_perm_p_recomputed.tsv", "Methods; Additional file 2")
add("D calibration", "perm_p_zero_additional_file2", int((af2["Per-cohort estimates"].p_permutation == 0).sum()), "permutation p equal to 0 in Additional file 2 (must be 0)", "Additional_file_2: Per-cohort estimates", "Additional file 2")

# ---------------------------------------------------------------- E. exclusion counts
ex = rd(f"{A46}/effect_size/exclusion_counts.tsv")
for _, r in ex.iterrows():
    g = r.group.replace(" ", "_")
    for col in ["excl_r30", "excl_r15", "excl_any_pubdir", "excl_r30_fe", "excl_r15_fe", "excl_r30_ancom", "excl_r15_ancom"]:
        add("E exclusions", f"{col}_{r.arm}_{g}", int(r[col]), f"{r.arm}, {r.group} ({r.targets} targets): {col} (CLR HKSJ unless _fe or _ancom)", f"{A46}/effect_size/exclusion_counts.tsv", "Abstract; Results; Discussion; Conclusions; Fig. 2")
def excl(v, a, thr):
    x = PLv[(PLv.version == v) & (PLv.k > 0)]; d = x.published_dir
    b = np.where(d.isna(), np.maximum(x.ci_lo_95.abs(), x.ci_hi_95.abs()), np.where(d > 0, x.ci_hi_95, -x.ci_lo_95))
    return int((b < thr).sum()), len(x)
for v, a, thr in (("RFI primary", "RFI", 0.30), ("RFI CLR", "RFI", 0.30), ("RFI LinDA (late)", "RFI", 0.30), ("RFI ALDEx2 (late)", "RFI", 0.30),
                  ("Methane primary", "Methane", 0.15), ("Methane CLR", "Methane", 0.15), ("Methane LinDA (late)", "Methane", 0.15), ("Methane ALDEx2 (late)", "Methane", 0.15)):
    if (PLv.version == v).any():
        k, nn = excl(v, a, thr); add("E exclusions", f"excl_{thr:.2f}_{v.replace(' ', '_')}", f"{k} of {nn}", f"{v}: Amendment 37-style intervals (W-derived r, REML-Wald) excluding r = {thr:.2f} in the published direction; comparison only", "Additional_file_2: Pooled results (ci_lo_95, ci_hi_95, published_dir)", "Results (estimator comparison)")

# ---------------------------------------------------------------- F. coverage (46.2)
cs = rd(f"{A46}/coverage_fix1/cohort_summary.tsv"); ps = rd(f"{A46}/coverage_fix1/pooled_summary.tsv")
add("F coverage", "coverage_cohort_range", f"{cs.coverage_mean.min():.3f} to {cs.coverage_mean.max():.3f}", "per-cohort coverage averaged over targets, all cohorts and rho", f"{A46}/coverage_fix1/cohort_summary.tsv", "Results ('94.4 to 96.6%')")
t1 = cs[cs.rho == 0].type1_mean; add("F coverage", "type1_range", f"{t1.min():.3f} to {t1.max():.3f}", "type I error of the per-cohort test at alpha 0.05", f"{A46}/coverage_fix1/cohort_summary.tsv", "Results")
for tau in (0, 0.05, 0.10):
    x = ps[ps.tau == tau]
    add("F coverage", f"coverage_hksj_tau{tau}", f"{x.coverage_hksj_mean.min():.3f} to {x.coverage_hksj_mean.max():.3f}", f"pooled HKSJ coverage, targets with k >= 2, tau = {tau}", f"{A46}/coverage_fix1/pooled_summary.tsv", "Results")
    add("F coverage", f"coverage_fe_tau{tau}", f"{x.coverage_fe_mean.min():.3f} to {x.coverage_fe_mean.max():.3f}", f"pooled fixed-effect coverage, tau = {tau}", f"{A46}/coverage_fix1/pooled_summary.tsv", "Results")
    add("F coverage", f"coverage_wald_tau{tau}", f"{x.coverage_reml_wald_mean.min():.3f} to {x.coverage_reml_wald_mean.max():.3f}", f"pooled REML-Wald coverage, tau = {tau}", f"{A46}/coverage_fix1/pooled_summary.tsv", "Results")
add("F coverage", "coverage_targets_k2plus", int(ps.targets.iloc[0]), "targets with k >= 2 in the pooled coverage check", f"{A46}/coverage_fix1/pooled_summary.tsv", "Results")
add("F coverage", "coverage_pass", open(P(A46, "coverage_fix1", "pass.txt")).read().split("46.2 overall: ")[1].split(" ")[0], "46.2 overall result", f"{A46}/coverage_fix1/pass.txt", "Methods; Results")

# ---------------------------------------------------------------- G. Freedman-Lane (46.3)
add("G Freedman-Lane", "clr_fl_p05", f"{int((ce.p_fl < 0.05).sum())} of {len(ce)}", "CLR per-cohort estimates with Freedman-Lane p < 0.05 (9,999 permutations)", f"{A46}/effect_size/cohort_estimates.tsv: p_fl", "Results")
add("G Freedman-Lane", "clr_fl_max_diff_vs_t", float((ce.p_fl - ce.p_t).abs().max()), "largest |Freedman-Lane p - parametric p|, CLR", f"{A46}/effect_size/cohort_estimates.tsv", "Methods; Results")
fl = rd(f"{A46}/effect_size/ancombc2_fl_methane.tsv").merge(q[q.arm == "methane"][["table", "id", "p_perm_b1"]], on=["table", "id"], how="left")
for _, r in fl.iterrows():
    add("G Freedman-Lane", f"ancombc2_fl_{r.table}_{r.id}", f"{r.p_fl:.3f} (within-design {r.p_perm_b1:.3f})" if pd.notna(r.p_perm_b1) else f"{r.p_fl:.3f}",
        f"ANCOM-BC2 Freedman-Lane p (200 permutations), {r.id} in {r.table}", f"{A46}/effect_size/ancombc2_fl_methane.tsv", "Results (Roseburia 0.11 to 0.02)")

# ---------------------------------------------------------------- H. scans
ss = rd(f"{A46}/scan/scan_summary.tsv")
for _, r in ss.iterrows():
    add("H scan", f"scan_{r.table}", f"{int(r.n_q05)} of {int(r.n_genera)} (min attainable {int(r.min_attainable_discoveries)}; {int(r.permutations)} permutations)", f"A46 scan, {r.table}: BH q < 0.05 discoveries", f"{A46}/scan/scan_summary.tsv", "Results; Methods")
    add("H scan", f"scan_null_{r.table}", f"mean {r.null_mean_q05:.2f}, max {int(r.null_max_q05)}", f"A46 scan null check (20 permuted traits), {r.table}", f"{A46}/scan/scan_summary.tsv", "Results")
wd = rd(f"{A46}/scan/wallace_difford_correlation.tsv")
for _, r in wd.iterrows():
    add("H scan", f"spearman_clr_{r.table_b}", f"{r.spearman_t:.3f} over {int(r.shared_genera)} genera", f"Wallace vs Difford, CLR scan t, descriptive", f"{A46}/scan/wallace_difford_correlation.tsv", "Results")
w = rd("results/denovo_scan/wallace_prokaryotes_scan.tsv"); d = rd("results/denovo_scan/difford_bact_scan.tsv"); g = sorted(set(w.genus) & set(d.genus))
x = w.set_index("genus").loc[g, "W_cal"]; y = d.set_index("genus").loc[g, "W_cal"]; ok = np.isfinite(x.values) & np.isfinite(y.values)
add("H scan", "spearman_ancombc2_calibrated", f"{spearmanr(x.values[ok], y.values[ok])[0]:.3f} over {int(ok.sum())} genera", "Wallace vs Difford, calibrated ANCOM-BC2 W, descriptive", "results/denovo_scan/*_scan.tsv: W_cal", "Abstract (-0.30); Fig. 5")

# ---------------------------------------------------------------- I. Wallace farms
fs = rd(f"{A46}/wallace_farms/summary.tsv")
add("I farms", "farm_genera", len(fs), "Wallace genera in the farm analysis (Amendment 40)", f"{A46}/wallace_farms/summary.tsv", "Results; README")
add("I farms", "farm_lofo_sign_kept", int(fs.lofo_sign_kept.sum()), "keep their sign in every leave-one-farm-out refit", f"{A46}/wallace_farms/summary.tsv", "Results; README")
add("I farms", "farm_lofo_stable", int(fs.lofo_stable_vs_full_data_threshold.sum()), "also above the full-data null 95th percentile ('stable', incl. Succinivibrionaceae M2)", f"{A46}/wallace_farms/summary.tsv", "Abstract ('16 genera and ... Succinivibrionaceae'); Results")
add("I farms", "farm_same_sign_all7", int((fs.farms_same_sign == 7).sum()), "same sign in all seven separately fitted farms", f"{A46}/wallace_farms/summary.tsv", "Abstract; Results; README")
add("I farms", "farm_same_sign_all_fitted", int((fs.farms_same_sign == fs.farms_fit).sum()), "same sign in every farm where fitted", f"{A46}/wallace_farms/summary.tsv", "Results")
add("I farms", "farm_I2_median", float(fs.I2.median()), "median I2 across farms (%)", f"{A46}/wallace_farms/summary.tsv", "Abstract; Results; README")
add("I farms", "farm_opposite_significant", int((fs.farms_interval_excludes_0_opposite_sign > 0).sum()), "genera with a farm significant in the opposite direction", f"{A46}/wallace_farms/summary.tsv", "Results")

# ---------------------------------------------------------------- J. key target estimates
pe = rd("results/power/power_equivalence.tsv")
for i in ["T4", "N13", "C6", "M1", "M2", "M3", "S1", "N30", "N32", "N38", "N40"]:
    a = pl[pl.id == i].iloc[0]; b = pe[(pe.id == i) & (pe.arm == a.arm)].iloc[0]
    add("J targets", f"est_{i}_clr", fmt_ci(a.r, a.lo, a.hi) + f"; k = {int(a.k)}; {a.interval}", f"{i} {a.target}: CLR partial r (A46 estimand)", f"{A46}/effect_size/pooled.tsv", "Results; Fig. 2")
    add("J targets", f"est_{i}_ancombc2", fmt_ci(b.r, b.ci_lo, b.ci_hi) + f"; q = {b.q_bh:.2f}", f"{i} {a.target}: calibrated ANCOM-BC2 r (Amendment 37, comparison)", "results/power/power_equivalence.tsv", "Results")
dm = rd("results/methane_pooled_ch4_dmi/estimates_calibrated.tsv"); r = dm[(dm.id == "M1") & (dm.cohort == "wallace")].iloc[0]
add("J targets", "est_M1_wallace_per_kg_dmi", fmt_ci(r.r, r.lo, r.hi) + f"; p = {(round(r.p_perm * r.k_null) + 1) / (r.k_null + 1):.3f}", "M1, Wallace, methane per kg DMI (sensitivity)", "results/methane_pooled_ch4_dmi/estimates_calibrated.tsv", "Abstract; Results; Fig. 4")
ns = rd("results/methane_pooled_noSE1/estimates_calibrated.tsv"); r = ns[(ns.id == "M1") & (ns.cohort == "wallace")].iloc[0]
add("J targets", "est_M1_wallace_noSE1", fmt_ci(r.r, r.lo, r.hi), "M1, Wallace without SE1 (overlap sensitivity)", "results/methane_pooled_noSE1/estimates_calibrated.tsv", "Results; Fig. 4; Additional file 1")

# ---------------------------------------------------------------- write
M = pd.DataFrame(rows)
try: head = subprocess.check_output(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], text=True).strip()
except Exception: head = "unknown"
M.insert(0, "repo_commit", head)
M.to_csv(os.path.join(od, "headline_numbers.tsv"), sep="\t", index=False)
T = pl.merge(pe[["arm", "id", "r", "ci_lo", "ci_hi", "q_bh", "label"]].rename(columns={"r": "r_ancombc2", "ci_lo": "lo_ancombc2", "ci_hi": "hi_ancombc2", "q_bh": "q_ancombc2", "label": "label_primary"}), on=["arm", "id"], how="left")
T = T.merge(led[["arm", "id", "evidence_category", "animals_total"]], on=["arm", "id"], how="left")
T = T[["arm", "id", "target", "evidence_category", "stratum", "published_dir", "k", "cohorts", "animals_total", "r", "lo", "hi", "interval", "excl_r30", "excl_r15", "excl_any",
       "r_ancombc2", "lo_ancombc2", "hi_ancombc2", "q_ancombc2", "label_primary"]]
T.insert(0, "repo_commit", head); T.to_csv(os.path.join(od, "target_estimates.tsv"), sep="\t", index=False)
print(f"manifest at commit {head}: {len(M)} headline numbers in {M.section.nunique()} sections; {len(T)} targets")
with pd.option_context("display.width", 250, "display.max_colwidth", 70, "display.max_rows", 500):
    print(M[["section", "key", "value"]].to_string(index=False))
