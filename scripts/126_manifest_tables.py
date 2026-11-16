#!/usr/bin/env python3
"""126: table-level numbers for the manifest (review item A6, completing 119 and 124). One row per number in the manuscript that neither
results/manifest/headline_numbers.tsv (119) nor results/late_a47/manifest_a47.tsv (124) records: Table 1 (detection checks), Table 2 and the
calibration paragraph (200 within-design permutations per table), the depth paragraph, the estimator-dependence paragraphs, the LinDA and ALDEx2
null spread, the transfer positive controls, the spike-in and the effective n of the methane intervals. Each value is computed from the committed
result file named in its row (no fit, no hand-typed value). The three values recorded only in the plan (the depth correlations and the calibrated CLR
W of Eubacterium) are read from the committed plan text, and the script stops if that text is missing.
Run: python3 scripts/126_manifest_tables.py   (pandas). Writes results/manifest/table_numbers.tsv (never overwrites).
"""
import glob, os, re, subprocess, sys
import numpy as np, pandas as pd

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort"); P = lambda *a: os.path.join(ROOT, *a)
out = P("results", "manifest", "table_numbers.tsv")
if os.path.exists(out): sys.exit(f"output exists, never overwritten: {out}")
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
NSEED = int(os.environ.get("MANIFEST_TEST_SEEDS", "200")); TEST = NSEED != 200
rows = []
def add(section, key, value, desc, src, used=""):
    if isinstance(value, (np.integer,)): value = int(value)
    elif isinstance(value, (float, np.floating)): value = round(float(value), 4)
    rows.append(dict(section=section, key=key, value=value, description=desc, source=src, used_in=used))
def one(df, **kw):
    x = df
    for k, v in kw.items(): x = x[x[k] == v]
    assert len(x) == 1, (kw, len(x)); return x.iloc[0]

# ---------------------------------------------------------------- Table 1: detection checks
pc = rd("results/positive_controls/positive_controls_summary.tsv"); pc = pc[pc.profile == "profiles_gtdb_full"]
for coh, term in [("nelore", "diet"), ("charolais_mrna", "diet"), ("charolais_mrna", "year")]:
    r = one(pc, cohort=coh, term=term)
    add("T1 detection", f"{coh}_{term}_q05", int(r.n_q05), f"{coh}, {term} contrast: genera at BH q < 0.05", "results/positive_controls/positive_controls_summary.tsv (profiles_gtdb_full)", "Table 1")
    add("T1 detection", f"{coh}_{term}_genera", int(r.n_genera), f"{coh}, {term} contrast: genera tested", "results/positive_controls/positive_controls_summary.tsv (profiles_gtdb_full)", "Table 1")
for t in ["alberta", "difford_bact", "difford_arch", "wallace_prokaryotes"]:
    f = f"results/positive_controls/16s/{t}_summary.tsv"; r = rd(f).iloc[0]
    for k in ["n_animals", "n_genera", "n_q05"]: add("T1 detection", f"{t}_{k}", int(r[k]), f"{t}, {r.control} contrast: {k}", f, "Table 1")
    add("T1 detection", f"{t}_pass", bool(r["pass"]), f"{t}: detection check passed", f, "Table 1; Results")

# ---------------------------------------------------------------- Table 2 and the calibration paragraph
ss = rd("results/denovo_scan/scan_summary.tsv")
T2 = ["nelore_genus", "charolais_mrna_genus", "alberta_genus", "difford_bact", "wallace_prokaryotes"]
for t in T2:
    r = one(ss, table=t)
    add("T2 calibration", f"{t}_genera", int(r.n_genera), f"{t}: genera in the scan", "results/denovo_scan/scan_summary.tsv: n_genera", "Table 2")
    add("T2 calibration", f"{t}_fpr", float(r.null_fpr_p05), f"{t}: share of genera at nominal p < 0.05 under permutation, mean over 200 shuffles", "results/denovo_scan/scan_summary.tsv: null_fpr_p05", "Table 2; Abstract (6 to 26%); Results")
    add("T2 calibration", f"{t}_max_q05", int(r.null_max_n_q05), f"{t}: most genera at q < 0.05 in one shuffle", "results/denovo_scan/scan_summary.tsv: null_max_n_q05", "Table 2; Fig. 3a")
    add("T2 calibration", f"{t}_null_sd_W", float(r.null_W_sd), f"{t}: null SD of W", "results/denovo_scan/scan_summary.tsv: null_W_sd", "Table 2; Results (1.27 to 2.37)")
add("T2 calibration", "null_sd_W_range", f"{ss[ss.table.isin(T2)].null_W_sd.min():.2f} to {ss[ss.table.isin(T2)].null_W_sd.max():.2f}", "null SD of W across the five tables", "results/denovo_scan/scan_summary.tsv", "Results")
# per-shuffle false positive rates, exactly as Figure 3a (103): same files, same genus filter
tabs = [("nelore_genus", "scan_null", 2001), ("charolais_mrna_genus", "scan_null", 2001), ("alberta_genus", "scan_null", 2001), ("difford_bact", "methane_dryrun", 1001), ("wallace_prokaryotes", "methane_dryrun", 1001)]
for t, d, s0 in tabs:
    fs = [P("results", d, f"{t}_seed{s}_ancombc2.csv") for s in range(s0, s0 + 200)]; fs = [f for f in fs if os.path.exists(f)][:NSEED]
    assert len(fs) == NSEED, (t, len(fs))
    fpr, nq = [], []
    for f in fs:
        x = pd.read_csv(f); nm = x["taxon_name"] if "taxon_name" in x else x["taxon"]
        pc_ = [c for c in x.columns if re.match(r"^p_(RFI|CH4)", c)]; assert len(pc_) == 1; lc = pc_[0].replace("p_", "lfc_", 1); sc = pc_[0].replace("p_", "se_", 1)
        keep = ~nm.astype(str).str.startswith("TGT_") & np.isfinite(x[lc] / x[sc]); p = x.loc[keep, pc_[0]].to_numpy()
        fpr.append(np.mean(p < 0.05)); o = np.argsort(p); m = len(p); q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; nq.append(int(np.sum(q < 0.05)))
    fpr = np.array(fpr); nq = np.array(nq); src = f"results/{d}/{t}_seed{s0}..{s0 + 199}_ancombc2.csv (as Fig. 3a)"
    add("T2 calibration", f"{t}_fpr_median_pct", round(100 * float(np.median(fpr)), 1), f"{t}: median share of genera at nominal p < 0.05 across shuffles (%)", src, "Results (3.7%, 26%, 25%)")
    add("T2 calibration", f"{t}_fpr_max_pct", round(100 * float(fpr.max()), 1), f"{t}: largest share in one shuffle (%)", src, "Results (35% in Nelore)")
    add("T2 calibration", f"{t}_q05_in_fpr_max_shuffle", int(nq[fpr.argmax()]), f"{t}: genera at q < 0.05 in the shuffle with the largest share", src, "Results (679 in Nelore)")
    if not TEST: assert int(nq.max()) == int(one(ss, table=t).null_max_n_q05), (t, nq.max())
# targeted RFI features: median null SD per table; the largest methane null SD
e = rd("results/replication_rfi3_calibrated/estimates_calibrated_all.tsv"); m = e.groupby("table").sd_null.median()
add("T2 calibration", "rfi_target_null_sd_median_range", f"{m.min():.2f} to {m.max():.2f}", "median null SD of W of the targeted RFI features, per table", "results/replication_rfi3_calibrated/estimates_calibrated_all.tsv: sd_null", "Results (0.75 to 1.16)")
mp = rd("results/methane_pooled/estimates_calibrated.tsv"); s1 = one(mp, id="S1", cohort="wallace")
add("T2 calibration", "sharpea_wallace_W_raw", round(float(s1.W_raw), 1), "Sharpea, Wallace: raw ANCOM-BC2 W", "results/methane_pooled/estimates_calibrated.tsv: W_raw", "Results (-9.4)")
add("T2 calibration", "sharpea_wallace_null_sd", round(float(s1.sd_null), 2), "Sharpea, Wallace: null SD of W", "results/methane_pooled/estimates_calibrated.tsv: sd_null", "Results (5.3); Fig. 3b (5.26)")

# ---------------------------------------------------------------- depth paragraph
w = rd("results/denovo_scan/wallace_prokaryotes_scan.tsv"); wd = rd("results/denovo_scan/wallace_prokaryotes_depth_scan.tsv")
add("D depth", "wallace_median_absW", f"{w.W.abs().median():.2f} to {wd.W.abs().median():.2f}", "Wallace, median |W| across genera, unadjusted to depth-adjusted (raw ANCOM-BC2 W)", "results/denovo_scan/wallace_prokaryotes{,_depth}_scan.tsv: W", "Results (2.74 to 1.20)")
fr = rd("results/wallace_farms/farm_consistency.tsv"); frd = rd("results/wallace_farms_depth/farm_consistency.tsv")
add("D depth", "farm_robust_removed", f"{int(fr.farm_robust.sum())} to {int(frd.farm_robust.sum())}", "Wallace farm-robust features, unadjusted to depth-adjusted", "results/wallace_farms{,_depth}/farm_consistency.tsv: farm_robust", "Results (removed 4 of the 17)")
pl = open(P("plan", "ANALYSIS_PLAN.md"), encoding="utf-8").read()
need = "pooled +0.155 (p 7e-7); FI1 +0.18, IT1 +0.33, IT2 +0.35, IT3 +0.09, SE1 +0.03, UK1 -0.03, UK2 -0.04. Difford, within herd, +0.052 (p 0.16). Alberta, RFI within diet, -0.072 (p 0.09)"
assert need in pl, "Amendment 41 finding (depth correlations) not found in the plan"
for k, v, dsc in [("wallace_pooled", 0.155, "Wallace, methane within farm, pooled"), ("wallace_IT1", 0.33, "Wallace farm IT1"), ("wallace_IT2", 0.35, "Wallace farm IT2"), ("wallace_IT3", 0.09, "Wallace farm IT3"),
                  ("difford", 0.052, "Difford, within herd"), ("alberta_rfi", -0.072, "Alberta, RFI within diet")]:
    add("D depth", f"depth_rho_{k}", v, f"Spearman rho of log10 library size with the trait: {dsc}", "plan/ANALYSIS_PLAN.md, Amendment 41, item 1 (Finding)", "Results (0.16; two of the three Italian farms 0.33 and 0.35; 0.05; -0.07)")

# ---------------------------------------------------------------- estimator dependence: Eubacterium (N38) and Sharpea (S1), Wallace
a = rd("results/methane/wallace_prokaryotes_targets.tsv"); ad = rd("results/methane/wallace_prokaryotes_depth_targets.tsv")
c = rd("results/clr/wallace_prokaryotes_clr_targets.tsv"); cd = rd("results/clr/wallace_prokaryotes_depth_clr_targets.tsv"); dc = rd("results/clr/difford_bact_clr_targets.tsv")
for key, val, dsc, src in [
    ("eubacterium_ancombc2_W_raw", one(a, id="N38").W, "Eubacterium, Wallace: raw ANCOM-BC2 W", "results/methane/wallace_prokaryotes_targets.tsv: W"),
    ("eubacterium_clr_W_raw", one(c, id="N38").W, "Eubacterium, Wallace: raw CLR W", "results/clr/wallace_prokaryotes_clr_targets.tsv: W"),
    ("eubacterium_clr_W_depth", one(cd, id="N38").W, "Eubacterium, Wallace: CLR W, depth-adjusted", "results/clr/wallace_prokaryotes_depth_clr_targets.tsv: W"),
    ("sharpea_ancombc2_W_raw", one(a, id="S1").W, "Sharpea, Wallace: raw ANCOM-BC2 W", "results/methane/wallace_prokaryotes_targets.tsv: W"),
    ("sharpea_ancombc2_W_depth", one(ad, id="S1").W, "Sharpea, Wallace: raw ANCOM-BC2 W, depth-adjusted", "results/methane/wallace_prokaryotes_depth_targets.tsv: W"),
    ("sharpea_clr_W", one(c, id="S1").W, "Sharpea, Wallace: CLR W", "results/clr/wallace_prokaryotes_clr_targets.tsv: W"),
    ("sharpea_clr_W_depth", one(cd, id="S1").W, "Sharpea, Wallace: CLR W, depth-adjusted", "results/clr/wallace_prokaryotes_depth_clr_targets.tsv: W"),
    ("sharpea_difford_clr_W", one(dc, id="S1").W, "Sharpea, Difford: CLR W", "results/clr/difford_bact_clr_targets.tsv: W"),
    ("sharpea_difford_ancombc2_W_cal", one(mp, id="S1", cohort="difford").W, "Sharpea, Difford: calibrated ANCOM-BC2 W", "results/methane_pooled/estimates_calibrated.tsv: W")]:
    add("E estimators", key, float(val), dsc, src, "Results")
need = "Calibrated CLR, Wallace: N38 Eubacterium W -3.80"
assert need in pl, "Amendment 42 calibrated CLR result for Eubacterium not found in the plan"
add("E estimators", "eubacterium_clr_W_calibrated", -3.80, "Eubacterium, Wallace: calibrated CLR W (second estimator)", "plan/ANALYSIS_PLAN.md, Amendment 42 result", "Results")
dm = rd("results/methane_pooled_ch4_dmi/estimates_calibrated.tsv"); r = one(dm, id="M1", cohort="wallace")
add("E estimators", "clade_ratio_ch4dmi_wallace_W", round(float(r.W), 2), "clade ratio, Wallace, methane per kg DMI: calibrated W", "results/methane_pooled_ch4_dmi/estimates_calibrated.tsv: W", "Abstract; Results; Fig. 4")
add("E estimators", "clade_ratio_ch4dmi_wallace_r", round(float(r.r), 2), "clade ratio, Wallace, methane per kg DMI: r", "results/methane_pooled_ch4_dmi/estimates_calibrated.tsv: r", "Results")

# ---------------------------------------------------------------- LinDA and ALDEx2 null spread; effective n
for meth, d in [("LinDA", "results/late_linda/pooled/null_sd_by_table.tsv"), ("ALDEx2", "results/late_aldex2/pooled/null_sd_by_table.tsv")]:
    x = rd(d); add("L late methods", f"{meth}_null_sd_range", f"{x.sd_min.min():.2f} to {x.sd_max.max():.2f}", f"{meth}: null SD of W over targets and tables", f"{d}: sd_min, sd_max", "Results")
pe = rd("results/power/power_equivalence.tsv"); m = pe[pe.arm == "Methane"]
add("L late methods", "methane_ancombc2_neff_range", f"{int(m.n_eff.min()):,} to {int(m.n_eff.max()):,}", "methane, primary ANCOM-BC2 pooled intervals: effective n", "results/power/power_equivalence.tsv: n_eff", "Results (198 to 1,729)")

# ---------------------------------------------------------------- transfer positive controls and spike-in
tp = rd("results/late_positive_controls/transfer_positive_controls.tsv")
for _, x in tp.iterrows():
    k = re.sub(r"[^a-z0-9]+", "_", f"{x.comparison[:12]} {x.term}".lower()).strip("_")
    add("P controls", f"{k}_spearman", f"{x.spearman_W:.2f} (p {x.p_spearman:.3g}; {int(x.n_shared_genera)} genera; same-sign hits {int(x.same_sign_in_both)} vs {x.expected_same_sign:.1f} expected, p {x.p_sign if pd.notna(x.p_sign) else 'NA'}; pass {x['pass']})",
        f"{x.comparison}: {x.term}", "results/late_positive_controls/transfer_positive_controls.tsv", "Results; Discussion")
sp = rd("results/late_spikein/spikein_summary.tsv")
for _, x in sp.iterrows():
    add("P controls", f"spikein_{x.table}_{x.b_lfc_per_sd}", f"recovery calibrated {100 * x.power_calibrated:.1f}%; fpr raw {100 * x.fpr_raw:.1f}%, calibrated {100 * x.fpr_calibrated:.1f}%; est/planted {x.median_est_over_planted:.2f}; median calibrated r {x.median_r_calibrated:.3f}",
        f"spike-in {x.table}, {x.b_lfc_per_sd} natural-log fold change per SD", "results/late_spikein/spikein_summary.tsv", "Results; Discussion")
add("P controls", "spikein_est_over_planted_range", f"{sp.median_est_over_planted.min():.2f} to {sp.median_est_over_planted.max():.2f}", "spike-in: median estimated over planted fold change, range", "results/late_spikein/spikein_summary.tsv", "Results (0.95 to 1.02)")

head = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
M = pd.DataFrame(rows); M.insert(0, "repo_commit", head)
os.makedirs(os.path.dirname(out), exist_ok=True); M.to_csv(out + ".tmp", sep="\t", index=False); os.replace(out + ".tmp", out)
print(f"table numbers at commit {head}: {len(M)} rows in {M.section.nunique()} sections{' (TEST: ' + str(NSEED) + ' seeds per table)' if TEST else ''}")
print(M[["section", "key", "value"]].to_string(index=False))
