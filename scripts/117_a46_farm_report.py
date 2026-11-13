#!/usr/bin/env python3
"""117: Amendment 46.5.4 (late analysis, run after all primary results were seen; reporting only, no new fit).
Wallace leave-one-farm-out (Amendment 40) reported with what stability against the full-data threshold does not show: per-farm estimates with
95% intervals, the number of farms fitted, sign agreement and I^2, for the genera fixed by Amendment 40 item 2.
Per-farm estimates are the uncalibrated ANCOM-BC2 within-farm fits of script 94 (results/wallace_farms/per_farm.tsv; trait CH4 g/d standardised
within farm); interval = lfc +/- 1.96 SE (model-based, not permutation-calibrated). I^2 is the REML value from script 94 (farms >= 3).
Run: python3 scripts/117_a46_farm_report.py. Writes results/late_a46/wallace_farms/ (never overwrites).
"""
import os
import pandas as pd

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort"); P = lambda *a: os.path.join(ROOT, *a)
LABEL = "late analysis, run after all primary results were seen (Amendment 46.5.4)"
od = P("results", "late_a46", "wallace_farms"); os.makedirs(od, exist_ok=True)
outs = [os.path.join(od, f) for f in ("per_farm_intervals.tsv", "summary.tsv")]
ex = [f for f in outs if os.path.exists(f)]; assert not ex, f"output exists, Amendment 46 never overwrites: {ex}"
fc = pd.read_csv(P("results", "wallace_farms", "farm_consistency.tsv"), sep="\t")
pf = pd.read_csv(P("results", "wallace_farms", "per_farm.tsv"), sep="\t")
full = pd.read_csv(P("results", "wallace_farms", "full.tsv"), sep="\t").set_index("genus")
farms = sorted(pf.farm.unique())
x = pf[pf.genus.isin(fc.genus)].copy()
x["lo_95"] = x.lfc - 1.96 * x.se; x["hi_95"] = x.lfc + 1.96 * x.se
x["sign_as_full"] = [ (w > 0) == (full.loc[g, "W"] > 0) for g, w in zip(x.genus, x.W)]
x["interval_excludes_0"] = (x.lo_95 > 0) | (x.hi_95 < 0)
x["analysis"] = LABEL
S = fc[["genus", "W_full", "null_q95", "W_lofo_min", "W_lofo_max", "lofo_sign_kept", "farm_robust", "farms_fit", "farms_same_sign", "I2", "weakest_without"]].copy()
S = S.rename(columns={"farm_robust": "lofo_stable_vs_full_data_threshold"})
g = x.groupby("genus")
S["farms_interval_excludes_0_same_sign"] = S.genus.map(g.apply(lambda d: int((d.interval_excludes_0 & d.sign_as_full).sum())))
S["farms_interval_excludes_0_opposite_sign"] = S.genus.map(g.apply(lambda d: int((d.interval_excludes_0 & ~d.sign_as_full).sum())))
S["farms_not_fitted"] = S.genus.map(g.apply(lambda d: ";".join(f for f in farms if f not in set(d.farm)))).fillna(";".join(farms))
S["note"] = ("leave-one-farm-out compares each refit with the full-data null 95th percentile of |W| (stability), not with its own null; "
             "per-farm fits are uncalibrated and model-based; a genus below 10% prevalence in a farm is not fitted there")
S["analysis"] = LABEL
x[["genus", "farm", "lfc", "se", "lo_95", "hi_95", "W", "sign_as_full", "interval_excludes_0", "analysis"]].sort_values(["genus", "farm"]).to_csv(outs[0], sep="\t", index=False)
S.to_csv(outs[1], sep="\t", index=False)
n = len(S)
print(LABEL)
print(f"genera {n}; sign kept in every leave-one-farm-out refit: {int(S.lofo_sign_kept.sum())}; also above the full-data null 95th percentile: {int(S.lofo_stable_vs_full_data_threshold.sum())}")
print(f"farms fitted per genus: median {S.farms_fit.median():.0f} (range {S.farms_fit.min()} to {S.farms_fit.max()}); same sign as the full fit in every fitted farm: {int((S.farms_same_sign == S.farms_fit).sum())}; in all 7 farms: {int((S.farms_same_sign == 7).sum())}")
print(f"I2 median {S.I2.median():.0f}% (range {S.I2.min():.0f} to {S.I2.max():.0f}; {int(S.I2.isna().sum())} genera with fewer than 3 farms)")
print(f"farms whose own interval excludes 0 in the full-fit direction, median per genus {S.farms_interval_excludes_0_same_sign.median():.0f}; genera with any farm significant in the opposite direction: {int((S.farms_interval_excludes_0_opposite_sign > 0).sum())}")
print(S[["genus", "W_full", "farms_fit", "farms_same_sign", "farms_interval_excludes_0_same_sign", "farms_interval_excludes_0_opposite_sign", "I2", "lofo_stable_vs_full_data_threshold"]].to_string(index=False))
