#!/usr/bin/env python3
"""127: the numbers draft 17 adds to the manuscript that neither 119, 124 nor 126 records. Each is computed from a committed result file (no fit,
no hand-typed value). (a) Charolais end-to-end coverage per target over the nine rho values (47.1 check 2); (b) the single-cohort intervals whose bound
in the published direction lies just under r = 0.15; (c) the claims matched to a published direction in rumen samples (stratum "directed rumen", both arms),
split by the number of cohorts, with their r = 0.30 exclusions; (d) the two most extreme Nelore null shuffles, read as in Figure 3a (103) and 126.
Run: python3 scripts/127_manifest_draft17.py   (pandas). Writes results/manifest/draft17_numbers.tsv (never overwrites).
"""
import os, re, subprocess, sys
import numpy as np, pandas as pd

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort"); P = lambda *a: os.path.join(ROOT, *a)
out = P("results", "manifest", "draft17_numbers.tsv")
if os.path.exists(out): sys.exit(f"output exists, never overwritten: {out}")
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
rows = []
def add(section, key, value, desc, src, used=""):
    if isinstance(value, (np.integer,)): value = int(value)
    elif isinstance(value, (float, np.floating)): value = round(float(value), 4)
    rows.append(dict(section=section, key=key, value=value, description=desc, source=src, used_in=used))

# (a) Charolais end-to-end coverage per target
src = "results/late_a47/intervals/check2_cohort.tsv"; c = rd(src); c = c[c.cohort == "charolais_mrna"]
for i in ["N02", "N20", "N22"]:
    x = c[c.id == i]; assert len(x) == 9 and x.rho.nunique() == 9, (i, len(x))
    add("A coverage", f"charolais_{i}_coverage_range", f"{x.coverage.min():.3f} to {x.coverage.max():.3f}",
        f"{i}, Charolais: end-to-end coverage of the 95% interval over rho in 0, +/-0.10, +/-0.15, +/-0.20, +/-0.30 (reps {int(x.reps.min())} to {int(x.reps.max())})",
        f"{src} (check 2, coverage)", "Results; Limitations")
add("A coverage", "charolais_coverage_min_all_targets", f"{c.coverage.min():.3f}", "Charolais: lowest end-to-end coverage over all its targets and rho", f"{src} (check 2, coverage)", "context")

# (b) bounds just under r = 0.15
src = "results/late_a47/intervals/pooled_intervals.tsv"; p = rd(src); p = p[p.version == "primary"]
for i in ["N11", "N20", "N22"]:
    r = p[p.id == i]; assert len(r) == 1; r = r.iloc[0]; assert r.k == 1 and bool(r.excl_15_mkh)
    b = -r.lo_mkh if r.published_dir < 0 else r.hi_mkh
    add("B bounds", f"{i}_bound_published_dir", f"{b:.3f}", f"{i} {r.target} ({r.cohorts}): interval bound in the published direction, oriented (excludes r = 0.15 by a margin under 0.01)",
        f"{src}: lo_mkh/hi_mkh, published_dir", "Results")

# (c) claims matched to a published direction in rumen samples
m = p[p.stratum == "directed rumen"]; ex = m[m.excl_30_mkh.astype(bool)]
add("C matched", "matched_claims", len(m), "claims matched to a published direction in rumen samples (both arms)", f"{src}: stratum", "Results")
add("C matched", "matched_claims_k_ge2", int((m.k >= 2).sum()), "matched claims tested in two or more cohorts", f"{src}: stratum, k", "Results")
add("C matched", "matched_claims_k1", int((m.k == 1).sum()), "matched claims tested in one cohort", f"{src}: stratum, k", "Results")
add("C matched", "matched_excl30", f"{len(ex)} ({(ex.arm == 'RFI').sum()} RFI, {(ex.arm == 'Methane').sum()} methane): {', '.join(ex.id)}", "matched claims whose interval excludes r = 0.30", f"{src}: excl_30_mkh", "Results")
add("C matched", "matched_excl30_k_ge2", int((ex.k >= 2).sum()), "of those, tested in two or more cohorts", f"{src}: excl_30_mkh, k", "Results")
assert len(m) == (m.k >= 2).sum() + (m.k == 1).sum()

# (d) Nelore null shuffles
fs = [P("results", "scan_null", f"nelore_genus_seed{s}_ancombc2.csv") for s in range(2001, 2201)]; assert all(map(os.path.exists, fs))
r = []
for s, f in zip(range(2001, 2201), fs):
    x = pd.read_csv(f); nm = x["taxon_name"] if "taxon_name" in x else x["taxon"]
    pc = [k for k in x.columns if re.match(r"^p_(RFI|CH4)", k)]; assert len(pc) == 1; pc = pc[0]
    keep = ~nm.astype(str).str.startswith("TGT_") & np.isfinite(x[pc.replace("p_", "lfc_", 1)] / x[pc.replace("p_", "se_", 1)])
    pv = np.sort(x.loc[keep, pc].to_numpy()); n = len(pv); q = np.minimum.accumulate((pv * n / np.arange(1, n + 1))[::-1])[::-1]
    r.append((s, 100 * np.mean(pv < 0.05), int((q < 0.05).sum())))
d = pd.DataFrame(r, columns=["seed", "fpr_pct", "n_q05"]); a = d.loc[d.fpr_pct.idxmax()]; b = d.loc[d.n_q05.idxmax()]
src = "results/scan_null/nelore_genus_seed2001..2200_ancombc2.csv (as Fig. 3a)"
add("D Nelore null", "nelore_max_fpr_shuffle", f"{a.fpr_pct:.1f}% at p < 0.05, {int(a.n_q05)} at q < 0.05 (seed {int(a.seed)})", "Nelore: the shuffle with the largest share of genera at nominal p < 0.05", src, "Results (up to 35%)")
add("D Nelore null", "nelore_max_q05_shuffle", f"{int(b.n_q05)} at q < 0.05, {b.fpr_pct:.1f}% at p < 0.05 (seed {int(b.seed)})", "Nelore: the shuffle with the most genera at q < 0.05", src, "Results (up to 679); Table 2")

head = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
M = pd.DataFrame(rows); M.insert(0, "repo_commit", head)
os.makedirs(os.path.dirname(out), exist_ok=True); M.to_csv(out + ".tmp", sep="\t", index=False); os.replace(out + ".tmp", out)
print(f"draft 17 numbers at commit {head}: {len(M)} rows"); print(M[["section", "key", "value"]].to_string(index=False))
