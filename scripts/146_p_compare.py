#!/usr/bin/env python3
"""146: parametric and empirical permutation p side by side for the ANCOM-BC2 target estimates (drafts 20 and 21; reporting only, no new fit or permutation).

For every primary target-cohort estimate the stored Amendment 37 values are W_raw, the permutation-null SD (sd_null) and the share of the k_null
within-design permutations with |W_null| >= |W_raw| (p_perm, b/B). From these:
  W_cal          = W_raw / max(1, sd_null), the rescaled W the labels use (it corrects the spread of the null, not its shape);
  p_param        = 2 * (1 - Phi(|W_cal|)), the normal-tail p of the rescaled W (draft 20);
  p_t            = 2 * pt(-|W_cal|, df), the Student-t p on the residual df (draft 21). This is the p the labels use for a target in one cohort,
                   including a leave-one-cohort-out analysis reduced to one cohort (scripts/lib_replication.R, pool_one, k = 1); for a target in two or
                   more cohorts the labels use a REML random-effects z test on Fisher z of W / sqrt(W^2 + df). Neither is an empirical permutation p;
  p_perm_emp     = (b + 1) / (B + 1), b = round(p_perm * B), B = k_null (200 for every estimate cited);
  agree_at_0.05, agree_t_at_0.05 = the normal-tail or Student-t p and the empirical p on the same side of 0.05.
Inputs (unchanged): results/late_a48/seven/seven_claims.tsv (48.1), results/replication_rfi3_calibrated/estimates_calibrated_all.tsv,
results/methane_pooled/estimates_calibrated.tsv, results/late_a46/effect_size/cohort_estimates.tsv (published directions).
Writes (refuses to overwrite): results/late_p_compare/{seven_claims_p.tsv, all_targets_p.tsv, summary.tsv}, results/manifest/draft20_numbers.tsv and
results/manifest/draft21_numbers.tsv.
Stops without writing if any value cited in draft 20 or draft 21 differs from what the files give. Usage: 146_p_compare.py [root]
"""
import os, sys, math, subprocess
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
OD = P("results", "late_p_compare")
OUT = [os.path.join(OD, f) for f in ("seven_claims_p.tsv", "all_targets_p.tsv", "summary.tsv")] + [P("results", "manifest", f"draft{d}_numbers.tsv") for d in (20, 21)]
ex = [f for f in OUT if os.path.exists(f)]
if ex: sys.exit(f"output exists, never overwritten: {ex}")
LABEL = "reporting only (drafts 20 and 21): computed from the stored Amendment 37 permutation results; no new fit or permutation"
rd = lambda *a: pd.read_csv(P(*a), sep="\t", dtype={"id": str})
p_norm = lambda w: math.erfc(abs(w) / math.sqrt(2))


def _betacf(a, b, x):
    # continued fraction for the regularized incomplete beta (modified Lentz)
    tiny, qab, qap, qam = 1e-300, a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny); h = d
    for m in range(1, 10001):
        m2 = 2 * m; aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > tiny else tiny); c = 1.0 + aa / c; c = c if abs(c) > tiny else tiny; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > tiny else tiny); c = 1.0 + aa / c; c = c if abs(c) > tiny else tiny; de = d * c; h *= de
        if abs(de - 1.0) < 1e-15: return h
    raise RuntimeError("betacf did not converge")
def betai(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    return bt * _betacf(a, b, x) / a if x < (a + 1) / (a + b + 2) else 1.0 - bt * _betacf(b, a, 1 - x) / b
def p_t2(w, df):
    """two-sided Student-t p, 2 * pt(-|w|, df), as lib_replication.R pool_one() for k = 1"""
    return betai(df / 2.0, 0.5, df / (df + w * w))


er = rd("results", "replication_rfi3_calibrated", "estimates_calibrated_all.tsv").assign(arm="RFI", reanalysis=False)
em = rd("results", "methane_pooled", "estimates_calibrated.tsv").assign(arm="Methane")
cols = ["arm", "id", "set", "cohort", "table", "rank", "n", "df", "W_raw", "sd_null", "k_null", "W", "p_perm", "reanalysis"]
E = pd.concat([er[cols], em[cols]], ignore_index=True)
assert not E.duplicated(["id", "cohort"]).any() and er.status.eq("ok").all() and em.status.eq("ok").all()
# W must be W_raw / max(1, sd_null), and p_perm a multiple of 1/k_null
assert (E.W - E.W_raw / E.sd_null.clip(lower=1)).abs().max() < 1e-9
assert ((E.p_perm * E.k_null - (E.p_perm * E.k_null).round()).abs() < 1e-6).all()
ce = rd("results", "late_a46", "effect_size", "cohort_estimates.tsv")
dirs = ce.groupby("id").published_dir.agg(lambda x: sorted(set(x.dropna())))   # empty for targets without a published direction (C1 to C6, T3, M3)
assert all(len(v) <= 1 for v in dirs), "a target has more than one published direction"
E["published_dir"] = E.id.map(lambda i: dirs[i][0] if i in dirs and len(dirs[i]) else float("nan"))
E["b"] = (E.p_perm * E.k_null).round().astype(int)
E["p_param_ancombc2"] = E.W.map(p_norm)
E["p_perm_emp_ancombc2"] = (E.b + 1) / (E.k_null + 1)
E["agree_at_0.05"] = (E.p_param_ancombc2 < 0.05) == (E.p_perm_emp_ancombc2 < 0.05)
E["p_t_ancombc2"] = [p_t2(w, d) for w, d in zip(E.W, E.df)]
E["agree_t_at_0.05"] = (E.p_t_ancombc2 < 0.05) == (E.p_perm_emp_ancombc2 < 0.05)
E["W_direction"] = [("no published direction" if pd.isna(d) else ("published" if math.copysign(1, w) == d else "against")) for w, d in zip(E.W, E.published_dir)]
E["in_61_tests"] = ~E.reanalysis.astype(bool)
E = E.rename(columns={"W": "W_ancombc2_cal", "p_perm": "p_perm_ancombc2_b_over_B (as stored)", "k_null": "B"})
E["analysis"] = LABEL

# ---- seven claims (48.1)
S = rd("results", "late_a48", "seven", "seven_claims.tsv")
K = S[["id", "target", "arm", "cohort", "published_dir", "W_ancombc2_cal", "p_perm_ancombc2"]].merge(
    E[["id", "cohort", "df", "B", "b", "W_ancombc2_cal", "p_param_ancombc2", "p_perm_emp_ancombc2", "agree_at_0.05", "p_t_ancombc2", "agree_t_at_0.05", "W_direction"]], on=["id", "cohort"], suffixes=("", "_e"), validate="1:1")
assert len(K) == 14 and (K.W_ancombc2_cal - K.W_ancombc2_cal_e).abs().max() < 1e-12 and (K.p_perm_ancombc2 - K.b / K.B).abs().max() < 1e-12
K = K.drop(columns=["W_ancombc2_cal_e"]).rename(columns={"p_perm_ancombc2": "p_perm_ancombc2_b_over_B (as stored)"}); K["analysis"] = LABEL

# ---- values cited in draft 20 (review chat, DRAFT20_NOTES_FOR_SERVER.md item 1): W, parametric p, empirical p to 2, 3 and 3 decimals
EXP = {("T2", "nelore"): (0.24, .808, .746), ("T2", "charolais_mrna"): (0.19, .846, .746), ("T4", "nelore"): (-1.36, .175, .194), ("T4", "charolais_mrna"): (-2.04, .042, .020),
       ("T5", "nelore"): (0.73, .467, .512), ("T5", "charolais_mrna"): (-0.14, .890, .876), ("T6", "nelore"): (-0.99, .324, .224), ("T6", "charolais_mrna"): (1.94, .052, .030),
       ("N15", "charolais_mrna"): (0.82, .411, .264), ("N15", "alberta_all"): (2.76, .006, .005), ("S1", "difford"): (-0.13, .900, .886), ("S1", "wallace"): (-1.78, .075, .075),
       ("M1", "wallace"): (0.82, .414, .408), ("M1", "difford"): (0.85, .397, .383)}
bad = []
for _, r in K.iterrows():
    w, pp, pe = EXP[(r.id, r.cohort)]
    got = (round(r.W_ancombc2_cal, 2), round(r.p_param_ancombc2, 3), round(r.p_perm_emp_ancombc2, 3))
    if got != (w, pp, pe): bad.append((r.id, r.cohort, got, (w, pp, pe)))
dis = K[~K["agree_at_0.05"]]
emp_pub = K[(K.p_perm_emp_ancombc2 < 0.05) & (K.W_direction == "published")]
both = emp_pub.groupby("id").size()
if bad: sys.exit(f"cited values differ from the files, nothing written: {bad}")
if not (K["agree_at_0.05"].sum() == 13 and len(dis) == 1 and dis.iloc[0][["id", "cohort", "W_direction"]].tolist() == ["T6", "charolais_mrna", "against"]):
    sys.exit(f"agreement count differs from draft 20, nothing written: {K['agree_at_0.05'].sum()} of 14; {dis[['id', 'cohort']].values.tolist()}")
if not (len(emp_pub) == 1 and emp_pub.iloc[0][["id", "cohort"]].tolist() == ["N15", "alberta_all"] and (both < 2).all()):
    sys.exit(f"published-direction empirical p < 0.05 differs from draft 20, nothing written: {emp_pub[['id', 'cohort']].values.tolist()}")

# ---- draft 21: Student-t reference (values cited in DRAFT21_NOTES_FOR_SERVER.md item 1)
EXPT = {("T4", "charolais_mrna"): (.052, .020), ("T6", "charolais_mrna"): (.063, .030), ("N15", "alberta_all"): (.006, .005)}
badt = [(i, c, round(K.loc[(K.id == i) & (K.cohort == c), "p_t_ancombc2"].iloc[0], 3), v) for (i, c), v in EXPT.items()
        if (round(K.loc[(K.id == i) & (K.cohort == c), "p_t_ancombc2"].iloc[0], 3), round(K.loc[(K.id == i) & (K.cohort == c), "p_perm_emp_ancombc2"].iloc[0], 3)) != v]
if badt: sys.exit(f"draft 21 Student-t values differ from the files, nothing written: {badt}")
disT = K[~K["agree_t_at_0.05"]]
if not (int(K["agree_t_at_0.05"].sum()) == 12 and sorted(map(tuple, disT[["id", "cohort"]].values)) == [("T4", "charolais_mrna"), ("T6", "charolais_mrna")] and (disT.W_direction == "against").all()):
    sys.exit(f"Student-t agreement differs from draft 21, nothing written: {int(K['agree_t_at_0.05'].sum())} of 14; {disT[['id', 'cohort']].values.tolist()}")

# ---- all primary target-cohort estimates (optional item; report only)
A = E[E.in_61_tests]; dA = A[~A["agree_at_0.05"]]; dT = A[~A["agree_t_at_0.05"]]
if not (len(dT) == 3 and sorted(map(tuple, dT[["id", "cohort"]].values)) == [("N37", "wallace"), ("T4", "charolais_mrna"), ("T6", "charolais_mrna")] and (dT.W_direction == "against").all()
        and round(dT.loc[dT.id == "N37", "p_t_ancombc2"].iloc[0], 3) == 0.041):
    sys.exit(f"Student-t disagreements over the 61 tests differ from draft 21, nothing written: {dT[['id', 'cohort']].values.tolist()}")
summ = pd.DataFrame([
    ("seven claims: cohort estimates", 14, ""), ("seven claims: agree at 0.05", int(K["agree_at_0.05"].sum()), ""),
    ("seven claims: disagree at 0.05", len(dis), "; ".join(f"{r.id} {r.cohort} (parametric {r.p_param_ancombc2:.3f}, empirical {r.p_perm_emp_ancombc2:.3f}, W {r.W_direction} the published direction)" for r in dis.itertuples())),
    ("seven claims: empirical p < 0.05 in the published direction", len(emp_pub), "; ".join(f"{r.id} {r.cohort} (empirical {r.p_perm_emp_ancombc2:.3f}, parametric {r.p_param_ancombc2:.3f})" for r in emp_pub.itertuples())),
    ("seven claims: empirical p < 0.05 in the published direction in both cohorts", int((both >= 2).sum()), ""),
    ("all primary target-cohort tests (61; the M2 Wallace reanalysis excluded)", len(A), ""),
    ("all primary tests: disagree at 0.05", len(dA), "; ".join(f"{r.id} {r.cohort} (parametric {r.p_param_ancombc2:.3f}, empirical {r.p_perm_emp_ancombc2:.3f}, W {r.W_direction})" for r in dA.itertuples())),
    ("all primary tests: disagreements in the published direction", int((dA.W_direction == "published").sum()), ""),
    ("all primary tests: parametric p < 0.05", int((A.p_param_ancombc2 < 0.05).sum()), ""),
    ("all primary tests: empirical p < 0.05", int((A.p_perm_emp_ancombc2 < 0.05).sum()), ""),
    ("all primary tests: empirical p < 0.05 in the published direction", int(((A.p_perm_emp_ancombc2 < 0.05) & (A.W_direction == "published")).sum()),
     "; ".join(f"{r.id} {r.cohort}" for r in A[(A.p_perm_emp_ancombc2 < 0.05) & (A.W_direction == "published")].itertuples())),
    ("all primary tests: without a published direction", int((A.W_direction == "no published direction").sum()), ""),
    ("seven claims: Student-t p and empirical p agree at 0.05 (draft 21)", int(K["agree_t_at_0.05"].sum()), "the labels use the Student-t p on the residual df for a target in one cohort"),
    ("seven claims: Student-t p and empirical p disagree at 0.05", len(disT), "; ".join(f"{r.id} {r.cohort} (Student-t {r.p_t_ancombc2:.3f}, empirical {r.p_perm_emp_ancombc2:.3f}, W {r.W_direction} the published direction)" for r in disT.itertuples())),
    ("all primary tests: Student-t p and empirical p disagree at 0.05", len(dT), "; ".join(f"{r.id} {r.cohort} (Student-t {r.p_t_ancombc2:.3f}, empirical {r.p_perm_emp_ancombc2:.3f}, W {r.W_direction})" for r in dT.itertuples())),
    ("all primary tests: Student-t disagreements in the published direction", int((dT.W_direction == "published").sum()), ""),
    ("all primary tests: Student-t p < 0.05", int((A.p_t_ancombc2 < 0.05).sum()), ""),
    ("B per estimate", "; ".join(f"{k}: {v}" for k, v in E.B.value_counts().sort_index().items()), "number of within-design permutations stored per estimate"),
], columns=["item", "value", "detail"]); summ["analysis"] = LABEL

head = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
t6 = K[(K.id == "T6") & (K.cohort == "charolais_mrna")].iloc[0]; n15 = K[(K.id == "N15") & (K.cohort == "alberta_all")].iloc[0]
M = pd.DataFrame([
    ("Results, ANCOM-BC2 paragraph", "parametric and empirical p agree at 0.05", "13 of 14", f"{int(K['agree_at_0.05'].sum())} of 14", "results/late_p_compare/seven_claims_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "M. smithii (T6) Charolais, empirical p", "0.030", f"{t6.p_perm_emp_ancombc2:.3f}", "results/late_p_compare/seven_claims_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "M. smithii (T6) Charolais, parametric p", "0.052", f"{t6.p_param_ancombc2:.3f}", "results/late_p_compare/seven_claims_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "T6 Charolais W relative to the published direction", "against", t6.W_direction, "results/late_p_compare/seven_claims_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "claims with empirical p < 0.05 in the published direction in both cohorts", "0", str(int((both >= 2).sum())), "results/late_p_compare/seven_claims_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "Selenomonas (N15) Alberta, empirical p", "0.005", f"{n15.p_perm_emp_ancombc2:.3f}", "results/late_p_compare/seven_claims_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "Selenomonas (N15) Alberta, parametric p", "0.006", f"{n15.p_param_ancombc2:.3f}", "results/late_p_compare/seven_claims_p.tsv"),
    ("Additional file 2, sheet 'A48 seven claims, per cohort'", "14 W, parametric and empirical p (review chat table)", "as tabulated", "all 14 match", "results/late_p_compare/seven_claims_p.tsv"),
], columns=["location", "quantity", "printed", "file_value", "source"])
assert (M.printed.where(M.printed != "as tabulated", "all 14 match") == M.file_value).all()
M["commit"] = head; M["analysis"] = LABEL

t4 = K[(K.id == "T4") & (K.cohort == "charolais_mrna")].iloc[0]; n37 = A[(A.id == "N37") & (A.cohort == "wallace")].iloc[0]
SRC = "results/late_p_compare/seven_claims_p.tsv"
M21 = pd.DataFrame([
    ("Results, ANCOM-BC2 paragraph", "Student-t p and empirical p agree at 0.05, seven claims", "12 of 14", f"{int(K['agree_t_at_0.05'].sum())} of 14", SRC),
    ("Results, ANCOM-BC2 paragraph", "T4 Charolais, Student-t p", "0.052", f"{t4.p_t_ancombc2:.3f}", SRC),
    ("Results, ANCOM-BC2 paragraph", "T4 Charolais, empirical p", "0.020", f"{t4.p_perm_emp_ancombc2:.3f}", SRC),
    ("Results, ANCOM-BC2 paragraph", "T6 Charolais, Student-t p", "0.063", f"{t6.p_t_ancombc2:.3f}", SRC),
    ("Results, ANCOM-BC2 paragraph", "T6 Charolais, empirical p", "0.030", f"{t6.p_perm_emp_ancombc2:.3f}", SRC),
    ("Results, ANCOM-BC2 paragraph", "T4 and T6 Charolais W relative to the published direction", "against", "against" if (disT.W_direction == "against").all() else "not all against", SRC),
    ("Results, ANCOM-BC2 paragraph", "N15 Alberta, Student-t p", "0.006", f"{n15.p_t_ancombc2:.3f}", SRC),
    ("Results, ANCOM-BC2 paragraph", "61 primary tests, Student-t vs empirical disagreements at 0.05", "3", str(len(dT)), "results/late_p_compare/all_targets_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "N37 Wallace, Student-t p", "0.041", f"{n37.p_t_ancombc2:.3f}", "results/late_p_compare/all_targets_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "N37 Wallace, empirical p", "0.055", f"{n37.p_perm_emp_ancombc2:.3f}", "results/late_p_compare/all_targets_p.tsv"),
    ("Results, ANCOM-BC2 paragraph", "Student-t disagreements in the published direction", "0", str(int((dT.W_direction == "published").sum())), "results/late_p_compare/all_targets_p.tsv"),
], columns=["location", "quantity", "printed", "file_value", "source"])
assert (M21.printed == M21.file_value).all(), M21[M21.printed != M21.file_value]
M21["commit"] = head; M21["analysis"] = LABEL

os.makedirs(OD, exist_ok=True)
for f, x in zip(OUT, (K, E, summ, M, M21)):
    tmp = f + ".tmp"; x.to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)
print(f"seven claims: {int(K['agree_at_0.05'].sum())} of 14 agree at 0.05; disagreement: " + "; ".join(f"{r.id} {r.cohort} {r.p_param_ancombc2:.3f} vs {r.p_perm_emp_ancombc2:.3f} ({r.W_direction})" for r in dis.itertuples()))
print(f"all {len(A)} primary tests: {len(dA)} disagree at 0.05 (" + "; ".join(f"{r.id} {r.cohort} {r.p_param_ancombc2:.3f} vs {r.p_perm_emp_ancombc2:.3f}, {r.W_direction}" for r in dA.itertuples()) + f"); {int((dA.W_direction == 'published').sum())} in the published direction")
print(f"Student-t reference: seven claims {int(K['agree_t_at_0.05'].sum())} of 14 agree; 61 tests: {len(dT)} disagree (" + "; ".join(f"{r.id} {r.cohort} {r.p_t_ancombc2:.3f} vs {r.p_perm_emp_ancombc2:.3f}, {r.W_direction}" for r in dT.itertuples()) + ")")
print(f"draft 20 and 21 numbers at commit {head}: {len(M)} + {len(M21)} rows; every printed value matches; wrote", ", ".join(os.path.relpath(f, ROOT) for f in OUT))
