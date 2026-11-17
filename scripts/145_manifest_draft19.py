#!/usr/bin/env python3
"""145: the numbers draft 19 adds or rewords (DRAFT19_NOTES_FOR_SERVER.md), each computed from a committed result file (no fit, no hand-typed
value) and checked against the value draft 19 prints; the script stops on any mismatch. Late analysis, run after all primary results were seen
(Amendment 48, wording note for draft 19; reporting only).
  (a) ANCOM-BC2 null: unfiltered and filtered medians on the same statistic, and the Table 2 means (ancombc2/null_summary.tsv);
  (b) Selenomonas (N15) in Alberta, CLR (seven/seven_claims.tsv);
  (c) robustness of the six Charolais r = 0.30 exclusions over the 12 transformation cells (clr_grid/cohort_cells.tsv) and the 30
      leave-one-animal-out refits (charolais/loao.tsv), with the exact-genus (MC1) re-tabulation (mapping/mc1_pooled.tsv), and the same at r = 0.15;
  (d) single-cohort share of the 14 benchmark exclusions (late_a47 pooled_intervals.tsv);
  (e) dairy comparison, descriptive (dairy/agreement.tsv);
  (f) independent target-cohort tests and possible overlap (late_a46/ledger/overlap_audit.tsv).
Run: python3 scripts/145_manifest_draft19.py [root]   (pandas). Writes results/manifest/draft19_numbers.tsv (never overwrites).
"""
import os, sys, subprocess
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
out = P("results", "manifest", "draft19_numbers.tsv")
if os.path.exists(out): sys.exit(f"output exists, never overwritten: {out}")
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
rows, bad = [], []
def add(sec, key, value, printed, desc, src):
    rows.append(dict(section=sec, key=key, value=value, printed_in_draft19=printed, description=desc, source=src))
    if printed is not None and str(value) != str(printed): bad.append((key, value, printed))
pct = lambda x: f"{100 * x:.0f}"

# (a) ANCOM-BC2 null summaries
src = "results/late_a48/ancombc2/null_summary.tsv"; N = rd(src); assert len(N) == 5
for _, r in N.iterrows():
    add("A ANCOM-BC2 null", f"{r.table}_median_unfiltered", f"{r.median_share_p05:.3f}", None, f"{r.table}: median share of genera with p < 0.05 over 200 shuffles", src)
    add("A ANCOM-BC2 null", f"{r.table}_median_filtered", f"{r.median_share_p05_passed:.3f}", None, f"{r.table}: median share with p < 0.05 and passed_ss", src)
    add("A ANCOM-BC2 null", f"{r.table}_mean_unfiltered", f"{r.mean_share_p05:.3f}", None, f"{r.table}: mean share with p < 0.05 (Table 2)", src)
add("A ANCOM-BC2 null", "median_unfiltered_range", f"{pct(N.median_share_p05.min())} to {pct(N.median_share_p05.max())}%", "4 to 26%", "unfiltered median, range over the five tables", src)
add("A ANCOM-BC2 null", "median_filtered_range", f"{pct(N.median_share_p05_passed.min())} to {pct(N.median_share_p05_passed.max())}%", "2 to 9%", "filtered (passed_ss) median, range over the five tables", src)
add("A ANCOM-BC2 null", "mean_unfiltered_range", f"{pct(N.mean_share_p05.min())} to {pct(N.mean_share_p05.max())}%", "6 to 26%", "unfiltered mean (Table 2), range over the five tables", src)

# (b) Selenomonas, Alberta
src = "results/late_a48/seven/seven_claims.tsv"; S = rd(src); x = S[(S.id == "N15") & (S.cohort == "alberta_all")].iloc[0]
add("B Selenomonas", "N15_alberta_clr", f"{x.r_clr:.3f} ({x.lo_clr:.3f} to {x.hi_clr:.3f})", "0.033 (-0.050 to 0.116)", f"Selenomonas (N15), Alberta (n {x.n}): CLR partial r (95% interval)", src)

# (c) Charolais robustness
C = rd("results/late_a48/clr_grid/cohort_cells.tsv"); L = rd("results/late_a48/charolais/loao.tsv"); M = rd("results/late_a48/mapping/mc1_pooled.tsv")
src = "results/late_a48/clr_grid/cohort_cells.tsv; results/late_a48/charolais/loao.tsv"
PD = {"N02": -1, "N04": 1, "N11": -1, "N13": 1, "N20": -1, "N22": -1}
EXP15 = {"N04": ("10/12", "30/30"), "N11": ("10/12", "11/30"), "N13": ("12/12", "30/30"), "N20": ("12/12", "14/30"), "N22": ("8/12", "11/30")}
all30 = True
for i, d in PD.items():
    c = C[(C.id == i) & (C.cohort == "charolais_mrna")]; l = L[L.id == i]; assert len(c) == 12 and len(l) == 30, i
    bc = c.hi if d > 0 else -c.lo; bl = l.hi if d > 0 else -l.lo
    c30, l30 = f"{int((bc < 0.30).sum())}/12", f"{int((bl < 0.30).sum())}/30"; all30 &= (c30 == "12/12" and l30 == "30/30")
    add("C Charolais r = 0.30", f"{i}_cells_excl30", c30, "12/12", f"{i}: transformation cells in which the Charolais interval excludes r = 0.30 in the published direction", src)
    add("C Charolais r = 0.30", f"{i}_refits_excl30", l30, "30/30", f"{i}: leave-one-animal-out refits in which it excludes r = 0.30", src)
    if i in EXP15:
        add("C Charolais r = 0.15", f"{i}_cells_excl15", f"{int((bc < 0.15).sum())}/12", EXP15[i][0], f"{i}: cells excluding r = 0.15", src)
        add("C Charolais r = 0.15", f"{i}_refits_excl15", f"{int((bl < 0.15).sum())}/30", EXP15[i][1], f"{i}: refits excluding r = 0.15", src)
    m = M[(M.arm == "RFI") & (M.id == i)].iloc[0]
    mc1 = "no exact-genus estimate" if m.k == 0 else ("excluded" if str(m.excl_30) == "True" else "not excluded")
    add("C Charolais r = 0.30", f"{i}_mc1_r030", mc1, "no exact-genus estimate" if i in ("N02", "N11", "N13") else "excluded",
        f"{i}: r = 0.30 with exact-genus (MC1) estimates only", "results/late_a48/mapping/mc1_pooled.tsv")
add("C Charolais r = 0.30", "all_six_hold", str(all30), "True", "all six Charolais r = 0.30 exclusions hold in 12/12 cells and 30/30 refits", src)
held15 = [i for i in EXP15 if (lambda c, l: (c.hi if PD[i] > 0 else -c.lo).lt(0.15).all() and (l.hi if PD[i] > 0 else -l.lo).lt(0.15).all())(
          C[(C.id == i) & (C.cohort == "charolais_mrna")], L[L.id == i])]
add("C Charolais r = 0.15", "held_every_check", ";".join(held15), "N13", "Charolais r = 0.15 exclusions that hold in every cell and refit (Eubacterium)", src)

# (d) single-cohort share of the benchmark exclusions
src = "results/late_a47/intervals/pooled_intervals.tsv"; Pp = rd(src); Pp = Pp[Pp.version == "primary"]
ex = Pp[(Pp.arm.eq("RFI") & Pp.excl_30_mkh.astype(bool)) | (Pp.arm.eq("Methane") & Pp.excl_15_mkh.astype(bool))]
add("D exclusions", "benchmark_exclusions", len(ex), 14, "benchmark exclusions (RFI r = 0.30, methane r = 0.15)", src)
add("D exclusions", "single_cohort", int((ex.k == 1).sum()), 12, "of those, resting on one cohort", src)
add("D exclusions", "charolais_only", int(((ex.k == 1) & (ex.cohorts == "charolais_mrna")).sum()), 6, "of those, on the 30 Charolais animals", src)

# (e) dairy comparison
src = "results/late_a48/dairy/agreement.tsv"; A = rd(src); A = A[A.subset == "all shared genera"]
for v, pr in (("CH4", ("-0.37", "24 of 61")), ("CH4, depth-adjusted", ("-0.36", "26 of 61")), ("days in milk (positive control)", ("0.21", "37 of 61"))):
    a = A[A.variant == v].iloc[0]
    add("E dairy", f"rho|{v}", f"{a.spearman_rho:.2f}", pr[0], f"Spearman rho of CLR partial r between Wallace and Difford, {v}; descriptive", src)
    add("E dairy", f"same_sign|{v}", f"{int(a.same_sign)} of {int(a.genera)}", pr[1], f"genera with the same sign, {v}", src)

# (f) independent tests and overlap
src = "results/late_a46/ledger/overlap_audit.tsv"; O = rd(src)
ind, ov = int((O.role == "Independent test").sum()), int((O.role == "Possible overlap").sum())
add("F tests", "independent_tests", ind, 60, "independent target-cohort tests", src)
add("F tests", "possible_overlap", ov, 1, "target-cohort tests with possible overlap that enter the estimates (Wallace SE1)", src)
add("F tests", "total_tests", ind + ov, 61, "target-cohort tests entering the estimates (Table 5 sum of k)", src)
assert ind + ov == int(rd("results/late_a48/ledger/evidence_ledger_a47_v2.tsv").k_independent.sum())

if bad:
    for b in bad: print("MISMATCH", b)
    sys.exit("draft 19 text and result files disagree; nothing written")
head = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
Mf = pd.DataFrame(rows); Mf.insert(0, "repo_commit", head)
Mf["analysis"] = "late analysis, run after all primary results were seen (Amendment 48, wording note for draft 19; reporting only)"
os.makedirs(os.path.dirname(out), exist_ok=True); Mf.to_csv(out + ".tmp", sep="\t", index=False); os.replace(out + ".tmp", out)
print(f"draft 19 numbers at commit {head}: {len(Mf)} rows; every printed value matches its result file")
