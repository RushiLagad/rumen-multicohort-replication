#!/usr/bin/env python3
"""143: the numbers of draft 18 Tables 4 and 5, each computed from a committed result file (no fit, no hand-typed value), with the value as the
table prints it. Late analysis, run after all primary results were seen (Amendment 48.11, reporting only).
  Table 4 (claims testable in two independent cohorts): results/late_a48/seven/seven_claims.tsv and seven_pooled.tsv (48.1).
  Table 5 (testable targets by kind of evidence and k): the draft 18 evidence ledger, results/late_a48/ledger/evidence_ledger_a47_v2.tsv (141),
  in which Sharpea (S1) is a cross-phenotype and cross-host transfer; Table 5 counts it with the cross-phenotype transfers.
Each printed value is checked against the draft 18 table text (PRINTED below); the script stops on any mismatch.
Run: python3 scripts/143_manifest_draft18.py [root]   (pandas). Writes results/manifest/draft18_numbers.tsv (never overwrites).
"""
import os, sys, subprocess
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
out = P("results", "manifest", "draft18_numbers.tsv")
if os.path.exists(out): sys.exit(f"output exists, never overwritten: {out}")
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
rows = []
def add(table, key, value, desc, src, col):
    rows.append(dict(table=table, key=key, printed=value, description=desc, source=src, column=col))

def f2(x):   # two decimals, typographic minus, as printed
    s = f"{x:.2f}"; s = "0.00" if s == "-0.00" else s
    return s.replace("-", "−")
COH = {"nelore": "Nelore", "charolais_mrna": "Charolais", "alberta_all": "Alberta", "difford": "Difford", "wallace": "Wallace"}
PAT = {"A": "one cohort", "B": "against", "C": "bounded", "D": "discordant", "E": "inconclusive"}

# ---------------- Table 4
S = rd("results/late_a48/seven/seven_claims.tsv"); Pp = rd("results/late_a48/seven/seven_pooled.tsv")
ORDER = ["T2", "T4", "T5", "T6", "N15", "M1", "S1"]
PRINTED4 = {   # draft 18 Table 4, columns 2 to 5
 "T2": ("Nelore 0.17 (−0.12 to 0.43); Charolais −0.18 (−0.52 to 0.21)", "0.02 (−0.98 to 0.98), 0.18", "discordant / against", "none"),
 "T4": ("Nelore −0.15 (−0.41 to 0.14); Charolais −0.55 (−0.77 to −0.22)", "−0.35 (−1.00 to 0.99), 0.28", "against / against", "none"),
 "T5": ("Nelore 0.13 (−0.16 to 0.39); Charolais −0.30 (−0.61 to 0.08)", "−0.07 (−0.99 to 0.99), 0.25", "discordant / discordant", "none"),
 "T6": ("Nelore −0.13 (−0.39 to 0.16); Charolais −0.06 (−0.42 to 0.32)", "−0.10 (−0.92 to 0.88), 0.00", "inconclusive / discordant", "none"),
 "N15": ("Charolais −0.04 (−0.41 to 0.34); Alberta 0.03 (−0.05 to 0.12)", "0.03 (−0.46 to 0.51), 0.00", "discordant / one cohort", "ANCOM-BC2"),
 "M1": ("Wallace 0.03 (−0.04 to 0.10); Difford 0.03 (−0.04 to 0.11)", "0.03 (−0.28 to 0.34), 0.00", "bounded / bounded", "none"),
 "S1": ("Difford 0.03 (−0.04 to 0.10); Wallace −0.17 (−0.23 to −0.11)", "−0.07 (−0.87 to 0.83), 0.14", "one cohort / bounded", "CLR")}
src_c, src_p = "results/late_a48/seven/seven_claims.tsv", "results/late_a48/seven/seven_pooled.tsv"
bad = []
for i in ORDER:
    c = S[S.id == i]; p = Pp[Pp.id == i]; assert len(c) == 2 and len(p) == 1, i; p = p.iloc[0]
    co = sorted(c.itertuples(), key=lambda r: p.cohorts.split(";").index(r.cohort))
    for r in co:
        add("Table 4", f"{i}_{r.cohort}_clr", f"{f2(r.r_clr)} ({f2(r.lo_clr)} to {f2(r.hi_clr)})", f"{i} {COH[r.cohort]}: CLR partial r (95% Fisher interval), n {r.n}", src_c, "r_clr, lo_clr, hi_clr")
    cell2 = "; ".join(f"{COH[r.cohort]} {f2(r.r_clr)} ({f2(r.lo_clr)} to {f2(r.hi_clr)})" for r in co)
    cell3 = f"{f2(p.r_mkh)} ({f2(p.lo_mkh)} to {f2(p.hi_mkh)}), {p.tau_hat:.2f}"
    add("Table 4", f"{i}_pooled", cell3, f"{i}: pooled r (modified Knapp-Hartung 95% interval), tau-hat", src_p, "r_mkh, lo_mkh, hi_mkh, tau_hat")
    cell4 = f"{PAT[p.pattern_clr[0]]} / {PAT[p.pattern_ancombc2[0]]}"
    add("Table 4", f"{i}_pattern", cell4, f"{i}: pattern under CLR / calibrated ANCOM-BC2 (48.1 rule)", src_p, "pattern_clr, pattern_ancombc2")
    sa = [n for n, v in (("CLR", p.secondary_Sa_clr), ("ANCOM-BC2", p.secondary_Sa_ancombc2)) if str(v) == "True"]
    cell5 = ", ".join(sa) if sa else "none"
    add("Table 4", f"{i}_secondary", cell5, f"{i}: secondary criterion S-a met (nominal single-cohort support, reported only)", src_p, "secondary_Sa_clr, secondary_Sa_ancombc2")
    add("Table 4", f"{i}_reasons", p.reasons_not_replicated, f"{i}: reasons for no replicated label (Section 10 conditions)", src_p, "reasons_not_replicated")
    for got, want, what in ((cell2, PRINTED4[i][0], "cohort CLR"), (cell3, PRINTED4[i][1], "pooled"), (cell4, PRINTED4[i][2], "pattern"), (cell5, PRINTED4[i][3], "secondary")):
        if got != want: bad.append((i, what, got, want))
allfail = all(("pooled q >= 0.05" in r and "leave-one-cohort-out" in r) for r in Pp.reasons_not_replicated)
add("Table 4", "all_fail_q_and_loo", str(allfail), "every claim failed on pooled q >= 0.05 and the leave-one-cohort-out condition (footnote)", src_p, "reasons_not_replicated")
if not allfail: bad.append(("all", "footnote", allfail, True))

# ---------------- Table 5
src_l = "results/late_a48/ledger/evidence_ledger_a47_v2.tsv"; L = rd(src_l)
assert len(L) == 38
def kind(c):
    for k in ("Matched claim", "Published claim without direction", "Candidate taxon", "Cross-site transfer", "Cross-phenotype"):
        if c.startswith(k): return "Cross-phenotype transfer" if k == "Cross-phenotype" else k
    raise ValueError(c)
L["kind"] = L.evidence_category.map(kind)
assert L.loc[L.id == "S1", "kind"].iloc[0] == "Cross-phenotype transfer", "S1 must be a transfer in the draft 18 ledger"
PRINTED5 = {("RFI", "Matched claim"): (13, 0 + 5, 0, 18), ("RFI", "Published claim without direction"): (0, 1, 0, 1), ("RFI", "Candidate taxon"): (0, 1, 5, 6),
            ("RFI", "Cross-site transfer"): (1, 1, 0, 2), ("Methane", "Matched claim"): (1, 1, 0, 2), ("Methane", "Published claim without direction"): (1, 0, 0, 1),
            ("Methane", "Cross-phenotype transfer"): (4, 4, 0, 8), ("All", ""): (20, 13, 5, 38)}
for (arm, k), want in PRINTED5.items():
    x = L if arm == "All" else L[(L.arm == arm) & (L.kind == k)]
    got = tuple(int((x.k_independent == j).sum()) for j in (1, 2, 3)) + (len(x),)
    ids = ";".join(x.sort_values(["k_independent", "id"]).id) if arm != "All" else ""
    for j, lab in enumerate(("k = 1", "k = 2", "k = 3", "Total")):
        add("Table 5", f"{arm}|{k}|{lab}", got[j], f"{arm}, {k or 'all targets'}: {lab}" + (f" (targets: {ids})" if lab == "Total" and ids else ""), src_l, "arm, evidence_category, k_independent")
    if got != want: bad.append((arm, k, got, want))
nel = L[(L.arm == "Methane") & (L.kind == "Cross-phenotype transfer")]
add("Table 5", "methane_transfer_composition", f"{int((nel.set == 'Nelore').sum())} Nelore genera + {int((nel.id == 'S1').sum())} Sharpea", "footnote: methane cross-phenotype transfers", src_l, "set, id")

if bad:
    for b in bad: print("MISMATCH", b)
    sys.exit("draft 18 table text and result files disagree; nothing written")
head = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
M = pd.DataFrame(rows); M.insert(0, "repo_commit", head)
M["analysis"] = "late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen (Amendment 48.11, reporting only)"
os.makedirs(os.path.dirname(out), exist_ok=True); M.to_csv(out + ".tmp", sep="\t", index=False); os.replace(out + ".tmp", out)
print(f"draft 18 Tables 4 and 5 at commit {head}: {len(M)} rows; every printed value matches its result file")
