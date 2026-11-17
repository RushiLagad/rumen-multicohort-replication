#!/usr/bin/env python3
"""141: Amendment 48.11 (reporting only, no new fit; late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen).
Rebuilds the Additional file 1 evidence ledger from the final Amendment 47 intervals and labels, and adds the 48.8 mapping-confidence classes.
The descriptive columns (source, site, phenotype, published direction, evidence category, phenotype note, cohorts) are carried from the
Amendment 46 ledger (results/late_a46/ledger/evidence_ledger.tsv, script 118), which stays in Additional file 1 as a sheet marked historical.
Estimates: results/late_a47/intervals/pooled_intervals.tsv (version primary): CLR partial r, modified Knapp-Hartung interval for k >= 2 and the
cohort's own interval for k = 1, with the method, role and validation statement. Labels: results/late_a47/labels/labels_{rfi,methane}.tsv
(label_a47) beside the Amendment 37 primary label. Mapping: results/late_a48/mapping/mapping_classes.tsv (MC1 to MC4 per target and cohort).
Writes results/late_a48/ledger/evidence_ledger_a47.tsv (never overwrites). Usage: 141_a48_ledger.py [root]
"""
import os, sys
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
LABEL = "late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen (Amendment 48.11, reporting only)"
od = P("results", "late_a48", "ledger"); out = os.path.join(od, "evidence_ledger_a47.tsv")
if os.path.exists(out): sys.exit(f"output exists, Amendment 48 never overwrites: {out}")
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
L46 = rd("results/late_a46/ledger/evidence_ledger.tsv")
p = rd("results/late_a47/intervals/pooled_intervals.tsv"); p = p[p.version == "primary"].copy()
lab = pd.concat([rd("results/late_a47/labels/labels_rfi.tsv"), rd("results/late_a47/labels/labels_methane.tsv")], ignore_index=True)
mc = rd("results/late_a48/mapping/mapping_classes.tsv")
assert len(L46) == 38 and len(p) == 38 and set(L46.id) == set(p.id), "ledger and Amendment 47 intervals must cover the same 38 targets"
assert (p.designated_method == "mkh").all() and p.validated_procedure.all()

COH = {"nelore": "Nelore", "charolais_mrna": "Charolais", "alberta_all": "Alberta", "difford": "Difford", "wallace": "Wallace"}
keep = ["arm", "id", "target", "set", "rank", "source_site", "source_phenotype", "published_direction", "source", "evidence_category",
        "phenotype_note", "k_independent", "cohorts", "animals_total"]
L = L46[keep].merge(p[["arm", "id", "k", "cohorts", "r", "lo_mkh", "hi_mkh", "p_mkh", "interval", "statement", "stratum", "claim",
                       "excl_15_mkh", "excl_30_mkh", "excl0_mkh", "tau_hat", "published_dir"]].rename(columns={"cohorts": "cohorts_a47"}),
                   on=["arm", "id"], how="left", validate="one_to_one")
assert (L.k == L.k_independent).all(), "k differs between the Amendment 46 ledger and Amendment 47"
lab1 = lab[["arm", "id", "label_a47", "label_a37_primary", "label_note"]]
L = L.merge(lab1, on=["arm", "id"], how="left", validate="one_to_one")
assert L.label_a47.notna().all()

L["method"] = L.interval.map({"pooled": "modified Knapp-Hartung (REML random effects, t on k - 1 df), validated in Amendment 47",
                              "cohort (k = 1)": "single cohort: Fisher z interval of the CLR partial r (Amendment 46.1)"})
L["role"] = L.interval.map({"pooled": "pooled over independent cohorts", "cohort (k = 1)": "one independent cohort"})
L["validation_statement"] = L.statement
L["in_published_direction"] = L.published_dir.map({1.0: "positive", -1.0: "negative"}).fillna("no published direction")
d = L.published_dir.fillna(0)
L["interval_against_published_direction"] = ((d > 0) & (L.hi_mkh < 0)) | ((d < 0) & (L.lo_mkh > 0))
L["excludes_r015"] = L.excl_15_mkh.astype(bool); L["excludes_r030"] = L.excl_30_mkh.astype(bool)
L["excludes_0"] = L.excl0_mkh.astype(bool)
# a descriptive pooled interval carries no exclusion statement (47.1 addendum item 3)
desc = L.statement.str.startswith("descriptive")
assert not (desc & (L.excludes_r015 | L.excludes_r030 | L.interval_against_published_direction)).any()

# 48.8 mapping confidence per cohort (independent cohorts of this target only)
mc["c"] = mc.cohort.map(COH).fillna(mc.cohort)
def mapping(row):
    cs = [c for c in str(row.cohorts_a47).split(";") if c and c != "nan"]
    x = mc[(mc.id == row.id) & (mc.arm == row.arm) & (mc.cohort.isin(cs))]
    assert len(x) == len(cs), f"mapping class missing for {row.arm} {row.id}"
    x = x.set_index("cohort").loc[cs]
    return pd.Series({"mapping_confidence": "; ".join(f"{COH.get(c, c)} {k}" for c, k in zip(cs, x["class"])),
                      "mapping_basis": "; ".join(f"{COH.get(c, c)}: {b}" for c, b in zip(cs, x["basis"])),
                      "mapping_lowest": max(x["class"]) if len(x) else "", "all_MC1": bool(len(x)) and (x["class"] == "MC1").all()})
L = pd.concat([L, L.apply(mapping, axis=1)], axis=1)

O = L[keep[:-2] + ["cohorts_a47", "animals_total", "stratum", "claim", "role", "method", "validation_statement", "r", "lo_mkh", "hi_mkh", "p_mkh",
                   "tau_hat", "in_published_direction", "excludes_0", "interval_against_published_direction", "excludes_r015", "excludes_r030",
                   "label_a47", "label_a37_primary", "label_note", "mapping_confidence", "mapping_lowest", "all_MC1", "mapping_basis"]].rename(columns={
    "cohorts_a47": "cohorts", "r": "r_clr", "lo_mkh": "lo_95", "hi_mkh": "hi_95", "p_mkh": "p", "label_a37_primary": "label_primary_a37"})
O["analysis"] = LABEL
os.makedirs(od, exist_ok=True); tmp = out + f".tmp{os.getpid()}"; O.to_csv(tmp, sep="\t", index=False); os.replace(tmp, out)
print(f"wrote {out}: {len(O)} targets; RFI r = 0.30 exclusions {int((O.arm.eq('RFI') & O.excludes_r030).sum())}, "
      f"methane r = 0.15 exclusions {int((O.arm.eq('Methane') & O.excludes_r015).sum())}; labels {O.label_a47.value_counts().to_dict()}; "
      f"all-MC1 targets {int(O.all_MC1.sum())}")
