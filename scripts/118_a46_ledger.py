#!/usr/bin/env python3
"""118: Amendment 46 (late analysis, run after all primary results were seen; reporting only, no new fit), review action items A1 and A4.
Evidence ledger: one row per testable target with its source claim, site, phenotype, published direction, an evidence category, the independent
cohorts that test it (k) and the Amendment 46 estimate. Overlap audit: one row per target and cohort, with the cohort's role (independent test,
discovery cohort, re-analysis, possible overlap) and whether any primary estimate is affected. Facts come from plan/target_keys.tsv, the fit
files and the amendments cited in each row (14, 19, 32, 39 and the M1 source-trait record). Writes results/late_a46/ledger/ (never overwrites).
Run: python3 scripts/118_a46_ledger.py
"""
import os, glob
import pandas as pd

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort"); P = lambda *a: os.path.join(ROOT, *a)
LABEL = "late analysis, run after all primary results were seen (Amendment 46, review items A1 and A4)"
od = P("results", "late_a46", "ledger"); os.makedirs(od, exist_ok=True)
outs = [os.path.join(od, f) for f in ("evidence_ledger.tsv", "overlap_audit.tsv")]
ex = [f for f in outs if os.path.exists(f)]; assert not ex, f"output exists, Amendment 46 never overwrites: {ex}"
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
tk = rd("plan/target_keys.tsv").set_index("id")
pool = rd("results/late_a46/effect_size/pooled.tsv"); coh = rd("results/late_a46/effect_size/cohort_estimates.tsv")

SOURCE = {   # as in scripts/104_additional_files.py
    "Li": "Li et al. 2019, Microbiome 7:6 (PRJNA448333 steers; 8 high vs 8 low RFI per breed, 16S, uncorrected t tests; plan 3.1)",
    "Nelore": "Conteville et al. 2023, Front Microbiol 14:1282851, Supplementary Table 7 (the Nelore cohort PRJNA987743 itself; continuous RFI_GC or RME_GC)",
    "Canonical": "Recurrent candidates named in the authors' 2026 review; no single source claim, no direction (plan 3.3)",
    "MethaneSecondary": "Methane secondary set, reported in a single region (Amendment 19)"}
SOURCE_ID = {
    "M1": "Danielsson et al. 2017 (daily output, Swedish dairy cows) and Shi et al. 2014 / Kittelmann et al. 2014 (yield, New Zealand sheep)",
    "M2": "Methane primary set, reported in more than one country (Amendment 19); includes Wallace et al. 2019 for the Wallace cohort",
    "M3": "Methane primary set, literature split on direction (Amendment 19)"}

def category(arm, i, t):
    site, ph, d = t["site"], t["phenotype"], str(t["direction"])
    if t["set"] == "Canonical": return "Candidate taxon (no single source claim, no direction)"
    if arm == "Methane" and ph == "RME_GC": return "Cross-phenotype transfer (residual methane claim tested against daily output)" + ("; fecal source" if site == "feces" else "")
    if site == "feces": return "Cross-site transfer (fecal claim tested in rumen)"
    if i == "M3" or not any(k in d for k in ("higher", "lower")): return "Published claim without direction"
    return "Matched claim (rumen, same trait type, published direction)"
def phen_note(arm, i, t):
    if arm == "RFI": return "RFI in source and test; Charolais tested as high vs low RFI class (extreme groups)" + ("; Li source used 8 vs 8 extreme-RFI steers per breed" if t["set"] == "Li" else "")
    if i == "M1": return "Claim published for daily output (Danielsson) and yield (Shi); tested as g per day (primary), yield as a sensitivity (M1 source-trait record)"
    if t["phenotype"] == "RME_GC": return "Source: residual methane emission (RME_GC, Nelore); tested: methane g per day"
    return "Methane trait as reported in the sources (Amendment 19); tested: methane g per day"

L = []
for _, p in pool.iterrows():
    t = tk.loc[p.id]; c = coh[(coh.id == p.id) & (coh.arm == p.arm)]
    L.append(dict(arm=p.arm, id=p.id, target=t["target"], set=t["set"], rank=t["target_rank"], source_site=t["site"], source_phenotype=t["phenotype"],
                  published_direction=t["direction"], source=SOURCE_ID.get(p.id) or SOURCE.get(t["set"], ""),
                  evidence_category=category(p.arm, p.id, t), phenotype_note=phen_note(p.arm, p.id, t),
                  k_independent=int(p.k), cohorts=p.cohorts, animals_total=int(c.n.sum()),
                  r_clr=p.r, lo_95=p.lo, hi_95=p.hi, interval=p.interval, excludes_r030=p.excl_r30, excludes_r015=p.excl_r15,
                  label_primary=p.label_a37, analysis=LABEL))
L = pd.DataFrame(L)

# ---- overlap audit: every testable target x cohort, plus the pairs excluded or re-analysed by rule
ASSAY = {"nelore": "shotgun metagenome, GTDB r207", "charolais_mrna": "metatranscriptome, GTDB r207", "alberta_all": "16S rRNA gene, SILVA 138.1",
         "difford": "16S rRNA gene, SILVA 138.1", "wallace": "16S rRNA gene, SILVA 138.1"}
noSE1 = rd("results/methane_pooled_noSE1/estimates_calibrated.tsv")
alb = pd.concat([rd(f) for f in glob.glob(P("results", "alberta", "alberta_*_all_targets.tsv"))])
li_in_alberta = sorted(alb[(alb.set == "Li") & (alb.status == "ok")].id)
O = []
def row(arm, i, cohort, role, in_primary, note, sens=""):
    O.append(dict(arm=arm, id=i, target=tk.loc[i, "target"], cohort=cohort, assay=ASSAY.get(cohort, ""), role=role, in_primary_estimate=in_primary,
                  basis=note, sensitivity=sens, analysis=LABEL))
for _, c in coh.iterrows():
    i, k = c.id, c.cohort; t = tk.loc[i]
    if c.arm == "RFI" and k == "alberta_all" and t["set"] == "Li":
        row(c.arm, i, k, "Possible overlap", "yes", "Up to 48 of 567 Alberta animals are Li et al. discovery animals and cannot be identified (Amendment 39)", "Alberta 529 version (superseded; Amendment 39)")
    elif c.arm == "Methane" and k == "wallace" and i == "M1":
        x = noSE1[(noSE1.id == "M1") & (noSE1.cohort == "wallace")].iloc[0]
        row(c.arm, i, k, "Possible overlap", "yes", "Swedish farm SE1 may be the Lovsta herd of Danielsson et al. 2017, an M1 source (Amendment 32 item 2)",
            f"Wallace without SE1: r = {x.r:.3f} (95% CI {x.lo:.3f} to {x.hi:.3f}), ANCOM-BC2 calibrated version")
    elif c.arm == "Methane" and k == "difford" and i in ("M2", "S1"):
        row(c.arm, i, k, "Independent test", "yes", "Difford et al. 2018 S1 Table lists estimates for this taxon, but it did not contribute to target choice (Amendment 32 item 1)")
    else:
        row(c.arm, i, k, "Independent test", "yes", "No shared animals with the source study; cohort inclusion rule (Amendment 19 item 1)")
# pairs not in any primary estimate, by rule
for i in sorted(L[(L.arm == "RFI") & (L.set == "Nelore")].id):
    row("RFI", i, "nelore", "Discovery cohort", "no", "The target was reported in this cohort (Conteville et al. 2023); Nelore-set targets are not tested in Nelore (plan 3.2)")
row("Methane", "M2", "wallace", "Re-analysis", "no", "Wallace et al. 2019 reported Succinivibrionaceae with methane in this cohort; re-analysis enters only a sensitivity pooling (Amendment 32 item 1)",
    "Methane with re-analysis version (Additional file 2)")
for i in sorted(L[(L.arm == "RFI") & (L.set == "Li")].id):
    if i not in li_in_alberta:
        row("RFI", i, "alberta_all", "Possible overlap, no estimate", "no", "Up to 48 Alberta animals are Li et al. discovery animals (Amendment 39), but this target has no tested feature in the Alberta tables, so no estimate is affected")
O = pd.DataFrame(O).sort_values(["arm", "id", "cohort"])
L.to_csv(outs[0], sep="\t", index=False); O.to_csv(outs[1], sep="\t", index=False)
print(LABEL)
print("== evidence categories (testable targets)"); print(L.groupby(["arm", "evidence_category"]).size().to_string())
print("== overlap audit roles"); print(O.groupby(["arm", "role", "in_primary_estimate"]).size().to_string())
print("Li targets with an Alberta estimate:", li_in_alberta or "none")
