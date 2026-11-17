#!/usr/bin/env python3
"""144: Amendment 48.11 source audit and draft 18 reporting files (late analysis, run after all primary results were seen; reporting only, no fit,
no estimate, interval or label changes). Writes new files only, never overwrites:

  results/late_a48/source_audit.tsv                 the review chat's source audit of the methane claims (SOURCE_AUDIT_48_11.md, 17 Nov 2026)
  results/late_a48/ledger/evidence_ledger_a47_v2.tsv the 141 ledger with
        - S1 Sharpea reclassed as a cross-phenotype and cross-host transfer (sheep methane-yield claim tested against cattle daily output),
        - the M1, M2 and S1 source and phenotype columns from the audit (M1 and M2 flagged as group/cluster contrasts),
        - "validated" replaced by "met the Amendment 47 coverage rule" (validation_statement, method) and by "did not meet the Amendment 47
          coverage rule" in label_note (C1);
  results/late_a48/ledger/category_counts_v2.tsv    targets by arm, kind of evidence and number of independent cohorts (draft 18 Table 5),
                                                     and the matched-claim counts (draft 17 vs after 48.11);
  results/late_a48/dairy/not_harmonized_sources.tsv  dairy/not_harmonized.tsv (48.3) with the entries filled from the source papers,
        Difford et al. 2018 PLoS Genet 14:e1007580 and Wallace et al. 2019 Sci Adv 5:eaav8391 (Materials and Methods, read 17 Nov 2026).
Usage: 144_a48_source_audit.py [root]
"""
import os, sys
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
LABEL = "late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen (Amendment 48.11 source audit, reporting only)"
OUT = {"audit": P("results", "late_a48", "source_audit.tsv"), "ledger": P("results", "late_a48", "ledger", "evidence_ledger_a47_v2.tsv"),
       "counts": P("results", "late_a48", "ledger", "category_counts_v2.tsv"), "nh": P("results", "late_a48", "dairy", "not_harmonized_sources.tsv")}
ex = [f for f in OUT.values() if os.path.exists(f)]
if ex: sys.exit(f"output exists, never overwritten: {ex}")
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
def write(df, f):
    df = df.copy(); df["analysis"] = LABEL; tmp = f + f".tmp{os.getpid()}"; df.to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)

# ---------------- source audit (SOURCE_AUDIT_48_11.md)
AUD = pd.DataFrame([
 ("S1", "Sharpea", "Kittelmann et al. 2014, PLoS One 9:e103171", "sheep (New Zealand), 118 high and low emitters",
  "high vs low emitters selected on mean CH4 yield; Sharpea azabuensis characterises the low-CH4 ruminotype S (community types Q, S vs H)",
  "CH4 yield, g per kg DMI", "lower with higher yield", "cross-phenotype and cross-host transfer (no cattle daily-output source)", "transfer; not a matched claim"),
 ("S1", "Sharpea", "Kamke et al. 2016, Microbiome", "sheep, low vs high methane yield", "group contrast", "CH4 yield",
  "enriched in low-yield sheep", "cross-phenotype and cross-host transfer", "transfer; not a matched claim"),
 ("M2", "Succinivibrionaceae", "Danielsson et al. 2017, Front Microbiol 8:226", "dairy cows (Sweden), 21 of 73 (12 Holstein, 9 Swedish Red)",
  "cows chosen as persistent low or high emitters on total mean CH4 production (g/day); Succinivibrionaceae 2.2% in bacterial cluster L vs 0.4% in cluster H; clusters differed in g/d (301 vs 321) and g/kg ECM (8.3 vs 9.7)",
  "CH4, g per day (selection and cluster comparison)", "lower with higher methane", "matched claim, kept; flagged as group/cluster contrast", "matched claim"),
 ("M2", "Succinivibrionaceae", "Wallace et al. 2015, BMC Genomics 16:839", "beef steers (UK), 4 extreme pairs from 72",
  "extremes on CH4 per kg DM intake, respiration chambers; Proteobacteria (88% Succinivibrionaceae) 2.7% vs 11.2% in high vs low",
  "CH4 yield, g per kg DMI", "lower with higher yield", "supporting, yield only", "not counted"),
 ("M2", "Succinivibrionaceae", "Wallace et al. 2019, Sci Adv 5:eaav8391", "the Wallace cohort itself", "", "", "",
  "excluded by the Amendment 32 independence rule (not counted where tested)", "not counted"),
 ("M1", "M. gottschalkii clade vs M. ruminantium clade", "Danielsson et al. 2017, Front Microbiol 8:226", "dairy cows (Sweden), as above",
  "persistent low vs high emitters on g/day: M. ruminantium clade 61.8% vs 40.7% of archaeal sequences in cluster L vs H; M. gottschalkii clade 25.4% vs 44.0%",
  "CH4, g per day", "ratio higher with higher methane", "matched claim, kept; flagged as group contrast", "matched claim"),
 ("M1", "M. gottschalkii clade vs M. ruminantium clade", "Shi et al. 2014, Genome Res 24:1517; Kittelmann et al. 2014", "sheep", "yield contrasts",
  "CH4 yield, g per kg DMI", "as above", "supporting, yield; tested as the Wallace yield sensitivity", "not counted"),
], columns=["id", "target", "source", "host_n", "contrast", "phenotype_units", "direction", "class_48_11", "counted_as"])
AUD.insert(0, "rule", "Amendment 48.11 (committed 4c43bdec): a methane target is a matched claim only if at least one source independent of the test cohort "
                      "reports it in cattle for daily output (g/d) with a stated direction; otherwise a cross-phenotype and, where relevant, cross-host transfer; "
                      "group or community-type contrasts keep their class but are flagged")
AUD["audited_by"] = "review chat, 17 Nov 2026 (SOURCE_AUDIT_48_11.md)"

# ---------------- ledger v2
L = rd("results/late_a48/ledger/evidence_ledger_a47.tsv"); assert len(L) == 38
n_valid = int((L.validation_statement == "validated").sum())
L["validation_statement"] = L.validation_statement.replace({"validated": "met the Amendment 47 coverage rule"})
L["method"] = L.method.str.replace("validated in Amendment 47", "which met the Amendment 47 coverage rule", regex=False)
L["label_note"] = L.label_note.fillna("").str.replace("configuration not validated", "configuration did not meet the Amendment 47 coverage rule", regex=False)
assert not L.validation_statement.eq("validated").any() and not L.method.str.contains("validated").any() and not L.label_note.str.contains("validated").any()
L["source_contrast_flag"] = ""
upd = {
 "S1": dict(evidence_category="Cross-phenotype transfer (methane yield claim in sheep tested against daily output)",
            source="Kittelmann et al. 2014, PLoS One 9:e103171 (sheep, New Zealand, 118 high and low emitters selected on CH4 yield; Sharpea azabuensis characterises the low-CH4 ruminotype S); Kamke et al. 2016, Microbiome (sheep, low vs high methane yield)",
            source_phenotype="CH4 yield (g per kg DMI), sheep",
            phenotype_note="Claim published for methane yield in sheep (high vs low emitter groups); tested against daily output (g per day) in dairy cows: cross-phenotype and cross-host transfer (Amendment 48.11 source audit); yield tested in Wallace as a sensitivity (48.1)",
            source_contrast_flag="group contrast (high vs low emitters; community types)"),
 "M2": dict(source="Danielsson et al. 2017, Front Microbiol 8:226 (dairy cows, Sweden; 21 of 73 persistent low or high emitters selected on CH4 g/day; Succinivibrionaceae 2.2% in bacterial cluster L vs 0.4% in cluster H); supporting, yield only: Wallace et al. 2015, BMC Genomics 16:839; Wallace et al. 2019 (the Wallace cohort) excluded by the Amendment 32 independence rule",
            source_phenotype="CH4, g per day (dairy cows)",
            phenotype_note="Claim published for daily output in dairy cows as a group/cluster contrast; tested as a linear association with g per day (a different estimand); matched claim kept, flagged (Amendment 48.11 source audit)",
            source_contrast_flag="group/cluster contrast (persistent low vs high emitters)"),
 "M1": dict(source="Danielsson et al. 2017, Front Microbiol 8:226 (dairy cows, Sweden; persistent low vs high emitters on CH4 g/day: M. ruminantium clade 61.8% vs 40.7% and M. gottschalkii clade 25.4% vs 44.0% of archaeal sequences in cluster L vs H); supporting, yield: Shi et al. 2014, Genome Res 24:1517 and Kittelmann et al. 2014 (sheep)",
            source_phenotype="CH4, g per day (dairy cows); yield in sheep (supporting)",
            phenotype_note="Claim published for daily output in dairy cows as a group contrast (Danielsson) and for yield in sheep (Shi, Kittelmann); tested as a linear association with g per day (primary), yield as the Wallace sensitivity; matched claim kept, flagged (Amendment 48.11 source audit)",
            source_contrast_flag="group contrast (persistent low vs high emitters)")}
for i, u in upd.items():
    m = L.id == i; assert m.sum() == 1, i
    for c, v in u.items(): L.loc[m, c] = v
L["source_audit"] = L.id.map(lambda i: "Amendment 48.11 source audit (results/late_a48/source_audit.tsv)" if i in upd else "")
assert L.drop(columns=["analysis"]).shape[1] == rd("results/late_a48/ledger/evidence_ledger_a47.tsv").drop(columns=["analysis"]).shape[1] + 2

def kind(c):
    for k in ("Matched claim", "Published claim without direction", "Candidate taxon", "Cross-site transfer", "Cross-phenotype"):
        if c.startswith(k): return "Cross-phenotype transfer" if k == "Cross-phenotype" else k
    raise ValueError(c)
L["kind"] = L.evidence_category.map(kind)
C = (L.groupby(["arm", "kind", "k_independent"]).size().unstack(fill_value=0).reindex(columns=[1, 2, 3], fill_value=0))
C.columns = ["k = 1", "k = 2", "k = 3"]; C["Total"] = C.sum(axis=1); C = C.reset_index()
tot = pd.DataFrame([{"arm": "All", "kind": "", **{c: int(C[c].sum()) for c in ["k = 1", "k = 2", "k = 3", "Total"]}}])
C = pd.concat([C, tot], ignore_index=True); C.insert(0, "table", "targets by kind of evidence and k (draft 18 Table 5)")
mat = L[L.kind == "Matched claim"]
CNT = pd.DataFrame([
 ("matched claims", 21, len(mat)), ("matched with k >= 2 independent cohorts", 7, int((mat.k_independent >= 2).sum())),
 ("matched with k = 1", 14, int((mat.k_independent == 1).sum())), ("other targets", 17, int(len(L) - len(mat))),
 ("cross-phenotype transfers (methane)", 7, int(((L.arm == "Methane") & (L.kind == "Cross-phenotype transfer")).sum()))],
 columns=["count", "draft_17", "after_48_11"])
assert list(CNT.after_48_11) == [20, 6, 14, 18, 8], CNT
CNT.insert(0, "table", "matched-claim counts (SOURCE_AUDIT_48_11.md)")
counts = pd.concat([C, CNT], ignore_index=True)

# ---------------- not_harmonized filled from the source papers
NH = rd("results/late_a48/dairy/not_harmonized.tsv")
D18 = "Difford et al. 2018, PLoS Genet 14:e1007580, Materials and methods"; W19 = "Wallace et al. 2019, Sci Adv 5:eaav8391, Materials and Methods"
FILL = {
 "primer set and amplified region": (
   "bacteria 27F/534R, V1-V3 (expected 508 bp); archaea S-D-Arch-0519-a-S-15/S-D-Arch-1041-a-A-18, V4-V6 (542 bp); Illumina 300 bp paired-end, about half on MiSeq and half on HiSeq",
   "primers designed in silico with ecoPrimers/OBITools (table S2 of the source); prokaryotic OTUs from the ProkA set (bacteria and archaea); MiSeq 250 bp paired-end (archaeal marker: HiSeq 100 bp paired-end)"),
 "read length and read use": (
   "300 bp paired-end sequenced; this re-analysis used forward reads only (DADA2, script 73)",
   "250 bp paired-end sequenced; processed as in script 79"),
 "breed composition": (
   "Holstein (750 lactating cows)",
   "Holstein (816 cows; UK1, UK2, IT1 to IT3) and Nordic Red (200 cows; SE1, FI1)"),
 "production system and diet": (
   "Denmark; five commercial herds and one research herd (DCRC); indoor feeding, ad libitum TMR (rolled barley, corn silage, grass-clover silage, rapeseed and soybean meal) plus up to 3 kg concentrate during automatic milking",
   "UK, Italy, Sweden, Finland; seven farms; diets standardized within farm: maize silage, grass silage or grass hay and concentrates (UK, IT), grass silage and concentrates (SE, FI); ad libitum TMR (IT, SE, FI) or PMR plus concentrate in robotic milking (UK); cows 10 to 40 weeks postpartum"),
 "methane measurement method": (
   "breath CH4 and CO2 concentrations (ppm) at the automatic milking system feed bin: portable FTIR (Gasmet DX-4000) for 7 days in the commercial herds, permanently installed NDIR (Guardian NG) in the research herd; daily CH4 (g/day) from the CH4:CO2 ratio and predicted CO2 output",
   "breath sampling during robotic milking (UK); GreenFeed bait station (IT, SE); respiration chambers for 5 days (FI); differs by country, so CH4 is z-scored within farm in this re-analysis"),
 "sampling (rumen fluid collection)": (
   "oral probe ('Flora Rumen Scoop'), about 40 mL of rumen liquid with particulate matter, immediately after the methane recording week",
   "oral ruminal probe (ruminator), first litre discarded, 0.5 L retained; day 1 of the sampling week, 2 to 5 hours after morning feeding; same method at all centres"),
}
assert set(FILL) == set(NH["item"]), set(NH["item"]) ^ set(FILL)
NH["difford"] = NH["item"].map(lambda k: FILL[k][0]); NH["wallace"] = NH["item"].map(lambda k: FILL[k][1])
NH["source_difford"] = D18; NH["source_wallace"] = W19
NH["note"] = "filled from the source papers for draft 18 (Amendment 48.3 'not harmonizable'); tabulated, not modelled"
NH = NH.drop(columns=["analysis"])

write(AUD, OUT["audit"]); write(L.drop(columns=["kind", "analysis"]), OUT["ledger"]); write(counts, OUT["counts"]); write(NH, OUT["nh"])
print(f"wrote {OUT['audit']} ({len(AUD)} rows)"); print(f"wrote {OUT['ledger']}: S1 reclassed; M1, M2, S1 sources from the audit; 'validated' replaced in {n_valid} rows")
print(f"wrote {OUT['counts']}"); print(C.drop(columns=["table"]).to_string(index=False)); print(CNT.drop(columns=["table"]).to_string(index=False))
print(f"wrote {OUT['nh']} ({len(NH)} items)")
