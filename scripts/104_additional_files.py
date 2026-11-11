#!/usr/bin/env python3
"""104: Additional files 1 and 2 for the manuscript, built from existing outputs only (no fits).

Additional file 1: every pre-specified target with its source, published direction, mapping, and per cohort the
mapped feature, prevalence and fit status. Prevalence is not stored in any results table: scripts 21, 65 and 81
print it to their run logs, so it is parsed from those logs and checked against each fit's status.
Additional file 2: per-cohort calibrated estimates (W / max(1, null SD), permutation p) and pooled results for every
analysis version found by scripts 97 and 99, with the primary versions marked.

Run: python3 scripts/104_additional_files.py   (needs pandas and openpyxl; RUMEN_ROOT defaults to /data/rumen_multicohort)
Writes results/additional_files/Additional_file_1_targets.xlsx and Additional_file_2_estimates.xlsx.
"""
import os, re, math
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort")
OD = os.path.join(ROOT, "results", "additional_files"); os.makedirs(OD, exist_ok=True)
P = lambda *a: os.path.join(ROOT, *a)
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})

# ---------------------------------------------------------------- coverage tables printed in the run logs
COLS = {"id", "target", "set", "rank", "n_keys", "n_present", "feature", "reads", "prevalence", "pass10"}
def _hdr(l):
    t = l.split(); return bool(t) and set(t) <= COLS
def parse_cov(path):
    """R print.data.frame output (right-aligned, wrapped into column blocks) following the 'targets N' line."""
    L = open(P(path)).read().splitlines()
    s = next(i for i, l in enumerate(L) if re.search(r"targets \d+", l)); i = s + 1
    assert _hdr(L[i]), (path, L[i])
    j = i + 1
    while j < len(L) and not _hdr(L[j]) and L[j].strip() and not L[j].startswith(("Loading", "ANCOM", "There were", "Warning")):
        j += 1
    n = j - i - 1; rec = {}
    while i < len(L) and _hdr(L[i]):
        ends = [m.end() for m in re.finditer(r"\S+", L[i])]; names = L[i].split()
        for r, row in enumerate(L[i + 1:i + 1 + n]):
            assert not _hdr(row), (path, row)
            d = rec.setdefault(r, {}); start = 0
            for nm, e in zip(names, ends): d[nm] = row[start:e].strip(); start = e
        i += 1 + n
    df = pd.DataFrame([rec[k] for k in sorted(rec)])
    df["prevalence"] = pd.to_numeric(df["prevalence"].replace({"NA": None, "NaN": None}), errors="coerce")
    df["feature"] = df["feature"].replace({"<NA>": None})
    return df

# per cohort: (column label, [(coverage log, targets file), ...]) ; Alberta = all-animal analysis (primary, Amendment 39)
RANKS = ["species", "genus", "family", "phylum"]
COHORTS = [
    ("Nelore (shotgun, RFI)", [(f"logs/replication/nelore_{r}_fit.log", f"results/replication/nelore_{r}_targets.tsv") for r in RANKS]),
    ("Charolais (mRNA, RFI)", [(f"logs/replication/charolais_mrna_{r}_fit.log", f"results/replication/charolais_mrna_{r}_targets.tsv") for r in RANKS]),
    ("Alberta (16S, RFI, all animals)", [("logs/alberta/genus_mode1.log", "results/alberta/alberta_genus_all_targets.tsv"),
                                         ("logs/alberta_family_silva_all.log", "results/alberta/alberta_family_all_targets.tsv")]),
    ("Difford (16S, CH4)", [("logs/methane/difford_bact.log", "results/methane/difford_bact_targets.tsv"),
                            ("logs/methane/difford_arch.log", "results/methane/difford_arch_targets.tsv")]),
    ("Wallace (16S, CH4)", [("logs/methane/wallace_prokaryotes.log", "results/methane/wallace_prokaryotes_targets.tsv")]),
    ("Wallace (16S, RFI NRC, sensitivity)", [("logs/wallace_rfi_real.log", "results/methane/wallace_prokaryotes_rfi_nrc_targets.tsv")]),
]
STATUS_TXT = {"ok": "tested", "filtered": "below 10% prevalence, not tested", "no_feature": "no matching feature",
              "untestable_here": "not resolvable in this table"}

def cohort_cov(pairs):
    """One row per target id: feature, prevalence, status. For a target in several tables, the tested one wins."""
    out = {}
    for lg, tf in pairs:
        t = rd(tf); c = parse_cov(lg).set_index("id")
        for _, r in t.iterrows():
            cv = c.loc[r["id"]] if r["id"] in c.index else None
            prev = cv["prevalence"] if cv is not None else float("nan")
            feat = cv["feature"] if cv is not None else None
            # consistency with the fit: tested => passed 10% prevalence; filtered => below it
            if r["status"] == "filtered": assert prev < 0.10, (tf, r["id"], prev)
            if r["status"] == "ok" and r["id"] != "M1": assert prev >= 0.10, (tf, r["id"], prev)
            rec = dict(feature=feat if feat else (r["feature"] if isinstance(r["feature"], str) else None),
                       prevalence=prev, status=r["status"], n=r["n"], table=os.path.basename(tf).replace("_targets.tsv", ""))
            old = out.get(r["id"])
            if old is None or (old["status"] != "ok" and r["status"] == "ok"): out[r["id"]] = rec
    return out

# ---------------------------------------------------------------- sources (from plan/ANALYSIS_PLAN.md)
SOURCE = {
    "Li": "Li et al. 2019, Microbiome 7:6 (PRJNA448333 steers); plan section 3.1",
    "Nelore": "Conteville et al. 2023, Front Microbiol 14:1282851 (PMC10754987), Supplementary Table 7; plan section 3.2",
    "Canonical": "Recurrent candidates named in the authors' 2026 review; no single source paper, no direction; plan section 3.3",
}
SOURCE_ID = {
    "M1": "Shi et al. 2014, Genome Res 24:1517; Kittelmann et al. 2014; Danielsson et al. 2017, Front Microbiol 8:226 (plan, M1 source traits, 2026-11-10)",
    "M2": "Methane primary set, reported in more than one country (Amendment 19); Wallace et al. 2019 reported it in the Wallace cohort, so it is re-analysis there (Amendment 32)",
    "M3": "Methane primary set, literature split on direction, tested two-sided (Amendment 19)",
}
SOURCE_SET_M = {"MethaneSecondary": "Methane secondary set, reported in a single region (Amendment 19)"}
SOURCES_CHECKED = ("Sources checked for the methane targets (Amendment 32): Shi et al. 2014 and Kittelmann et al. 2014 (New Zealand sheep), "
                   "Danielsson et al. 2017 (Swedish dairy cows), Wallace et al. 2015 (UK beef), Kamke et al. 2016 (New Zealand sheep), "
                   "Wallace et al. 2019, Difford et al. 2018 and Zhang et al. 2020.")

def build_file1():
    tk = rd("plan/target_keys.tsv")
    cov = [(lab, cohort_cov(pairs)) for lab, pairs in COHORTS]
    rows = []
    for _, t in tk.iterrows():
        src = SOURCE_ID.get(t["id"]) or SOURCE.get(t["set"]) or SOURCE_SET_M.get(t["set"], "")
        r = {"Set": t["set"], "ID": t["id"], "Target": t["target"], "Rank": t["target_rank"], "Site": t["site"],
             "Published trait": t["phenotype"], "Published direction": t["direction"],
             "Published coefficient": pd.to_numeric(t["coef"], errors="coerce"), "Published SE": pd.to_numeric(t["stderr"], errors="coerce"),
             "Published q or p": pd.to_numeric(t["qval"], errors="coerce"), "Source": src,
             "Mapping status": t["status"], "Reference keys": t["gtdb_keys"] if isinstance(t["gtdb_keys"], str) else "",
             "Mapping note": t["note"] if isinstance(t["note"], str) else ""}
        for lab, c in cov:
            x = c.get(t["id"])
            if x is None:
                r[f"{lab}: feature"] = ""; r[f"{lab}: prevalence"] = None; r[f"{lab}: status"] = "not in this cohort's target list"
            else:
                feat = x["feature"] or ""
                if t["id"] == "M1" and x["status"] == "ok": feat = f"clade log-ratio (cows with >=100 assigned reads used)"
                r[f"{lab}: feature"] = feat
                r[f"{lab}: prevalence"] = None if (isinstance(x["prevalence"], float) and math.isnan(x["prevalence"])) else x["prevalence"]
                r[f"{lab}: status"] = STATUS_TXT.get(x["status"], x["status"])
        rows.append(r)
    T = pd.DataFrame(rows)
    notes = pd.DataFrame({"Item": [
        "Contents", "Prevalence", "Feature", "Status", "Cohorts", "Reference keys", "Published values",
        "Methane sources", "Built by", "Built from"],
        "Description": [
        "All pre-specified targets (plan/target_keys.tsv): Li set (6), Nelore set (64), Canonical set (6), Methane set (3), Methane secondary set (3).",
        "Fraction of animals with at least one read of the mapped feature, before fitting. Targets below 0.10 were not tested (ANCOM-BC2 prv_cut 0.10). Blank: no matching feature.",
        "Name of the feature in the cohort table. TGT_<id> is a composite formed by summing the reference keys listed for that target.",
        "tested; below 10% prevalence, not tested; no matching feature; not resolvable in this table (M1 and M3 where the amplified region or table cannot resolve them); not in this cohort's target list (for example Nelore-set targets in the Nelore discovery cohort).",
        "Nelore PRJNA987743 (52 animals), Charolais PRJEB70895 mRNA (30), Alberta PRJNA393057 16S, all-animal analysis (567; primary per Amendment 39), Difford PRJEB28065 16S (749 cows; M1 from the archaeal table), Wallace PRJNA517480 16S (1,009 cows), and Wallace RFI NRC (1,004 cows, sensitivity only, Amendment 43).",
        "GTDB r207 keys for the shotgun and mRNA cohorts; the 16S cohorts match SILVA 138.1 genus names exactly (Amendments 17, 35, 36).",
        "As reported in the source (plan/target_keys.tsv). The Nelore set has coefficients and q-values from Supplementary Table 7; the Li set has p-values only.",
        SOURCES_CHECKED + " Per-target citations for M2, M3 and the secondary set are given at set level in the plan.",
        "scripts/104_additional_files.py",
        "plan/target_keys.tsv; coverage tables printed by scripts 21, 65 and 81 in logs/replication, logs/alberta, logs/methane and logs/wallace_rfi_real.log; fit status from results/replication, results/alberta and results/methane *_targets.tsv."]})
    f = os.path.join(OD, "Additional_file_1_targets.xlsx")
    with pd.ExcelWriter(f, engine="openpyxl") as w:
        notes.to_excel(w, "Notes", index=False); T.to_excel(w, "Targets", index=False)
    style(f, {"Targets": 3}); return f, T

# ---------------------------------------------------------------- Additional file 2
def partial_r(W, df): return W / math.sqrt(W * W + df)
def ci(r, df):
    z = math.atanh(r); s = math.sqrt(1 / (df - 1)); return math.tanh(z - 1.96 * s), math.tanh(z + 1.96 * s)

VERSIONS = [  # (version, arm, role, description, pooled file, per-cohort file or None)
    ("RFI primary", "RFI", "primary", "Calibrated ANCOM-BC2 (Amendment 37); Nelore, Charolais mRNA and Alberta all animals (567; Amendment 39).",
     "results/replication_rfi3_calibrated/classification_all.tsv", "results/replication_rfi3_calibrated/estimates_calibrated_all.tsv"),
    ("RFI, Alberta 529", "RFI", "superseded", "As primary, Alberta without the 35 PRJNA448333 ID matches (529); superseded by Amendment 39.",
     "results/replication_rfi3_calibrated/classification_primary.tsv", "results/replication_rfi3_calibrated/estimates_calibrated.tsv"),
    ("RFI uncalibrated", "RFI", "superseded", "ANCOM-BC2 W before permutation calibration, Alberta all animals (script 66; superseded by Amendment 37).",
     "results/replication_rfi3/classification_all.tsv", None),
    ("RFI uncalibrated, Alberta 529", "RFI", "superseded", "ANCOM-BC2 W before calibration, Alberta 529 (script 66).",
     "results/replication_rfi3/classification_primary.tsv", None),
    ("RFI + Wallace RFI NRC", "RFI", "sensitivity", "Primary RFI pooling with Wallace RFI NRC added as a fourth cohort (Amendment 43; intake partly alkane-estimated).",
     "results/replication_rfi4_wallace_sensitivity/classification.tsv", "results/replication_rfi4_wallace_sensitivity/estimates_calibrated.tsv"),
    ("RFI CLR", "RFI", "sensitivity", "Second estimator: CLR linear model, same animals, covariates and targets (Amendment 42); pooled without calibration.",
     "results/clr_pooled/classification_rfi_clr.tsv", "CLR_RFI"),
    ("Methane primary", "Methane", "primary", "Calibrated ANCOM-BC2; Difford and Wallace, CH4 g/d; Wallace M2 is re-analysis and excluded (Amendment 32).",
     "results/methane_pooled/classification_primary.tsv", "results/methane_pooled/estimates_calibrated.tsv"),
    ("Methane with re-analysis", "Methane", "sensitivity", "As primary, Wallace M2 re-analysis included (Amendment 32 item 3).",
     "results/methane_pooled/classification_with_reanalysis.tsv", None),
    ("Methane, Wallace without SE1", "Methane", "sensitivity", "Wallace without the Swedish farm SE1 (possible overlap with Danielsson et al. 2017; Amendment 32 item 2).",
     "results/methane_pooled_noSE1/classification_primary.tsv", "results/methane_pooled_noSE1/estimates_calibrated.tsv"),
    ("Methane, Wallace without SE1, with re-analysis", "Methane", "sensitivity", "As above, Wallace M2 re-analysis included.",
     "results/methane_pooled_noSE1/classification_with_reanalysis.tsv", None),
    ("Methane, Wallace depth-adjusted", "Methane", "sensitivity", "Wallace model with log sequencing depth (Amendment 41).",
     "results/methane_pooled_depth/classification_primary.tsv", "results/methane_pooled_depth/estimates_calibrated.tsv"),
    ("Methane, Wallace depth-adjusted, with re-analysis", "Methane", "sensitivity", "As above, Wallace M2 re-analysis included.",
     "results/methane_pooled_depth/classification_with_reanalysis.tsv", None),
    ("Methane, Wallace per kg DMI", "Methane", "sensitivity", "Wallace CH4 g/kg DMI (intake partly alkane-estimated; Amendment 34 item 2); Difford g/d.",
     "results/methane_pooled_ch4_dmi/classification_primary.tsv", "results/methane_pooled_ch4_dmi/estimates_calibrated.tsv"),
    ("Methane, Wallace per kg DMI, with re-analysis", "Methane", "sensitivity", "As above, Wallace M2 re-analysis included.",
     "results/methane_pooled_ch4_dmi/classification_with_reanalysis.tsv", None),
    ("Methane CLR", "Methane", "sensitivity", "Second estimator: CLR linear model (Amendment 42); M1 is a linear model in both; pooled without calibration.",
     "results/clr_pooled/classification_methane_clr.tsv", "CLR_CH4"),
]
POOLED_COLS = ["set", "id", "target", "rank", "stratum", "published_dir", "k", "r", "ci_lo", "ci_hi", "p", "q_bh", "family_size",
               "tau2", "p_hk", "p_fe", "concordant", "loo_all_sig", "loo_nominal", "label"]
COH_COLS = ["id", "set", "cohort", "table", "rank", "feature", "n", "df", "lfc", "se", "W_raw", "sd_null", "k_null", "W", "p_perm", "r", "lo", "hi", "reanalysis"]

def clr_cohort(arm):
    if arm == "CLR_RFI":
        fs = [f"results/clr/{c}_{r}_clr_targets.tsv" for c in ("nelore", "charolais_mrna") for r in RANKS] + \
             [f"results/clr/alberta_{r}_all_clr_targets.tsv" for r in ("genus", "family")]
    else:
        fs = ["results/clr/difford_bact_clr_targets.tsv", "results/clr/wallace_prokaryotes_clr_targets.tsv"]
    E = pd.concat([rd(f).assign(table=os.path.basename(f).replace("_clr_targets.tsv", "")) for f in fs]); E = E[E.status == "ok"].copy()
    if arm == "CLR_CH4":
        E["cohort"] = E.cohort.str.replace(r"^difford.*", "difford", regex=True).str.replace(r"^wallace.*", "wallace", regex=True)
        E["reanalysis"] = (E.id == "M2") & (E.cohort == "wallace")
    E["W_raw"] = E["W"]; E["r"] = [partial_r(w, d) for w, d in zip(E.W, E.df)]
    lohi = [ci(r, d) for r, d in zip(E.r, E.df)]; E["lo"] = [a for a, b in lohi]; E["hi"] = [b for a, b in lohi]
    return E

def build_file2():
    pooled, coh = [], []
    for v, arm, role, desc, pf, cf in VERSIONS:
        d = rd(pf)
        d = d[~d.set.isin(["Methane", "MethaneSecondary"])] if arm == "RFI" else d[d.set.isin(["Methane", "MethaneSecondary", "Nelore"])]
        if arm == "Methane": d = d[(d.set != "Nelore") | d.id.isin(rd("plan/target_keys.tsv").query("phenotype == 'RME_GC'").id)]
        d = d.reindex(columns=POOLED_COLS); d.insert(0, "role", role); d.insert(0, "arm", arm); d.insert(0, "version", v); pooled.append(d)
        if cf is None: continue
        e = clr_cohort(cf) if cf.startswith("CLR") else rd(cf)
        e = e.reindex(columns=COH_COLS)
        e.insert(0, "calibration", "none (CLR; pooled uncalibrated)" if cf.startswith("CLR") else "W = W_raw / max(1, sd_null); 200 within-design permutations")
        e.insert(0, "role", role); e.insert(0, "arm", arm); e.insert(0, "version", v); coh.append(e)
    PL = pd.concat(pooled, ignore_index=True); CO = pd.concat(coh, ignore_index=True)
    tk = rd("plan/target_keys.tsv").set_index("id")
    CO.insert(CO.columns.get_loc("set") + 1, "target", CO.id.map(tk["target"]))
    CO["set"] = CO["set"].fillna(CO.id.map(tk["set"]))
    PL = PL.rename(columns={"ci_lo": "ci_lo_95", "ci_hi": "ci_hi_95", "q_bh": "q_BH_within_set", "p_hk": "p_Hartung_Knapp", "p_fe": "p_fixed_effect"})
    CO = CO.rename(columns={"lo": "ci_lo_95", "hi": "ci_hi_95", "W": "W_calibrated", "sd_null": "null_SD_of_W", "k_null": "n_permutations", "p_perm": "p_permutation"})
    V = pd.DataFrame([dict(version=v, arm=a, role=ro, description=de, pooled_source=pf,
                           per_cohort_source=("results/clr/*_clr_targets.tsv" if (cf or "").startswith("CLR") else (cf or "(same per-cohort estimates as the version above it)")))
                      for v, a, ro, de, pf, cf in VERSIONS])
    notes = pd.DataFrame({"Item": ["Versions", "Primary", "Per-cohort sheet", "Pooled sheet", "Labels", "Untestable targets", "Built by"],
        "Description": [
        "Every analysis version in results/robustness (script 97) plus the Wallace RFI NRC sensitivity pooling (script 99). See the Versions sheet.",
        "RFI primary and Methane primary (role = primary). All other versions are sensitivity analyses or superseded and set no label.",
        "One row per target, cohort and version. W_raw is the ANCOM-BC2 (or CLR) estimate / SE; null_SD_of_W is its SD over 200 within-design phenotype permutations; W_calibrated = W_raw / max(1, null SD); p_permutation is the share of permutations with |W| at least |W_raw|; r is the partial correlation W / sqrt(W^2 + df) with a Fisher-z 95% CI. Versions sharing per-cohort estimates with another version are not repeated.",
        "One row per target and version. Random-effects (REML) pooling of Fisher z over cohorts; k = cohorts pooled; q is Benjamini-Hochberg within each target set; loo_* are the leave-one-cohort-out checks; label is the pre-registered classification (Section 10, Amendment 14).",
        "No target is labelled replicated in any version.",
        "Targets with k = 0 are listed with label 'untestable' and blank estimates.",
        "scripts/104_additional_files.py, from the files named in the Versions sheet."]})
    f = os.path.join(OD, "Additional_file_2_estimates.xlsx")
    with pd.ExcelWriter(f, engine="openpyxl") as w:
        notes.to_excel(w, "Notes", index=False); V.to_excel(w, "Versions", index=False)
        CO.to_excel(w, "Per-cohort estimates", index=False); PL.to_excel(w, "Pooled results", index=False)
    style(f, {"Per-cohort estimates": 4, "Pooled results": 4}); return f, CO, PL

# ---------------------------------------------------------------- formatting
def style(f, freeze):
    wb = load_workbook(f); hdr = Font(name="Arial", bold=True, size=10); body = Font(name="Arial", size=10)
    fill = PatternFill("solid", fgColor="E8E6E0")
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row:
                c.font = body
                if isinstance(c.value, float): c.number_format = "0.0000"
        for c in ws[1]: c.font = hdr; c.fill = fill; c.alignment = Alignment(wrap_text=True, vertical="top")
        for i, col in enumerate(ws.columns, 1):
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
            ws.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 60 if ws.title != "Notes" else 140)
        if ws.title == "Notes":
            for c in ws["B"]: c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = ws.cell(row=2, column=freeze.get(ws.title, 1))
        ws.auto_filter.ref = ws.dimensions
    wb.save(f)

if __name__ == "__main__":
    f1, T = build_file1(); f2, CO, PL = build_file2()
    print("wrote", f1, "| targets", len(T)); print("wrote", f2, "| per-cohort rows", len(CO), "| pooled rows", len(PL))
    print(PL.groupby(["version"]).label.value_counts().unstack(fill_value=0).to_string())
