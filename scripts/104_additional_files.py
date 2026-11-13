#!/usr/bin/env python3
"""104: Additional files 1 and 2 for the manuscript, built from existing outputs only (no fits).

Additional file 1: every pre-specified target with its source, published direction, mapping, and per cohort the
mapped feature, prevalence and fit status. Prevalence is not stored in any results table: scripts 21, 65 and 81
print it to their run logs, so it is parsed from those logs and checked against each fit's status.
Additional file 2: per-cohort calibrated estimates (W / max(1, null SD), permutation p) and pooled results for every
analysis version found by scripts 97 and 99, with the primary versions marked, and the Amendment 46 effect sizes, coverage and Freedman-Lane
results (late analysis, run after all primary results were seen). Every permutation p is (b + 1)/(B + 1) (Amendment 46.5.1).

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
    # Amendment 46 (review items A1, A4; script 118): evidence ledger and overlap audit, added when their outputs exist
    extra = []
    for nm, fn in (("Evidence ledger", "evidence_ledger.tsv"), ("Overlap audit", "overlap_audit.tsv")):
        fp = P("results", "late_a46", "ledger", fn)
        if os.path.exists(fp): extra.append((nm, pd.read_csv(fp, sep="\t", dtype={"id": str})))
    if extra:
        notes = pd.concat([notes, pd.DataFrame({"Item": ["Evidence ledger", "Overlap audit"], "Description": [
            "One row per testable target (27 RFI, 11 methane): source claim, site, phenotype, published direction, evidence category (matched claim; published claim without direction; candidate taxon; cross-site transfer, fecal claim tested in rumen; cross-phenotype transfer, residual methane claim tested against daily output), independent cohorts k, animals, and the Amendment 46 CLR partial r with its HKSJ (k >= 2) or cohort (k = 1) interval. Late analysis, run after all primary results were seen (Amendment 46; script 118).",
            "One row per target and cohort: independent test, discovery cohort (not tested), re-analysis (sensitivity pooling only) or possible overlap, with the amendment that sets the rule and whether a primary estimate is affected. Possible overlaps: up to 48 Li et al. discovery animals among the 567 Alberta animals (no Li target has a tested feature in Alberta, so no estimate is affected); Wallace farm SE1 and the M1 source Danielsson et al. 2017 (sensitivity without SE1)."]})], ignore_index=True)
    f = os.path.join(OD, "Additional_file_1_targets.xlsx")
    write_book(f, [("Notes", notes), ("Targets", T)] + extra, {"Targets": 3, "Evidence ledger": 3, "Overlap audit": 3}); return f, T

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
    ("Methane CLR", "Methane", "sensitivity", "Second estimator: CLR linear model (Amendment 42); M1 is a linear model in both. Pooled file (script 96) uses uncalibrated CLR W; per-cohort rows here are calibrated (addendum; null SD 0.95 to 1.08), which leaves every label unchanged.",
     "results/clr_pooled/classification_methane_clr.tsv", "CLR_CH4"),
]
# Amendment 44: late analyses, run after all primary results were seen (pre-specified in plan Section 6, not run at the time). Added when their outputs exist.
def _kperm(m):
    f = P(f"results/late_{m}/pooled/n_permutations.txt"); return open(f).read().strip() if os.path.exists(f) else "?"
LATE = "Late analysis, run after all primary results were seen (Amendment 44; pre-specified in plan Section 6, not run at the time). "
for m, nm in (("linda", "LinDA (MicrobiomeStat, feature.dat.type = count, other settings default)"), ("aldex2", "ALDEx2 1.34.0 (aldex.clr, 128 Monte Carlo instances, denom all; aldex.glm)")):
    if os.path.exists(P(f"results/late_{m}/pooled/classification_rfi.tsv")):
        k = _kperm(m); cal = f"same tables, animals, covariates and targets as the primary; calibrated by {k} within-design permutations as in Amendment 37"
        VERSIONS += [(f"RFI {nm.split(' ')[0]} (late)", "RFI", "late", LATE + nm + "; " + cal + "; Alberta all animals.",
                      f"results/late_{m}/pooled/classification_rfi.tsv", f"results/late_{m}/pooled/estimates_calibrated_rfi.tsv"),
                     (f"Methane {nm.split(' ')[0]} (late)", "Methane", "late", LATE + nm + "; " + cal + "; M1 is the primary calibrated linear model; Wallace M2 excluded as re-analysis.",
                      f"results/late_{m}/pooled/classification_methane_primary.tsv", f"results/late_{m}/pooled/estimates_calibrated_methane.tsv")]
POOLED_COLS = ["set", "id", "target", "rank", "stratum", "published_dir", "k", "r", "ci_lo", "ci_hi", "p", "q_bh", "family_size",
               "tau2", "p_hk", "p_fe", "concordant", "loo_all_sig", "loo_nominal", "label"]
COH_COLS = ["id", "set", "cohort", "table", "rank", "feature", "n", "df", "lfc", "se", "W_raw", "sd_null", "k_null", "W", "p_perm", "r", "lo", "hi", "reanalysis"]

def clr_cohort(arm):
    if arm == "CLR_RFI":
        fs = [f"results/clr/{c}_{r}_clr_targets.tsv" for c in ("nelore", "charolais_mrna") for r in RANKS] + \
             [f"results/clr/alberta_{r}_all_clr_targets.tsv" for r in ("genus", "family")]
    else:
        fs = ["results/clr/difford_bact_clr_targets.tsv", "results/clr/wallace_prokaryotes_clr_targets.tsv"]
    E = pd.concat([rd(f).assign(table=os.path.basename(f).replace("_clr_targets.tsv", "")) for f in fs], ignore_index=True); E = E[E.status == "ok"].copy()
    E["W_raw"] = E["W"]
    if arm == "CLR_CH4":
        # Amendment 42 addendum: CLR methane W calibrated by its 200 within-herd or within-farm permutations,
        # W = W_raw / max(1, null SD) as in script 83; p = (b + 1) / (n + 1), the rule behind the CLR p values in the Amendment 42 result.
        cal = []
        for tab, g in E.groupby("table"):
            fs = [f"results/clr/{tab}_seed{s}_clr_targets.tsv" for s in range(1001, 1201)]
            fs = [f for f in fs if os.path.exists(P(f))]; assert len(fs) == 200, (tab, len(fs))
            N = pd.concat([rd(f) for f in fs]); N = N[N.status == "ok"]
            for i, row in g.iterrows():
                w0 = N.W[N.id == row.id].to_numpy(); assert len(w0) >= 190, (tab, row.id, len(w0))
                cal.append((i, w0.std(ddof=1), len(w0), float(((abs(w0) >= abs(row.W_raw)).sum() + 1) / (len(w0) + 1))))
        for i, sd, k, p in cal:
            E.loc[i, "sd_null"] = sd; E.loc[i, "k_null"] = k; E.loc[i, "p_perm"] = p
            E.loc[i, "W"] = E.loc[i, "W_raw"] / max(1.0, sd)
        E["cohort"] = E.cohort.str.replace(r"^difford.*", "difford", regex=True).str.replace(r"^wallace.*", "wallace", regex=True)
        E["reanalysis"] = (E.id == "M2") & (E.cohort == "wallace")
    E["r"] = [partial_r(w, d) for w, d in zip(E.W, E.df)]
    lohi = [ci(r, d) for r, d in zip(E.r, E.df)]; E["lo"] = [a for a, b in lohi]; E["hi"] = [b for a, b in lohi]
    return E

def replicated_note(PL):
    x = PL[PL.label == "replicated"]
    if not len(x): return "No target is labelled replicated in any version."
    return "Targets labelled replicated: " + "; ".join(f"{i} in {v}" for i, v in zip(x.id, x.version)) + ". No other target is labelled replicated in any version."

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
        lm = re.match(r"results/late_(linda|aldex2)/", cf or "")
        e.insert(0, "calibration", {"CLR_RFI": "none (CLR, RFI arm: no CLR permutations, Amendment 42)",
                                    "CLR_CH4": "W = W_raw / max(1, sd_null); 200 CLR permutations (Amendment 42 addendum)"}.get(cf,
                                    f"W = W_raw / max(1, sd_null); {_kperm(lm.group(1))} within-design permutations (late, Amendment 44)" if lm else
                                    "W = W_raw / max(1, sd_null); 200 within-design permutations"))
        e.insert(0, "role", role); e.insert(0, "arm", arm); e.insert(0, "version", v); coh.append(e)
    PL = pd.concat(pooled, ignore_index=True); CO = pd.concat(coh, ignore_index=True)
    # Amendment 46.5.1: every permutation p as (b + 1)/(B + 1). Versions whose p was b/B (scripts 89, 83, 106) are converted from b = p * B,
    # which must be a whole number; script 99 (RFI + Wallace RFI NRC) and the calibrated CLR methane rows already use (b + 1)/(B + 1).
    conv = CO.p_perm.notna() & ~CO.version.isin(B1_VERSIONS)
    b = CO.loc[conv, "p_perm"] * CO.loc[conv, "k_null"]; assert (b - b.round()).abs().max() < 1e-6, "p_perm * n_permutations is not a whole number"
    n0 = int((CO.p_perm == 0).sum()); CO.loc[conv, "p_perm"] = (b.round() + 1) / (CO.loc[conv, "k_null"] + 1)
    assert not (CO.p_perm == 0).any(); print(f"46.5.1: {int(conv.sum())} permutation p converted to (b + 1)/(B + 1); {n0} were 0, now none")
    # Amendment 46.5.6: published direction (and stratum) filled in for the CLR rows, from plan/target_keys.tsv and the primary of the same arm
    tkd = rd("plan/target_keys.tsv").set_index("id")["direction"]
    pdir = lambda arm, i: (-1 if re.search(r"higher in L-RFI|lower with higher", str(tkd.get(i, ""))) else 1 if re.search(r"higher in H-RFI|higher with higher", str(tkd.get(i, ""))) else float("nan"))
    clr = PL.version.isin(["RFI CLR", "Methane CLR"])
    PL.loc[clr, "published_dir"] = [pdir(a, i) for a, i in zip(PL.arm[clr], PL.id[clr])]
    for arm, prim in (("RFI", "RFI primary"), ("Methane", "Methane primary")):
        smap = PL[PL.version == prim].set_index("id")["stratum"]; m = clr & (PL.arm == arm); PL.loc[m, "stratum"] = PL.id[m].map(smap)
    tk = rd("plan/target_keys.tsv").set_index("id")
    CO.insert(CO.columns.get_loc("set") + 1, "target", CO.id.map(tk["target"]))
    CO["set"] = CO["set"].fillna(CO.id.map(tk["set"]))
    PL = PL.rename(columns={"ci_lo": "ci_lo_95", "ci_hi": "ci_hi_95", "q_bh": "q_BH_within_set", "p_hk": "p_Hartung_Knapp", "p_fe": "p_fixed_effect"})
    CO = CO.rename(columns={"lo": "ci_lo_95", "hi": "ci_hi_95", "W": "W_calibrated", "sd_null": "null_SD_of_W", "k_null": "n_permutations", "p_perm": "p_permutation"})
    V = pd.DataFrame([dict(version=v, arm=a, role=ro, description=de, pooled_source=pf,
                           per_cohort_source=("results/clr/*_clr_targets.tsv" if (cf or "").startswith("CLR") else (cf or "(same per-cohort estimates as the version above it)")))
                      for v, a, ro, de, pf, cf in VERSIONS])
    notes = pd.DataFrame({"Item": ["Versions", "Primary", "Per-cohort sheet", "Pooled sheet", "Labels", "Untestable targets",
                                   "A46 sheets", "A46 effect size", "A46 coverage", "A46 Freedman-Lane", "Built by"],
        "Description": [
        "Every analysis version in results/robustness (script 97), the Wallace RFI NRC sensitivity pooling (script 99) and, with role = late, the LinDA and ALDEx2 analyses run after all primary results were seen (Amendment 44, script 106). See the Versions sheet.",
        "RFI primary and Methane primary (role = primary). All other versions are sensitivity analyses, superseded or late analyses and set no label.",
        "One row per target, cohort and version. W_raw is the ANCOM-BC2 (or CLR) estimate / SE; null_SD_of_W is its SD over 200 within-design phenotype permutations; W_calibrated = W_raw / max(1, null SD); p_permutation is (b + 1) / (n + 1), b = permutations with |W| at least |W_raw|, n = permutations (Amendment 46.5.1; earlier versions of this file gave b / n for most versions, so 17 values were 0); r is the partial correlation W / sqrt(W^2 + df) with a Fisher-z 95% CI. Versions sharing per-cohort estimates with another version are not repeated.",
        "One row per target and version. Random-effects (REML) pooling of Fisher z over cohorts; k = cohorts pooled; q is Benjamini-Hochberg within each target set; loo_* are the leave-one-cohort-out checks; label is the pre-specified classification (Section 10, Amendment 14; the plan was not registered with a third party). published_dir: +1 higher with higher trait (or in H-RFI), -1 lower (or higher in L-RFI), blank where no direction was published; filled in for the CLR rows from plan/target_keys.tsv.",
        replicated_note(PL),
        "Targets with k = 0 are listed with label 'untestable' and blank estimates.",
        "Sheets starting 'A46' are late analyses, run after all primary results were seen (Amendment 46, in response to an external statistical review); they set no label. Source files in results/late_a46/.",
        "Partial correlation of the target's CLR abundance with the trait given the cohort's covariates, by OLS (r = t / sqrt(t^2 + df), Fisher z with SE 1 / sqrt(df - 1)); for Charolais a point-biserial partial correlation (extreme-group design). Pooled by REML with Hartung-Knapp-Sidik-Jonkman intervals (t on k - 1 df), the cohort's own interval where k = 1; fixed-effect and REML-Wald intervals beside it. Benchmarks r = 0.30 and 0.15 in the published direction (both bounds where none). The exclusion statements in the manuscript use this estimand and interval. Per-cohort p_fl: Freedman-Lane permutation, 9,999 permutations, (b + 1)/(B + 1).",
        "Coverage of the intervals by simulation of the CLR outcome with each cohort's real trait and covariates (1,000 replicates per table and rho; rho 0 to 0.30; tau 0 to 0.10). Pass rules: (a) per-cohort coverage averaged over targets 0.93 to 0.97; (b) pooled HKSJ coverage averaged over targets with k >= 2 at least 0.93. Both pass (results/late_a46/coverage_fix1/; the first run in results/late_a46/coverage/ had a sign error and is kept unchanged).",
        "ANCOM-BC2 on the methane primary targets with 200 Freedman-Lane permutations (trait residualised on the covariates, residuals permuted within herd or farm), p = (b + 1)/(B + 1). The 'A46 Wallace farms' sheets report the leave-one-farm-out genera (Amendment 40) with per-farm estimates and 95% intervals (uncalibrated within-farm ANCOM-BC2 fits), farms fitted, sign agreement and I^2 (script 117).",
        "scripts/104_additional_files.py, from the files named in the Versions sheet."]})
    f = os.path.join(OD, "Additional_file_2_estimates.xlsx")
    A = a46_sheets(); extra = [(k, v) for k, v in A.items() if v is not None]
    write_book(f, [("Notes", notes), ("Versions", V), ("Per-cohort estimates", CO), ("Pooled results", PL)] + extra,
               {"Per-cohort estimates": 4, "Pooled results": 4, "A46 effect size, per cohort": 3, "A46 effect size, pooled": 3}); return f, CO, PL

# ---------------------------------------------------------------- Amendment 46 (late analysis, run after all primary results were seen)
B1_VERSIONS = {"RFI + Wallace RFI NRC", "Methane CLR"}
A46 = "results/late_a46"
def a46_sheets():
    """Amendment 46 outputs as sheets; the coverage sheet uses coverage_fix1 (the first run in coverage/ had a sign bug, Amendment 46 result item 2)."""
    out = {}
    def opt(f): return rd(f) if os.path.exists(P(f)) else None
    out["A46 effect size, per cohort"] = opt(f"{A46}/effect_size/cohort_estimates.tsv")
    pl = opt(f"{A46}/effect_size/pooled.tsv")
    if pl is not None:
        pl = pl.rename(columns={"lo": "lo_95", "hi": "hi_95", "p": "p_pooled", "interval": "interval_k2plus_or_cohort", "bound": "largest_r_in_published_direction_within_interval",
                                "excl_r30": "excludes_r030_benchmark", "excl_r15": "excludes_r015_benchmark", "excl_any": "interval_wholly_against_published_direction",
                                "r_ancom": "r_ANCOMBC2_calibrated", "lo_ancom": "lo_ANCOMBC2", "hi_ancom": "hi_ANCOMBC2", "excl_r30_ancom": "excludes_r030_ANCOMBC2_derived",
                                "excl_r15_ancom": "excludes_r015_ANCOMBC2_derived", "label_a37": "label_primary_Amendment37"})
    out["A46 effect size, pooled"] = pl
    out["A46 exclusion counts"] = opt(f"{A46}/effect_size/exclusion_counts.tsv")
    cs, ps = opt(f"{A46}/coverage_fix1/cohort_summary.tsv"), opt(f"{A46}/coverage_fix1/pooled_summary.tsv")
    if cs is not None and ps is not None:
        cs.insert(0, "level", "per cohort (rule a)"); ps.insert(0, "level", "pooled HKSJ, targets with k >= 2 (rule b)")
        out["A46 coverage"] = pd.concat([cs, ps], ignore_index=True)
    out["A46 Freedman-Lane ANCOM-BC2"] = opt(f"{A46}/effect_size/ancombc2_fl_methane.tsv")
    out["A46 Wallace farms, summary"] = opt(f"{A46}/wallace_farms/summary.tsv")          # script 117 (46.5.4)
    out["A46 Wallace farms, per farm"] = opt(f"{A46}/wallace_farms/per_farm_intervals.tsv")
    return out

# ---------------------------------------------------------------- formatting
def write_book(f, sheets, freeze):
    """Write to a temporary file and move it into place only when complete, so a failure never leaves a broken workbook."""
    tmp = f + ".tmp.xlsx"
    with pd.ExcelWriter(tmp, engine="openpyxl") as w:
        for name, df in sheets: df.to_excel(w, sheet_name=name, index=False)
    style(tmp, freeze); os.replace(tmp, f)

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
