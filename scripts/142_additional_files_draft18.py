#!/usr/bin/env python3
"""142: Additional files 1 to 4 for draft 18 (Amendment 48.10 and 48.11; reporting only, no fit). Writes results/additional_files_draft18/
and refuses to overwrite; the draft 17 files in results/additional_files/ are left as they are.

Additional file 1: the draft 17 workbook with
  - 'Evidence ledger' rebuilt from the final Amendment 47 intervals and labels with the 48.8 mapping-confidence column and the 48.11 source
    audit (results/late_a48/ledger/evidence_ledger_a47_v2.tsv, scripts 141 and 144: S1 a cross-phenotype and cross-host transfer; M1, M2 and S1
    sources from the audit; "met the Amendment 47 coverage rule" in place of "validated"); 'Source audit' (results/late_a48/source_audit.tsv);
    the 'Targets' source and published-trait columns of M1, M2 and S1 from the audit;
  - the Amendment 46 ledger kept as 'Ledger A46 (historical)';
  - a 'Mapping confidence' sheet (48.8, one row per target and cohort) and a 'Targets' column giving each cohort's class.
Additional file 2: the draft 17 workbook plus one sheet per Amendment 48 output (prefix 'A48'); the clade-ratio grid carries the analysis
  decoded from its seed (Amendment 48 result, item 1). The Figure 4 row that pooled Difford g per day with Wallace per kg intake stays in the
  existing 'Pooled results' sheet, version 'Methane, Wallace per kg DMI'.
Additional file 3: docs/ANALYSIS_CLASSES.md (48.10), docs/DEVIATIONS.md and plan/ANALYSIS_PLAN.md, as 116 does, styled from the draft 17 file.
Additional file 4: the 48.1 forest plot (results/figures_draft18/af4_seven_claims.pdf, script 139).
Draft 20 (reporting only, no new fit): the 'A48 seven claims, per cohort' sheet carries the parametric p of the rescaled W (used by the labels), the
  empirical permutation p (b + 1)/(B + 1) and their agreement at 0.05, from results/late_p_compare/seven_claims_p.tsv (script 146); the stored b/B
  column is renamed 'p_perm_ancombc2_b_over_B (as stored)'. Run 146 first.
Usage: 142_additional_files_draft18.py [root]   (needs pandoc >= 2.11, openpyxl; imports write_book from scripts/104_additional_files.py)
"""
import os, re, sys, shutil, subprocess, tempfile, zipfile, importlib.util
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
os.environ["RUMEN_ROOT"] = ROOT
spec = importlib.util.spec_from_file_location("af104", P("scripts", "104_additional_files.py")); af104 = importlib.util.module_from_spec(spec); spec.loader.exec_module(af104)
OD = P("results", "additional_files_draft18"); D17 = P("results", "additional_files")
OUT = {k: os.path.join(OD, f) for k, f in {1: "Additional_file_1_targets.xlsx", 2: "Additional_file_2_estimates.xlsx",
       3: "Additional_file_3_analysis_plan.docx", 4: "Additional_file_4_seven_claims.pdf"}.items()}
ex = [f for f in OUT.values() if os.path.exists(f)]
if ex: sys.exit(f"output exists, never overwritten: {ex}")
os.makedirs(OD, exist_ok=True)
rd = lambda f: pd.read_csv(P(f), sep="\t", dtype={"id": str})
A48 = "results/late_a48"
COH = {"nelore": "Nelore", "charolais_mrna": "Charolais", "alberta_all": "Alberta", "difford": "Difford", "wallace": "Wallace"}

# ---------------------------------------------------------------- Additional file 1
b1 = pd.read_excel(os.path.join(D17, "Additional_file_1_targets.xlsx"), sheet_name=None, dtype={"ID": str, "id": str})
assert list(b1) == ["Notes", "Targets", "Evidence ledger", "Overlap audit"], list(b1)
led = rd(f"{A48}/ledger/evidence_ledger_a47_v2.tsv"); mc = rd(f"{A48}/mapping/mapping_classes.tsv"); aud = rd(f"{A48}/source_audit.tsv")
assert len(led) == 38 and len(mc) == 61
T = b1["Targets"].copy()
cls = mc.assign(c=mc.cohort.map(COH)).groupby("id").apply(lambda x: "; ".join(f"{c} {k}" for c, k in zip(x.c, x["class"])), include_groups=False)
T.insert(T.columns.get_loc("Target") + 1, "Mapping confidence (48.8)", T["ID"].map(cls).fillna(""))
for i in ("M1", "M2", "S1"):   # 48.11 source audit: the ledger's source and phenotype for these three
    m = T["ID"] == i; assert m.sum() == 1, i; r = led[led.id == i].iloc[0]
    T.loc[m, "Source"] = r.source; T.loc[m, "Published trait"] = r.source_phenotype
M = mc.assign(cohort=mc.cohort.map(COH).fillna(mc.cohort))[["id", "target", "arm", "cohort", "rank", "class", "basis", "charolais_mrna_vs_16s_spearman", "analysis"]]
n1 = b1["Notes"].copy()
n1.loc[n1.Item == "Evidence ledger", "Description"] = ("One row per testable target (27 RFI, 11 methane), rebuilt for draft 18 from the final Amendment 47 results (script 141): "
    "source claim, site, phenotype, published direction and evidence category as in Amendment 46; the CLR partial r with the modified Knapp-Hartung interval "
    "(k >= 2) or the cohort's own interval (k = 1), the method, the role (pooled or one independent cohort) and the validation statement ('descriptive' "
    "intervals carry no exclusion statement); exclusion of 0, r = 0.15 and r = 0.30 in the published direction; the Amendment 47 CLR label beside the "
    "Amendment 37 primary label; the mapping-confidence class of each cohort (48.8); and, for M1, M2 and S1, the source, phenotype and contrast from "
    "the 48.11 source audit (S1 Sharpea is a cross-phenotype and cross-host transfer: a methane-yield claim in sheep tested against daily output in cows; "
    "M1 and M2 keep the matched class, flagged because the source contrast was a group or cluster comparison). 'met the Amendment 47 coverage rule' marks "
    "pooled intervals whose procedure met the coverage rule of Amendment 47. Late analysis, run after all primary results were seen.")
n1 = pd.concat([n1, pd.DataFrame({"Item": ["Ledger A46 (historical)", "Mapping confidence", "Mapping confidence (Targets)", "Source audit"], "Description": [
    "The Amendment 46 ledger (script 118) as published with draft 17, kept for the record. Its intervals (HKSJ) are superseded by the Amendment 47 "
    "modified Knapp-Hartung intervals in 'Evidence ledger'.",
    "Amendment 48.8: one row per target and independent cohort. MC1, one key with the exact genus name in that cohort's taxonomy; MC2, a GTDB-suffixed "
    "genus or several keys summed; MC3, family or phylum aggregate (and the clade ratio); MC4, species level. Assigned by rule from plan/target_keys.tsv "
    "and results/late_mapping/target_mapping_by_cohort.tsv (script 135). With MC1 estimates only, 5 of the 8 RFI r = 0.30 exclusions and 4 of the 6 "
    "methane r = 0.15 exclusions hold (Amendment 48 result, item 9, as corrected).",
    "The 'Targets' sheet gains the column 'Mapping confidence (48.8)', the class per independent cohort; blank for untestable targets.",
    "Amendment 48.11 source audit of the methane claims (results/late_a48/source_audit.tsv): per source, host and n, contrast, phenotype and units, "
    "direction and class. Rule: a methane target is a matched claim only if a source independent of the test cohort reports it in cattle for daily "
    "output with a stated direction. Only S1 changes class; labels, estimates and intervals are unchanged."]})], ignore_index=True)
S1 = ["Notes", "Targets", "Evidence ledger", "Ledger A46 (historical)", "Overlap audit", "Mapping confidence", "Source audit"]; assert all(len(n) <= 31 for n in S1)
af104.write_book(OUT[1], [("Notes", n1), ("Targets", T), ("Evidence ledger", led), ("Ledger A46 (historical)", b1["Evidence ledger"]),
                          ("Overlap audit", b1["Overlap audit"]), ("Mapping confidence", M), ("Source audit", aud)],
                 {"Targets": 4, "Evidence ledger": 3, "Ledger A46 (historical)": 3, "Overlap audit": 3, "Mapping confidence": 2, "Source audit": 3})
print("wrote", OUT[1])

# ---------------------------------------------------------------- Additional file 2
b2 = pd.read_excel(os.path.join(D17, "Additional_file_2_estimates.xlsx"), sheet_name=None, dtype={"id": str})
assert "Notes" in b2 and "Pooled results" in b2 and not any(k.startswith("A48") for k in b2)
g = rd(f"{A48}/clade_ratio/grid.tsv"); g["analysis_no"] = (g.seed - 4890000) // 100
assert ((g.seed - 4890000) % 100 == g.cell).all()
g.insert(1, "clade_ratio_analysis", g.analysis_no.map({1: "Wallace CH4 per kg DMI", 2: "Wallace CH4 g/d", 3: "Difford CH4 g/d", 4: "Wallace CH4 per kg DMI without SE1"}))
# draft 19 (wording; no number changes): descriptive-only labels, hypothetical attenuation scenarios, source-verified dairy table
DESC = "descriptive only; estimates not independent"
ag = rd(f"{A48}/dairy/agreement.tsv")
ag["bootstrap"] = ag.bootstrap.str.replace(r"^(\d+) genus resamples, seed (\d+)$", r"genus resampling (\1 resamples, seed \2); descriptive only; does not preserve dependence through shared animals and CLR reference", regex=True)
assert ag.bootstrap.str.contains("does not preserve dependence").all()
ag = ag.rename(columns={"rho_lo": "rho_lo (genus resampling; descriptive only)", "rho_hi": "rho_hi (genus resampling; descriptive only)", "p_sign_binom": f"p_sign_binom ({DESC})"})
cc = rd(f"{A48}/seven/concordance.tsv").rename(columns={"p_binom": f"p_binom ({DESC})"})
at = rd(f"{A48}/power/attenuation.tsv")
at["note"] = ("hypothetical reliability scenarios: true r = observed r / sqrt(trait x feature reliability), capped at 1; reliabilities 0.5, 0.7, 0.9 are "
              "assumed values, not estimates; the Charolais GTDB vs SILVA correlation (Amendment 45) is not a measure of feature reliability")
pe = rd(f"{A48}/power/published_effects.tsv")
pe.insert(len(pe.columns) - 1, "use", "not on the CLR scale (Li et al.: extreme-group r from a p value; Nelore: coefficient/SE with an assumed df); not compared with the CLR intervals in the manuscript")
nh = rd(f"{A48}/dairy/not_harmonized_sources.tsv")
# draft 20 (reporting only): parametric p of the rescaled W beside the empirical permutation p (script 146)
pcf = P("results", "late_p_compare", "seven_claims_p.tsv")
if not os.path.exists(pcf): sys.exit(f"missing {pcf}: run scripts/146_p_compare.py first")
s7 = rd(f"{A48}/seven/seven_claims.tsv"); pc = pd.read_csv(pcf, sep="\t", dtype={"id": str})
s7 = s7.merge(pc[["id", "cohort", "B", "p_param_ancombc2", "p_perm_emp_ancombc2", "agree_at_0.05"]], on=["id", "cohort"], how="left", validate="1:1")
assert len(s7) == 14 and s7.p_param_ancombc2.notna().all() and (s7.p_perm_ancombc2 * s7.B).sub((s7.p_perm_ancombc2 * s7.B).round()).abs().max() < 1e-9
s7 = s7.rename(columns={"p_perm_ancombc2": "p_perm_ancombc2_b_over_B (as stored)", "B": "B_permutations"})
k = list(s7.columns).index("p_perm_ancombc2_b_over_B (as stored)"); new = ["B_permutations", "p_perm_emp_ancombc2", "p_param_ancombc2", "agree_at_0.05"]
s7 = s7[[c for c in s7.columns[:k + 1] if c not in new] + new + [c for c in s7.columns[k + 1:] if c not in new]]
assert int(s7["agree_at_0.05"].sum()) == 13
P48 = ("48.1: p_param_ancombc2 is the normal-theory p of the W rescaled by its permutation-null SD (used by the labels); it corrects the spread of the null "
       "but not its shape. p_perm_emp_ancombc2 is the empirical permutation p, (b + 1)/(B + 1). Draft 20 reports their agreement for the seven claims (13 of 14 at 0.05).")
S = [  # (sheet name <= 31 chars, frame, note)
 ("A48 seven claims, per cohort", s7, "48.1: per cohort n, df, CLR partial r with 95% Fisher interval, parametric and Freedman-Lane p; ANCOM-BC2 W rescaled by its permutation-null SD when that exceeds 1 (W_ancombc2_cal), r and interval; the stored permutation share b/B, the empirical permutation p (b + 1)/(B + 1), the normal-theory p of the rescaled W and their agreement at 0.05 (draft 20, script 146; results/late_p_compare/); minimum detectable r at 80% power."),
 ("A48 seven claims, pooled", rd(f"{A48}/seven/seven_pooled.tsv"), "48.1: mKH pooled CLR r and interval, tau-hat, I-squared, Q p; fixed-effect comparator; Amendment 37 pooled r, q and label; Amendment 47 label; reasons for no replicated label; pattern codes A to E per estimator; secondary criteria S-a and S-b (reported only, never replication)."),
 ("A48 seven claims, secondary", rd(f"{A48}/seven/secondary_matched.tsv"), "48.1: secondary criteria for every matched claim (k = 1 and k >= 2)."),
 ("A48 concordance", cc, "48.1: share of cohort estimates with the published sign. The binomial p is descriptive only: the estimates are not independent across targets (draft 19 reports the sign counts only)."),
 ("A48 counts", rd(f"{A48}/seven/counts_by_classification.tsv"), "48.1: matched claims and secondary-criterion counts with S1 as a matched claim and as a phenotype-transfer test."),
 ("A48 S1 against yield", rd(f"{A48}/seven/s1_yield.tsv"), "48.1: S1 Sharpea against Wallace CH4 per kg DMI (k = 1)."),
 ("A48 ANCOM-BC2 configuration", rd(f"{A48}/ancombc2/config.tsv"), "48.2(a): every ancombc2() argument as executed, defaults resolved from the installed ANCOMBC version."),
 ("A48 ANCOM-BC2 label p check", rd(f"{A48}/ancombc2/label_p_check.tsv"), "48.2(a): labels use the normal-theory p of the rescaled W (W divided by its permutation-null SD when that exceeds 1), not the 200-permutation p."),
 ("A48 ANCOM-BC2 filtered null", rd(f"{A48}/ancombc2/null_summary.tsv"), "48.2(b): 200 saved within-design permutation fits per table; share with p < 0.05 unfiltered and with passed_ss; BH and diff_robust counts; P(any) with Clopper-Pearson 95% intervals."),
 ("A48 ANCOM-BC2 null, per shuffle", rd(f"{A48}/ancombc2/null_per_shuffle.tsv"), "48.2(b): the per-shuffle values behind the summary."),
 ("A48 ANCOM-BC2 real fits", rd(f"{A48}/ancombc2/real_fits.tsv"), "48.2(b): the same quantities in the real fits."),
 ("A48 ANCOM-BC2 positive controls", rd(f"{A48}/ancombc2/positive_controls.tsv"), "48.2(b): BH and diff_robust counts for the positive-control terms."),
 ("A48 Charolais structural zeros", rd(f"{A48}/ancombc2/charolais_struc_zero.tsv"), "48.2(c): Charolais refit with group = RFI_class, struc_zero = TRUE; sensitivity only, sets nothing."),
 ("A48 ANCOM-BC2 reading", rd(f"{A48}/ancombc2/reading.tsv"), "48.2(d): the pre-fixed reading."),
 ("A48 dairy genera", rd(f"{A48}/dairy/genera.tsv"), "48.3: harmonized shared genera, CLR partial r with 95% interval per cohort, Freedman-Lane p (9,999 within herd or farm), Q and I-squared; CH4, depth-adjusted and days-in-milk variants."),
 ("A48 dairy agreement", ag, "48.3: Spearman rho (descriptive), Pearson r of Fisher z, sign agreement, I-squared; per variant and subset. The rho interval comes from resampling genera and the sign p from a binomial test: both are descriptive only, because genus estimates share animals and the CLR reference (draft 19 reports rho and the sign counts only)."),
 ("A48 dairy not harmonized", nh, "48.3: differences between the cohorts that cannot be harmonized (breed, production system and diet, methane method, rumen sampling, primers and read length), from the Materials and Methods of Difford et al. 2018 and Wallace et al. 2019 (results/late_a48/dairy/not_harmonized_sources.tsv, script 144); tabulated, not modelled."),
 ("A48 farm random effects", rd(f"{A48}/farms/random_effects.tsv"), "48.4: REML random effects over Wallace farms per feature and estimand, mKH interval, 95% prediction interval, leave-one-farm-out."),
 ("A48 farm leave-one-out", rd(f"{A48}/farms/leave_one_farm_out.tsv"), "48.4: metafor leave1out per feature and estimand."),
 ("A48 farms, per farm CLR", rd(f"{A48}/farms/per_farm_clr.tsv"), "48.4: per-farm CLR partial r behind the random effects."),
 ("A48 CLR grid, robustness", rd(f"{A48}/clr_grid/robustness.tsv"), "48.5: per target, whether its reference statement holds in all 12 cells (zero handling x reference) and the cells where it fails. The statement is the strongest one in the reference cell: r = 0.15 for N04, N11, N13, N20, N22 and the methane targets, r = 0.30 for C5, C6 and N02. Every Charolais r = 0.30 exclusion holds in all 12 cells (results/manifest/draft19_numbers.tsv); the failures listed concern r = 0.15 or zero."),
 ("A48 CLR grid, changes", rd(f"{A48}/clr_grid/changes.tsv"), "48.5: every change from the reference cell (pseudo-count 0.5, R0)."),
 ("A48 CLR grid, pooled", rd(f"{A48}/clr_grid/pooled_cells.tsv"), "48.5: pooled estimates per cell."),
 ("A48 CLR grid, per cohort", rd(f"{A48}/clr_grid/cohort_cells.tsv"), "48.5: cohort estimates per cell (the Charolais intervals excluding 0 are read from here)."),
 ("A48 Charolais influence", rd(f"{A48}/charolais/summary.tsv"), "48.6: 30 leave-one-animal-out refits of the whole CLR pipeline; refits in which the reference statement (r = 0.15 except N02, r = 0.30) holds; DFBETAS; the 48.5 grid beside it. All six Charolais exclusions hold at r = 0.30 in 30 of 30 refits (results/manifest/draft19_numbers.tsv); the failures listed concern r = 0.15 or zero."),
 ("A48 Charolais refits", rd(f"{A48}/charolais/loao.tsv"), "48.6: every refit."),
 ("A48 MDE, per cohort", rd(f"{A48}/power/mde_cohort.tsv"), "48.7: minimum detectable r at 80% power, two-sided alpha 0.05, per cohort and table."),
 ("A48 MDE, pooled", rd(f"{A48}/power/mde_pooled.tsv"), "48.7: fixed-effect minimum detectable r and expected mKH half-width at tau = 0 per pooled configuration."),
 ("A48 published effects", pe, "48.7: published effects converted to r (Li et al. extreme-group values; Nelore coefficient/SE with an assumed df) beside each cohort's interval. These conversions are not on the CLR scale, so the comparison is not used in the manuscript (draft 19)."),
 ("A48 attenuation", at, "48.7: hypothetical reliability scenarios: observed r = true r x sqrt(trait x feature reliability) with assumed reliabilities 0.5, 0.7 and 0.9. No reliability is estimated; the Charolais GTDB vs SILVA correlation is not a reliability estimate."),
 ("A48 mapping, MC1 only", rd(f"{A48}/mapping/mc1_pooled.tsv"), "48.8: labels and exclusions with MC1 cohort estimates only (no refit). excl_*_full: the statement with all cohorts."),
 ("A48 clade ratio grid", g.drop(columns=["analysis_no"]), "48.9: 18 cells (assignment margin x pseudo-count x minimum assigned share) x 4 analyses; clade_ratio_analysis decoded from the seed (seed = 4890000 + 100 x analysis + cell)."),
 ("A48 clade ratio influence", rd(f"{A48}/clade_ratio/influence.tsv"), "48.9: Wallace CH4 per kg DMI, current cell: Cook's distance refit and MM regression."),
 ("A48 clade ratio assignment", rd(f"{A48}/clade_ratio/assignment.tsv"), "48.9: share of ASVs assigned at each margin."),
]
assert all(len(n) <= 31 for n, _, _ in S), [n for n, _, _ in S if len(n) > 31]
nb = b2["Notes"].copy(); nb["Description"] = nb.Description.str.replace("the validated procedure", "the procedure that met the Amendment 47 coverage rule", regex=False).str.replace("the per-target validated range", "the per-target range over which the coverage rule was met", regex=False)
n2 = pd.concat([nb, pd.DataFrame({"Item": ["A48 sheets", "Figure 4, mixed outcomes", "A48 seven claims, p values"] + [n for n, _, _ in S],
    "Description": ["Sheets starting 'A48' are late analyses, run after all primary, late, Amendment 46 and Amendment 47 results were seen "
                    "(Amendment 48 and its addenda, in response to pre-submission assessments of the draft); they set no label. Source files in "
                    "results/late_a48/; every number the manuscript takes from them is in results/late_a48/manifest_a48.tsv.",
                    "The pooled row combining Difford g per day with Wallace per kg intake (different outcomes), removed from Figure 4 in draft 18, "
                    "is in 'Pooled results', version 'Methane, Wallace per kg DMI', target M1.", P48] + [t for _, _, t in S]})], ignore_index=True)
sheets2 = [("Notes", n2)] + [(k, v) for k, v in b2.items() if k != "Notes"] + [(n, f) for n, f, _ in S]
af104.write_book(OUT[2], sheets2, {"Per-cohort estimates": 3, "Pooled results": 3})
print("wrote", OUT[2], "| sheets", len(sheets2))

# ---------------------------------------------------------------- Additional file 3 (as 116, with the classification table first)
TITLE = "Additional file 3. Analysis plan, amendments, deviation table and analysis classification"
cls_md = open(P("docs", "ANALYSIS_CLASSES.md"), encoding="utf-8").read()
cls_md = re.sub(r"^# Analysis classification \(Amendment 48\.10\)", "# Part 0. Analysis classification (Amendment 48.10)", cls_md, count=1, flags=re.M)
dev = open(P("docs", "DEVIATIONS.md"), encoding="utf-8").read(); part1 = dev[dev.index("# How to read this file"):].rstrip() + "\n"
plan = open(P("plan", "ANALYSIS_PLAN.md"), encoding="utf-8").read()
part2 = "\n".join(("#" + l) if re.match(r"^#{1,5} ", l) else l for l in plan.splitlines()) + "\n"
md = cls_md.rstrip() + "\n\n" + part1 + "\n\n# Part 2. Analysis plan with all amendments (unedited)\n\n" + part2
with tempfile.TemporaryDirectory() as td:
    src = os.path.join(td, "af3.md"); open(src, "w", encoding="utf-8").write(md)
    ref = os.path.join(td, "reference.docx"); shutil.copy(os.path.join(D17, "Additional_file_3_analysis_plan.docx"), ref)
    tmp = os.path.join(td, "af3.docx")
    subprocess.run(["pandoc", src, "-f", "markdown-raw_html", "-t", "docx", "--reference-doc", ref, "-o", tmp], check=True)
    fixed = os.path.join(td, "af3_fixed.docx")
    with zipfile.ZipFile(tmp) as zi, zipfile.ZipFile(fixed, "w", zipfile.ZIP_DEFLATED) as zo:
        for it in zi.infolist():
            b = zi.read(it.filename)
            if it.filename == "docProps/core.xml":
                s = b.decode("utf-8")
                s = re.sub(r"<dc:title>.*?</dc:title>|<dc:title/>", f"<dc:title>{TITLE}</dc:title>", s) if "dc:title" in s else s.replace("</cp:coreProperties>", f"<dc:title>{TITLE}</dc:title></cp:coreProperties>")
                b = s.encode("utf-8")
            zo.writestr(it, b)
    shutil.move(fixed, OUT[3])
print("wrote", OUT[3], "| Part 0 docs/ANALYSIS_CLASSES.md, Part 1 docs/DEVIATIONS.md, Part 2 plan/ANALYSIS_PLAN.md (", len(plan.splitlines()), "lines )")

# ---------------------------------------------------------------- Additional file 4
src4 = P("results", "figures_draft18", "af4_seven_claims.pdf")
if not os.path.exists(src4): sys.exit(f"missing {src4}: run scripts/139_figures_draft18.R first")
shutil.copyfile(src4, OUT[4]); print("wrote", OUT[4])
