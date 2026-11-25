#!/usr/bin/env python3
"""148: sample map, deposit run accession -> BioSample -> animal -> column of the analysed count table (draft 21 item 8; reporting only).

One row per run in the five public deposits used (PRJNA987743 Nelore, PRJEB70895 Charolais, PRJNA393057 Alberta, PRJEB28065 Difford,
PRJNA517480 Wallace), whether or not the run was analysed, with the reason when it was not. Sample identifiers in the count tables:
  Nelore, Charolais: run accession with the suffix "_clean" (the host-depleted reads); the two runs of Charolais animal 8574 are separate
                     columns and are summed into one animal profile before modelling (scripts/21_estimate_targets.R; Amendment 13);
  Alberta:           run accession (one run per animal);
  Difford:           the animal ID of the authors' phenotype table (Difford et al. 2018, S4), with all runs of the animal's BioSample summed;
  Wallace:           the deposit's sample alias (BREED_FARM_Sn), with the two sequencing rounds summed.
Inputs (all in the repository): metadata/ deposit run tables, harmonized/ run maps, profiles_16s/ count tables, results/16s/PRJNA393057/
excluded_runs.tsv; the GTDB count tables (profiles_gtdb_full/, Zenodo archive profiles_gtdb_full_counts.tar.gz) give the Nelore and
Charolais mRNA columns. Writes (refuses to overwrite): results/sample_map/sample_map.tsv and results/sample_map/summary.tsv.
Stops without writing unless every column of every analysed count table maps to exactly one deposit run or animal. Usage: 148_sample_map.py [root]
"""
import os, re, sys
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/data/rumen_multicohort"; P = lambda *a: os.path.join(ROOT, *a)
OD = P("results", "sample_map"); OUT = [os.path.join(OD, f) for f in ("sample_map.tsv", "summary.tsv")]
ex = [f for f in OUT if os.path.exists(f)]
if ex: sys.exit(f"output exists, never overwritten: {ex}")
rd = lambda *a: pd.read_csv(P(*a), sep="\t", dtype=str, keep_default_na=False)
cols = lambda *a: open(P(*a)).readline().rstrip("\n").split("\t")[1:]
T = {"nelore": ("profiles_gtdb_full", "nelore", "nelore_genus_counts.tsv"), "charolais_mrna": ("profiles_gtdb_full", "charolais_mrna", "charolais_mrna_genus_counts.tsv"),
     "charolais_16s": ("profiles_16s", "charolais_16s", "charolais_16s_genus_counts.tsv"), "alberta": ("profiles_16s", "alberta", "alberta_genus_counts.tsv"),
     "difford_bact": ("profiles_16s", "difford_bact", "difford_bact_genus_counts.tsv"), "difford_arch": ("profiles_16s", "difford_arch", "difford_arch_genus_counts.tsv"),
     "wallace_prokaryotes": ("profiles_16s", "wallace_prokaryotes", "wallace_prokaryotes_genus_counts.tsv"),
     "wallace_bacteria": ("profiles_16s", "wallace_bacteria", "wallace_bacteria_genus_counts.tsv"), "wallace_archaea": ("profiles_16s", "wallace_archaea", "wallace_archaea_genus_counts.tsv")}
C = {k: cols(*v) for k, v in T.items()}
for k, v in C.items(): assert len(v) == len(set(v)), f"duplicate columns in {k}"
tpath = lambda k: "/".join(T[k])
rows = []
def add(**r): rows.append(r)

# ---- Nelore and Charolais (shotgun metagenome, metatranscriptome, 16S rRNA gene)
hr = rd("harmonized", "runs.tsv"); hmap = hr.set_index("run_accession")
for bp, coh in (("PRJNA987743", "Nelore"), ("PRJEB70895", "Charolais")):
    for r in rd("metadata", f"{bp}_runs.tsv").itertuples():
        a = hmap.loc[r.run_accession] if r.run_accession in hmap.index else None
        aid = a.animal_id if a is not None else re.split(r"[-_]", r.sample_alias if coh == "Charolais" else r.sample_title)[0]
        if coh == "Nelore":
            dtype = "shotgun metagenome (" + ("rumen" if r.sample_title.endswith("Rumen") else "faeces") + ")"
            tab, col = "nelore", r.run_accession + "_clean"
        elif r.library_strategy == "RNA-Seq": dtype, tab, col = "metatranscriptome (rumen)", "charolais_mrna", r.run_accession + "_clean"
        else: dtype, tab, col = "16S rRNA gene amplicon (rumen)", "charolais_16s", r.run_accession
        used = col in C[tab]
        if used: note = "two runs of one BioSample summed into one animal profile before modelling (Amendment 13)" if (coh == "Charolais" and tab == "charolais_mrna" and (hr[(hr.layer == a.layer) & (hr.animal_id == aid)].shape[0] > 1)) else ""
        elif coh == "Nelore" and "faeces" in dtype: note = "faecal sample; the analysis uses rumen samples only"
        elif a is None: note = "not in the analysed run table (Amendment 13: animal 8508 has no 16S rRNA gene data; the 2020 sample of animal 8513 has contradictory sample fields)"
        else: note = "not in the count table"
        add(cohort=coh, bioproject=bp, run_accession=r.run_accession, biosample_accession=r.sample_accession, deposit_sample_alias=r.sample_alias,
            data_type=dtype, animal_id=aid, count_table=tpath(tab) if used else "", count_table_column=col if used else "", analysed=used, note=note)

# ---- Alberta (16S rRNA gene)
ar = rd("harmonized", "alberta_runs.tsv").set_index("run_accession"); exr = rd("results", "16s", "PRJNA393057", "excluded_runs.tsv").set_index("run")
for r in rd("metadata", "candidates", "PRJNA393057.runs.tsv").itertuples():
    a = ar.loc[r.run_accession] if r.run_accession in ar.index else None; used = r.run_accession in C["alberta"]
    note = "" if used else (f"excluded at read processing (rule {exr.loc[r.run_accession, 'rule']}; results/16s/PRJNA393057/excluded_runs.tsv)" if r.run_accession in exr.index else ("no animal ID in the deposit metadata" if a is None else "not in the count table"))
    add(cohort="Alberta", bioproject="PRJNA393057", run_accession=r.run_accession, biosample_accession=r.sample_accession, deposit_sample_alias=r.sample_alias,
        data_type="16S rRNA gene amplicon (rumen)", animal_id=a.animal_id if a is not None else "", count_table=tpath("alberta") if used else "",
        count_table_column=r.run_accession if used else "", analysed=used, note=note)

# ---- Difford (16S rRNA gene, bacterial and archaeal primers)
dr = rd("harmonized", "difford_runs.tsv").set_index("run")
for r in rd("metadata", "candidates", "PRJEB28065.runs.tsv").itertuples():
    d = dr.loc[r.run_accession]; tab = "difford_" + d.marker if d.marker in ("bact", "arch") else None
    used = d.in_s4 == "True" and tab is not None and d.cow in C[tab]
    if used: note = "all runs of the animal's BioSample summed before chimera removal" if (dr[(dr.cow == d.cow) & (dr.marker == d.marker) & (dr.in_s4 == "True")].shape[0] > 1) else ""
    elif d.in_s4 != "True": note = "no matching animal in the authors' phenotype table (Difford et al. 2018, S4); not analysed"
    elif tab is None: note = "primer pair not identified; not analysed"
    else: note = "animal not in the count table: excluded at read filtering (Amendment 22)"
    add(cohort="Difford", bioproject="PRJEB28065", run_accession=r.run_accession, biosample_accession=r.sample_accession, deposit_sample_alias=r.sample_alias,
        data_type=f"16S rRNA gene amplicon, {'bacterial' if d.marker == 'bact' else 'archaeal' if d.marker == 'arch' else 'unidentified'} primers (rumen)",
        animal_id=d.cow if d.in_s4 == "True" else "", count_table=tpath(tab) if used else "", count_table_column=d.cow if used else "", analysed=used, note=note)

# ---- Wallace (16S rRNA gene and other markers)
wr = rd("harmonized", "wallace_runs.tsv").set_index("run"); wb = rd("metadata", "candidates", "PRJNA517480.runs.tsv").set_index("run_accession")
MK = {"Prokaryotes": "wallace_prokaryotes", "Bacteria": "wallace_bacteria", "Archaea": "wallace_archaea"}
for r in rd("metadata", "candidates", "PRJNA517480.runs_full.tsv").itertuples():
    w = wr.loc[r.run_accession]; tab = MK.get(w.marker); used = tab is not None and w.cow in C[tab]
    two = (wr[(wr.cow == w.cow) & (wr.marker == w.marker)].shape[0] > 1)
    note = (("sequencing rounds SR1 and SR2 summed per animal before chimera removal" if two else "") if used else
            (f"{w.marker} marker; not analysed" if tab is None else "animal not in the count table: excluded at read filtering (Amendment 22)"))
    add(cohort="Wallace", bioproject="PRJNA517480", run_accession=r.run_accession, biosample_accession=wb.loc[r.run_accession, "sample_accession"] if r.run_accession in wb.index else "",
        deposit_sample_alias=r.sample_alias, data_type=f"amplicon, {w.marker} marker, round {w['round']} (rumen)", animal_id=w.cow,
        count_table=tpath(tab) if used else "", count_table_column=w.cow if used else "", analysed=used, note=note)

M = pd.DataFrame(rows)
# ---- every column of every analysed count table maps to the deposit; no animal appears in two columns of one table
bad = []
for k, v in C.items():
    got = set(M.loc[M.count_table == tpath(k), "count_table_column"])
    if got != set(v): bad.append((k, len(v), len(got), sorted(set(v) - got)[:5]))
    a = M[(M.count_table == tpath(k))].drop_duplicates(["count_table_column"]).groupby("animal_id").size()
    if (a > 1).sum() and k != "charolais_mrna": bad.append((k, "animal in more than one column", a[a > 1].index.tolist()[:5]))
if bad: sys.exit(f"count-table columns do not map one to one, nothing written: {bad}")
dupA = M[(M.count_table == tpath("charolais_mrna"))].drop_duplicates("count_table_column").groupby("animal_id").size()
assert dupA[dupA > 1].index.tolist() == ["8574"], dupA[dupA > 1]
assert not M.run_accession.duplicated().any()

S = []
for k in T:
    m = M[M.count_table == tpath(k)]
    S.append(dict(count_table=tpath(k), column_identifier={"nelore": "run accession + _clean", "charolais_mrna": "run accession + _clean", "charolais_16s": "run accession",
                  "alberta": "run accession"}.get(k, "animal ID (Difford S4)" if k.startswith("difford") else "deposit sample alias"),
                  columns=len(C[k]), analysed_runs=len(m), animals=m.animal_id.nunique(), biosamples=m.biosample_accession.nunique()))
for coh, g in M.groupby("cohort", sort=False):
    S.append(dict(count_table=f"{coh}: all deposit runs", column_identifier="", columns="", analysed_runs=f"{int(g.analysed.sum())} of {len(g)}",
                  animals=g[g.analysed].animal_id.nunique(), biosamples=g.biosample_accession.nunique()))
S = pd.DataFrame(S)

os.makedirs(OD, exist_ok=True)
for f, x in zip(OUT, (M, S)):
    tmp = f + ".tmp"; x.to_csv(tmp, sep="\t", index=False); os.replace(tmp, f)
print(S.to_string(index=False))
print("every count-table column maps to the deposit; wrote", ", ".join(os.path.relpath(f, ROOT) for f in OUT))
