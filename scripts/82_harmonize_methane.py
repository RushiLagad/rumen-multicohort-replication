# Amendments 28, 34, 36: cow-level phenotype tables for Difford and Wallace. Prints counts only.
import csv, collections, pandas as pd
s4 = pd.read_csv("raw/phenotypes/Difford2018_S4.bin", sep=None, engine="python", dtype=str)
s4.columns = ["cow", "instrument", "batch", "herd", "parity_raw", "dim", "ch4_gd"]
par = pd.to_numeric(s4.parity_raw, errors="coerce"); s4["parity"] = par.clip(upper=3).astype("Int64").astype(str).replace("3", "3+")
d = s4[["cow", "herd", "batch", "instrument", "parity", "dim", "ch4_gd"]]
d.to_csv("harmonized/difford_animals.tsv", sep="\t", index=False, lineterminator="\n")
print("difford: cows", len(d), "| with CH4", pd.to_numeric(d.ch4_gd, errors="coerce").notna().sum(), "| parity", dict(collections.Counter(d.parity)), "| herds", d.herd.nunique(), "| batches", d.batch.nunique())
x = pd.read_excel("raw/phenotypes/Wallace2019_Data_S1.xlsx", sheet_name="Phenotype data", header=None, dtype=str).iloc[3:]
x = x[x[1].notna()]
w = pd.DataFrame({"cow": ["%s_%s_S%s" % ("NR" if str(b).strip() == "2" else "HOL", str(f).strip(), str(i).strip()) for b, f, i in zip(x[6], x[4], x[1])],
                  "farm": x[4].str.strip(), "breed": x[6].str.strip(), "dim": x[7], "ch4_gd": x[68], "ch4_dmi": x[69], "rfi_nrc": x[52], "method": x[67].str.strip()})
par = pd.to_numeric(x[8], errors="coerce"); w["parity"] = par.clip(upper=3).astype("Int64").astype(str).replace("3", "3+").values
runs = {r["cow"] for r in csv.DictReader(open("harmonized/wallace_runs.tsv"), delimiter="\t")}
w["sequenced"] = w.cow.isin(runs); w.to_csv("harmonized/wallace_animals.tsv", sep="\t", index=False, lineterminator="\n")
print("wallace: rows", len(w), "| sequenced cows matched", w.sequenced.sum(), "of", len(runs), "| with CH4", pd.to_numeric(w.ch4_gd[w.sequenced], errors="coerce").notna().sum(), "| parity", dict(collections.Counter(w.parity[w.sequenced])))
