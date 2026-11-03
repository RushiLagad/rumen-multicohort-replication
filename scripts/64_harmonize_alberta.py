#!/usr/bin/env python3
# Amendment 19: harmonise Alberta (PRJNA393057) runs and animals after its genus table exists. Prints counts and design cross-tabs only, no phenotype values.
import csv, os, re, sys, zipfile, collections, xml.etree.ElementTree as ET
root = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort"); os.chdir(root)
ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}; z = zipfile.ZipFile("raw/phenotypes/Li2019_709_MOESM1.xlsx")
ss = ["".join(t.text or "" for t in si.iter("{%s}t" % ns["m"])) for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns)] if "xl/sharedStrings.xml" in z.namelist() else []
def val(c):
    if c.get("t") == "inlineStr": return "".join(t.text or "" for t in c.iter("{%s}t" % ns["m"]))
    v = c.find("m:v", ns); return "" if v is None else (ss[int(v.text)] if c.get("t") == "s" else v.text)
rows = [[val(c) for c in r.findall("m:c", ns)] for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).find("m:sheetData", ns).findall("m:row", ns)]
h = rows[2]; col = {k: h.index(k) for k in h}
need = {"animal_id": "Animal ID", "breed": "Breed", "sex": "Sex", "diet": "Diet", "age_d": "Age (day)", "dmi": "Dry matter intake (DMI; kg/day)", "adg": "Average daily gain (ADG; kg/day)", "rfi": "Residual feed intake (RFI; kg/day)"}
miss = [v for v in need.values() if v not in col]; assert not miss, f"columns not found: {miss}"
pheno = {}
for r in rows[3:]:
    if not r or not r[0]: continue
    r = r + [""] * (len(h) - len(r)); d = {k: r[col[v]].strip() for k, v in need.items()}
    d["breed"] = re.sub(r".*\((\w+)\)$", r"\1", d["breed"]); pheno[d["animal_id"]] = d
samp = {s["accession"]: s for s in csv.DictReader(open("metadata/candidates/PRJNA393057.samples.tsv"), delimiter="\t")}
excl = {r["run"] for r in csv.DictReader(open("results/16s/PRJNA393057/excluded_runs.tsv"), delimiter="\t")}
G = open("profiles_16s/alberta/alberta_genus_counts.tsv").readline().rstrip("\n").split("\t")[1:]
runs = []
for r in csv.DictReader(open("metadata/candidates/PRJNA393057.runs.tsv"), delimiter="\t"):
    s = samp[r["sample_accession"]]; runs.append({"run_accession": r["run_accession"], "animal_id": s["source_material_id"], "collection_date": s["collection_date"], "excluded": r["run_accession"] in excl, "in_genus_table": r["run_accession"] in G})
os.makedirs("harmonized", exist_ok=True)
with open("harmonized/alberta_runs.tsv", "w") as f:
    w = csv.DictWriter(f, list(runs[0]), delimiter="\t"); w.writeheader(); w.writerows(runs)
ids = collections.Counter(r["animal_id"] for r in runs); dup = [a for a, n in ids.items() if n > 1]
date = {r["animal_id"]: r["collection_date"] for r in runs}
use = [r for r in runs if r["in_genus_table"] and not r["excluded"]]
an = []
for r in use:
    p = pheno.get(r["animal_id"])
    if p is None: continue
    an.append(dict(p, collection_date=date[r["animal_id"]], has_rfi=p["rfi"] not in ("", "NA")))
with open("harmonized/alberta_animals.tsv", "w") as f:
    w = csv.DictWriter(f, list(an[0]), delimiter="\t"); w.writeheader(); w.writerows(an)
print("runs", len(runs), "| in genus table and not excluded", len(use), "| animal IDs repeated across runs", len(dup), "| runs with no phenotype row", sum(1 for r in use if r["animal_id"] not in pheno))
a = [x for x in an if x["has_rfi"]]; print("animals with microbiome and RFI:", len(a), "of", len(an))
def tab(k1, k2):
    t = collections.Counter((x[k1], x[k2]) for x in a); r1 = sorted({x[k1] for x in a}); r2 = sorted({x[k2] for x in a})
    print(f"  {k1} x {k2}:", "  ".join(f"{u}[" + " ".join(f"{t[(u, v)]}" for v in r2) + "]" for u in r1), " cols:", r2)
for k1, k2 in [("diet", "sex"), ("diet", "breed"), ("sex", "breed"), ("collection_date", "diet")]: tab(k1, k2)
