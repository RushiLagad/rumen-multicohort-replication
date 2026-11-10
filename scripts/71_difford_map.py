import csv, gzip, re, collections, pandas as pd
R = list(csv.DictReader(open("metadata/candidates/PRJEB28065.runs.tsv"), delimiter="\t"))
cen = {x["file"].split("_")[0]: x["prefix"] for x in csv.DictReader(open("results/markers/PRJEB28065_prefix_census.tsv"), delimiter="\t") if x["file"].endswith("_1.fastq.gz")}
s4 = pd.read_csv("raw/phenotypes/Difford2018_S4.bin", sep=None, engine="python", dtype=str,
                 usecols=["sample_id", "Sequencing_Instrument", "Sequencing_Batch", "herd"])
S4 = s4.set_index("sample_id")
def marker(p):
    if p.startswith("GAGTTTGATCC"): return "bact"
    if re.match(r"AGC[AC]GCCGCGGTAA", p): return "arch"
    return "other"
rows = []
for x in R:
    r = x["run_accession"]
    with gzip.open(f"raw/PRJEB28065/{r}_1.fastq.gz", "rt") as h: f = h.readline().split()
    hd = f[1] if len(f) > 1 else ""; p = hd.split(":")
    inst = p[0] if len(p) > 2 else "NA"; fc = p[2] if len(p) > 2 else "NA"
    plat = "MiSeq" if inst.startswith("M") else ("HiSeq" if inst.startswith(("HWI-D", "D")) else "NA")
    cow = re.sub(r"_(lib\d+|Bact|Arch|V\d+)$", "", re.sub(r"^(IC|NG)-\d+_", "", x["sample_alias"]))
    rows.append(dict(run=r, sample=x["sample_accession"], alias=x["sample_alias"], cow=cow, in_s4=cow in S4.index,
                     marker=marker(cen.get(r, "")), platform=plat, instrument=inst, flowcell=fc, reads=x["read_count"]))
d = pd.DataFrame(rows); d.to_csv("harmonized/difford_runs.tsv", sep="\t", index=False)
print("runs:", len(d), "| marker:", dict(d.marker.value_counts()), "| platform:", dict(d.platform.value_counts()))
print("runs whose cow ID is in S4:", d.in_s4.sum(), "| unmapped examples:", list(d.alias[~d.in_s4].head(8)))
m = d[d.in_s4]
for mk in ["bact", "arch"]:
    c = m[m.marker == mk].groupby("cow").size()
    print(mk, "| S4 cows with >=1 run:", c.size, "of", len(S4), "| runs per cow:", dict(collections.Counter(c.values)))
print("flowcells:", d.flowcell.nunique(), "| flowcells per marker:", dict(d.groupby("marker").flowcell.nunique()))
pc = m.groupby("cow").platform.agg(lambda v: "+".join(sorted(set(v))))
print("S4 instrument code vs header platform (cows):"); print(pd.crosstab(S4.loc[pc.index, "Sequencing_Instrument"], pc))
fb = m.groupby("cow").flowcell.agg(lambda v: "+".join(sorted(set(v))))
print("S4 batch -> distinct flowcell sets:", dict(pd.DataFrame({"b": S4.loc[fb.index, "Sequencing_Batch"], "f": fb}).groupby("b").f.nunique()))
print("herd x Sequencing_Batch (cow counts):"); print(pd.crosstab(s4.herd, s4.Sequencing_Batch))
