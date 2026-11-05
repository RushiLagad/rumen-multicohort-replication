# Amendment 28 items 6-8: Wallace (PRJNA517480) run map from NCBI library names. No phenotype data.
import csv, collections, statistics
lab = list(csv.DictReader(open("metadata/candidates/PRJNA517480.runs_full.tsv"), delimiter="\t"))
cen = {x["file"].split(".")[0].split("_")[0]: (x["file"], int(x["medlen"])) for x in csv.DictReader(open("results/markers/PRJNA517480_prefix_census.tsv"), delimiter="\t")}
rows, ml, bad = [], collections.defaultdict(list), 0
for x in lab:
    p = x["library_name"].split("_"); r = x["run_accession"]; fn, L = cen[r]; ml[p[3]].append(L)
    bad += x["sample_alias"] != "_".join(p[:3])
    rows.append(dict(run=r, cow=x["sample_alias"], breed=p[0], farm=p[1], marker=p[3], round=p[4], file="raw/PRJNA517480/" + fn, reads=x["read_count"]))
mm = {m: int(statistics.median(v)) for m, v in ml.items()}
for d in rows: d["marker_median_len"] = mm[d["marker"]]
with open("harmonized/wallace_runs.tsv", "w", newline="") as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
print("runs:", len(rows), "| alias != library-name cow:", bad, "| census median length per marker:", mm)
print("farms:", dict(collections.Counter(d["farm"] for d in rows if d["marker"] == "Prokaryotes")))
