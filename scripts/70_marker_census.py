import gzip, sys, os, glob, collections, statistics
from multiprocessing import Pool
P = sys.argv[1]; N = 2000
def one(f):
    c = collections.Counter(); L = []
    with gzip.open(f, "rt") as h:
        for i, l in enumerate(h):
            if i % 4 == 1:
                s = l.strip(); c[s[1:16]] += 1; L.append(len(s))
                if len(L) >= N: break
    if not L: return (os.path.basename(f), "EMPTY", 0, 0, 0)
    p, n = c.most_common(1)[0]
    return (os.path.basename(f), p, round(n/len(L), 3), int(statistics.median(L)), len(L))
fs = sorted(glob.glob(f"raw/{P}/*.fastq.gz"))
with Pool(32) as pool: res = pool.map(one, fs, chunksize=8)
with open(f"results/markers/{P}_prefix_census.tsv", "w") as o:
    o.write("file\tprefix\tfrac\tmedlen\tn\n")
    for x in res: o.write("\t".join(map(str, x)) + "\n")
mate = lambda f: "R2" if f.endswith("_2.fastq.gz") else ("R1" if f.endswith("_1.fastq.gz") else "SE")
print(P, "files:", len(res))
for (m, p), n in collections.Counter((mate(x[0]), x[1]) for x in res).most_common(20): print(m, p, n)
runs = collections.defaultdict(dict)
for x in res: runs[x[0].split("_")[0].split(".")[0]][mate(x[0])] = x[1]
print("R1/R2 prefix pairs:")
for k, n in collections.Counter((d.get("R1", "-"), d.get("R2", "-")) for d in runs.values() if "R1" in d).most_common(10): print(k, n)
print("files with modal-prefix fraction < 0.5:", sum(x[2] < 0.5 for x in res))
