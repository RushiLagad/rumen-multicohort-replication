# Amendment 29: random subsample of read pairs, mates kept together, seed fixed by run accession.
# Usage: 75_subsample_pairs.py in_1 in_2 out_1 out_2 k run
import gzip, random, sys, hashlib
r1, r2, o1, o2, k, run = sys.argv[1:]; k = int(k)
rng = random.Random(int(hashlib.md5(run.encode()).hexdigest()[:8], 16)); keep = []; n = 0
with gzip.open(r1, "rt") as a, gzip.open(r2, "rt") as b:
    while True:
        x = [a.readline() for _ in range(4)]; y = [b.readline() for _ in range(4)]
        if not x[0]: break
        assert y[0], f"{run}: R2 shorter than R1"
        if n < k: keep.append((x, y))
        else:
            j = rng.randrange(n + 1)
            if j < k: keep[j] = (x, y)
        n += 1
with gzip.open(o1, "wt", compresslevel=3) as a, gzip.open(o2, "wt", compresslevel=3) as b:
    for x, y in keep: a.writelines(x); b.writelines(y)
print(run, n, len(keep), sep="\t")
