#!/usr/bin/env python3
"""Merge per-sample Bracken outputs into taxon x sample count matrices (new_est_reads)."""
import sys, glob, os, pandas as pd
d, c = sys.argv[1], sys.argv[2]
for lvl, name in (('G', 'genus'), ('S', 'species')):
    cols = {}
    for f in sorted(glob.glob(f"{d}/*.bracken_{lvl}")):
        s = os.path.basename(f).replace(f".bracken_{lvl}", "")
        t = pd.read_csv(f, sep='\t'); cols[s] = t.set_index('name')['new_est_reads']
    if not cols: continue
    m = pd.DataFrame(cols).fillna(0).astype(int); m.index.name = 'taxon'
    m.to_csv(f"{d}/{c}_{name}_counts.tsv", sep='\t')
    print(f"{c} {name}: {m.shape[0]} taxa x {m.shape[1]} samples -> {c}_{name}_counts.tsv")
