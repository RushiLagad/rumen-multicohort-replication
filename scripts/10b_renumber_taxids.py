#!/usr/bin/env python3
"""Renumber Struo2/GTDB taxids to small sequential ints (root=1). Reason: kraken2 build_db
dropped every sequence whose leaf taxid > 2^31 (mini-test 2026-10-30; 31,330/62,443 genomes affected).
Keeps taxonomy/orig/* and seqid2taxid.map.orig_taxids; writes taxonomy/taxid_remap.tsv (old<TAB>new)."""
import sys, os, shutil, subprocess
db = sys.argv[1]; tax = os.path.join(db, 'taxonomy'); orig = os.path.join(tax, 'orig')
if os.path.exists(os.path.join(tax, '.renumbered')): print('already renumbered'); sys.exit(0)
os.makedirs(orig, exist_ok=True)
for f in ('nodes', 'names', 'merged', 'delnodes'):
    src = os.path.join(tax, f + '.dmp')
    if os.path.exists(src) and not os.path.exists(os.path.join(orig, f + '.dmp')): shutil.copy2(src, orig)
m_orig = os.path.join(db, 'seqid2taxid.map.orig_taxids')
if not os.path.exists(m_orig): shutil.copy2(os.path.join(db, 'seqid2taxid.map'), m_orig)
remap = {'1': '1'}; n = 1; nodes = []
for line in open(os.path.join(orig, 'nodes.dmp')):
    f = line.rstrip('\n').split('\t|\t')
    if f[0] not in remap: n += 1; remap[f[0]] = str(n)
    nodes.append(f)
with open(os.path.join(tax, 'taxid_remap.tsv'), 'w') as o:
    o.write('old_taxid\tnew_taxid\n'); [o.write(f'{k}\t{v}\n') for k, v in remap.items()]
with open(os.path.join(tax, 'nodes.dmp'), 'w') as o:
    for f in nodes: f[0] = remap[f[0]]; f[1] = remap[f[1]]; o.write('\t|\t'.join(f) + '\n')
skipped_names = 0
with open(os.path.join(tax, 'names.dmp'), 'w') as o:
    for line in open(os.path.join(orig, 'names.dmp')):
        f = line.split('\t|\t', 1)
        if f[0] in remap: f[0] = remap[f[0]]; o.write('\t|\t'.join(f))
        else: skipped_names += 1
open(os.path.join(tax, 'merged.dmp'), 'w').close(); open(os.path.join(tax, 'delnodes.dmp'), 'w').close()
unmapped = 0; seen = set(); written = 0
with open(os.path.join(db, 'seqid2taxid.map'), 'w') as o:
    for line in open(m_orig):
        s, t = line.rstrip('\n').split('\t')
        if s in seen: continue
        seen.add(s)
        if t in remap: o.write(f'{s}\t{remap[t]}\n'); written += 1
        else: unmapped += 1
open(os.path.join(tax, '.renumbered'), 'w').close()
print(f'nodes {len(nodes)}, max new taxid {n}, names skipped {skipped_names}, map: {written} unique contigs written, {unmapped} unmapped')
