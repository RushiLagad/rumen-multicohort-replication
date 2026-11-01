#!/usr/bin/env python3
"""Build plan/target_keys.tsv: map the pre-specified targets (plan Section 3) to GTDB r207 names present in
k2_gtdb207, by the name/suffix rule of Amendment 9. Inputs: /tmp/k2_inspect.txt (kraken2-inspect), Table_7.XLSX."""
import re, collections, pandas as pd
names = collections.defaultdict(set)
for line in open('/tmp/k2_inspect.txt'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6 and f[3].strip() in ('P', 'F', 'G', 'S'): names[f[3].strip()].add(f[5].strip())
def keys(t, rank, alt=None):
    pat = re.compile(r'^(' + '|'.join(re.escape(x) for x in [t] + (alt or [])) + r')(_[A-Z]+)?$')
    return sorted(n for n in names[rank] if pat.match(n))
def skeys(genus, epithet):
    pat = re.compile(r'^' + re.escape(genus) + r'(_[A-Z]+)? ' + re.escape(epithet) + r'(_[A-Z]+)?$')
    return sorted(n for n in names['S'] if pat.match(n))
rows = []
def add(set_, id_, target, trank, site, pheno, direction, coef, se, q, mrank, ks, status, note=''):
    rows.append(dict(set=set_, id=id_, target=target, target_rank=trank, site=site, phenotype=pheno, direction=direction,
        coef=coef, stderr=se, qval=q, match_rank=mrank, gtdb_keys=';'.join(ks), n_keys=len(ks),
        status=status if ks or status != 'mapped' else 'unmapped', note=note))
# 3.1 Li et al. 2019
add('Li', 'T1', 'unnamed RF19 (text: RF16)', '16S clade', 'rumen', 'RFI', 'higher in L-RFI (ANG)', '', '', 0.045, '', [], 'untestable_by_design', 'Greengenes clade label; no genome-based counterpart')
add('Li', 'T2', 'Firmicutes', 'phylum', 'rumen', 'RFI', 'higher in L-RFI (CHAR)', '', '', 0.042, 'P', keys('Firmicutes', 'P', ['Bacillota']), 'mapped', 'NCBI name now Bacillota; GTDB splits into Firmicutes + suffixed phyla, summed')
add('Li', 'T3', 'Chloroflexi', 'phylum', 'rumen', 'RFI', 'not stated (CHAR)', '', '', 0.046, 'P', keys('Chloroflexota', 'P', ['Chloroflexi']), 'mapped', 'NCBI/GTDB name Chloroflexota')
add('Li', 'T4', 'unclassified Methanomassiliicoccaceae', 'family', 'rumen', 'RFI', 'higher in H-RFI (CHAR)', '', '', 0.014, 'F', keys('Methanomassiliicoccaceae', 'F'), 'mapped', 'tested as the whole family; flagged a priori (few reference genomes)')
add('Li', 'T5', 'Blautia', 'genus', 'rumen', 'RFI', 'higher in L-RFI (HYB)', '', '', 0.035, 'G', keys('Blautia', 'G'), 'mapped')
add('Li', 'T6', 'Methanobrevibacter smithii', 'species', 'rumen', 'RFI', 'higher in L-RFI (HYB)', '', '', 0.040, 'S', skeys('Methanobrevibacter', 'smithii'), 'mapped', 'species-level sensitivity space')
# 3.2 Nelore Table 7
t7 = pd.read_excel('PRJNA987743/metadata/Table_7.XLSX', header=1)
t7.columns = [c.strip() for c in t7.columns]
for i, r in t7.iterrows():
    g = str(r['Feature']).strip(); ph = str(r['Metadata']).strip(); c = float(r['coef'])
    d = ('higher' if c > 0 else 'lower') + ' with higher ' + ph.replace('_GC', '')
    note = 'faecal association; cross-site test in rumen data' if str(r['Site']).strip().lower().startswith('fec') else ''
    add('Nelore', f'N{i+1:02d}', g, 'genus', str(r['Site']).strip().lower(), ph, d, c, float(r['stderr']), float(r['qval']), 'G', keys(g, 'G'), 'mapped', note)
# 3.3 canonical
for i, (t, rk) in enumerate([('Methanobrevibacter', 'G'), ('Succinivibrionaceae', 'F'), ('Fibrobacter', 'G'), ('Ruminococcus', 'G'), ('Prevotella', 'G'), ('Butyrivibrio', 'G')], 1):
    add('Canonical', f'C{i}', t, 'family' if rk == 'F' else 'genus', 'rumen', 'RFI', 'literature candidate', '', '', '', rk, keys(t, rk), 'mapped')
df = pd.DataFrame(rows); df.to_csv('plan/target_keys.tsv', sep='\t', index=False)
print(f"{len(df)} targets written to plan/target_keys.tsv")
for s, g in df.groupby('set', sort=False):
    print(f"\n{s}: {(g.status=='mapped').sum()} mapped, {(g.status=='unmapped').sum()} unmapped, {(g.status=='untestable_by_design').sum()} untestable by design")
    for _, r in g.iterrows():
        if r.status == 'unmapped': print(f"  unmapped: {r.id} {r.target} ({r.site}, {r.phenotype})")
        elif r.n_keys > 1: print(f"  {r.id} {r.target} -> {r.n_keys} keys: {r.gtdb_keys}")
print("\nfaecal-site Nelore targets:", (df.site == 'feces').sum())
