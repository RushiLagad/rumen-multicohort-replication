#!/usr/bin/env bash
# =============================================================================
# 04_find_crosswalk.sh
# Find a real SampleID <-> run_accession key for PRJNA987743.
#
# Read counts cannot do it: S2 #Raw_reads / 2 / ENA read_count is a tight
# 1.036-1.043, so ENA holds read pairs of a lightly filtered version of the
# same libraries. Same samples, ~3.9% offset, no exact key. Rank-matching is
# unsafe because adjacent samples differ by 0.1% while the ratio wobbles 0.7%.
#
# Strategy, in order:
#   A dump all 52 S2 SampleIDs beside all 52 ENA aliases and look for structure
#   B pull every BioSample attribute for all 52 samples and search for the
#     S2 SampleID strings. Submitters usually deposit their own lab ID.
#   C optimal assignment on scaled read counts WITH a margin report, so we
#     know whether it is safe rather than assuming
#
# Metadata only.
# =============================================================================
set -uo pipefail

ROOT=/data/rumen_multicohort
NEL="$ROOT/PRJNA987743/metadata"
cd "$NEL"

# ---- pull full BioSample XML for every rumen sample -------------------------
if [ ! -s samples_full.xml ]; then
  ACC=$(awk -F'\t' 'NR>1 {print $2}' rumen_runs.tsv | sort -u | paste -sd,)
  echo "fetching BioSample XML for $(echo "$ACC" | tr ',' '\n' | wc -l) samples"
  curl -sSL --retry 4 --retry-delay 5 -o samples_full.xml \
    "https://www.ebi.ac.uk/ena/browser/api/xml/${ACC}"
fi
echo "samples_full.xml: $(stat -c%s samples_full.xml) bytes"
echo

python3 -I - << 'PY'
import pandas as pd, re, xml.etree.ElementTree as ET, collections, itertools
import numpy as np

NEL="/data/rumen_multicohort/PRJNA987743/metadata"

# ---------------------------------------------------------------- S2 ----
raw=pd.read_excel(f"{NEL}/Table_2.XLSX",header=None)
h=next(i for i in range(6) if raw.iloc[i].astype(str)
       .str.contains("SampleID",case=False,na=False).any())
s2=pd.read_excel(f"{NEL}/Table_2.XLSX",header=h)
s2.columns=[str(c).strip() for c in s2.columns]
s2["SampleID"]=s2.SampleID.astype(str).str.strip()
rum=s2[s2.Site.astype(str).str.strip().str.lower()=="rumen"].copy()
rum["_raw"]=pd.to_numeric(rum["#Raw_reads"],errors="coerce")

print("="*70); print("A  S2 SampleID structure, rumen rows"); print("="*70)
suf=rum.SampleID.str.extract(r'^(\d+)\s*([A-Za-z]*)$')
suf.columns=["num","suffix"]
print("  suffixes:",dict(collections.Counter(suf.suffix)))
nums=pd.to_numeric(suf.num,errors="coerce")
print(f"  numeric part: {nums.min():.0f} .. {nums.max():.0f}, {nums.nunique()} unique")
print("  all 52 rumen SampleIDs:")
print("   ", " ".join(sorted(rum.SampleID)))
print("\n  SampleID x Diet:")
print(pd.crosstab(suf.suffix,rum.Diet.values).to_string().replace("\n","\n    "))

# --------------------------------------------------------------- ENA ----
ena=pd.read_csv(f"{NEL}/rumen_runs.tsv",sep="\t",dtype=str).fillna("")
ena["read_count"]=pd.to_numeric(ena.read_count,errors="coerce")
print("\n"+"="*70); print("   ENA alias / title structure"); print("="*70)
print("  all 52 aliases:")
print("   ", " ".join(sorted(ena.sample_alias)))
print("  all 52 titles:")
print("   ", " ".join(sorted(set(ena.sample_title))))

# ------------------------------------------------- B BioSample attrs ----
print("\n"+"="*70); print("B  BioSample attributes, searching for S2 SampleIDs"); print("="*70)
try:
    tree=ET.parse(f"{NEL}/samples_full.xml"); root=tree.getroot()
except Exception as e:
    print("  XML unreadable:",e); root=None

recs=[]
if root is not None:
    for sm in root.iter("SAMPLE"):
        acc=sm.get("accession",""); alias=sm.get("alias","")
        d={"sample_accession":acc,"xml_alias":alias}
        for t in sm.iter("TITLE"):
            d["xml_title"]=(t.text or "").strip()
        for idt in sm.iter("EXTERNAL_ID"):
            d["external_"+(idt.get("namespace") or "x")]=(idt.text or "").strip()
        for sa in sm.iter("SAMPLE_ATTRIBUTE"):
            tg=sa.findtext("TAG","") or ""; vl=sa.findtext("VALUE","") or ""
            if tg: d["attr:"+tg.strip()]=vl.strip()
        recs.append(d)

if recs:
    X=pd.DataFrame(recs).fillna("")
    print(f"  parsed {len(X)} samples, {X.shape[1]} fields")
    print("\n  every field, with cardinality and examples:")
    for c in X.columns:
        v=X[c][X[c]!=""]
        if not len(v): continue
        print(f"    {c[:44]:46s} n={len(v):3d} uniq={v.nunique():3d} "
              f"{[str(x)[:22] for x in v.unique()[:3]]}")

    s2set=set(rum.SampleID)
    s2low={x.lower() for x in s2set}
    print("\n  --- does any field contain the S2 SampleIDs ---")
    hit=False
    for c in X.columns:
        vals=X[c].astype(str)
        exact=sum(1 for v in vals if v.strip().lower() in s2low)
        sub=sum(1 for v in vals if any(k in v.lower() for k in s2low))
        if exact or sub:
            hit=True
            print(f"    {c:46s} exact {exact:3d}  substring {sub:3d}  <== CROSSWALK")
            for v in vals[:4]: print(f"        {v[:60]}")
    if not hit:
        print("    no field contains any S2 SampleID, exact or substring")
else:
    X=pd.DataFrame()
    print("  no SAMPLE records parsed")

# ------------------------------- C assignment with margin check ---------
print("\n"+"="*70); print("C  scaled read-count assignment, with safety margin"); print("="*70)
a=rum[["SampleID","_raw","Diet","Weight"]].dropna(subset=["_raw"]).reset_index(drop=True)
b=ena[["run_accession","sample_alias","read_count"]].dropna(subset=["read_count"]).reset_index(drop=True)
ratio=np.median(np.sort(a._raw.values/2)/np.sort(b.read_count.values))
print(f"  median scale factor (S2/2 over ENA): {ratio:.5f}")
pred=a._raw.values/2/ratio                     # S2 mapped onto ENA scale
obs=b.read_count.values
D=np.abs(pred[:,None]-obs[None,:])/obs[None,:] # relative distance matrix

order=np.argsort(D,axis=1)
best=D[np.arange(len(a)),order[:,0]]
second=D[np.arange(len(a)),order[:,1]]
margin=second/np.maximum(best,1e-12)

print(f"  best-match relative error : median {np.median(best)*100:.3f}%  "
      f"max {best.max()*100:.3f}%")
print(f"  margin (2nd best / best)  : min {margin.min():.2f}x  "
      f"median {np.median(margin):.1f}x")
amb=(margin<3).sum()
print(f"  samples where 2nd best is within 3x of best: {amb}")
print(f"  assignment is 1:1 and consistent: "
      f"{len(set(order[:,0]))==len(a)}")

if amb==0 and len(set(order[:,0]))==len(a):
    print("\n  SAFE. Every S2 sample has one ENA run at least 3x closer than any")
    print("  alternative, and the assignment is a permutation. Usable as a join.")
else:
    print(f"\n  NOT SAFE on its own. {amb} ambiguous samples. Worst cases:")
    for i in np.argsort(margin)[:6]:
        j0,j1=order[i,0],order[i,1]
        print(f"    {a.SampleID[i]:>8s}  best {b.run_accession[j0]} "
              f"err {best[i]*100:.3f}%   runner-up {b.run_accession[j1]} "
              f"err {second[i]*100:.3f}%  margin {margin[i]:.2f}x")

pd.DataFrame({
  "SampleID":a.SampleID,
  "Diet":a.Diet,"Weight":a.Weight,
  "s2_raw_reads":a._raw.astype("Int64"),
  "proposed_run":b.run_accession.values[order[:,0]],
  "proposed_alias":b.sample_alias.values[order[:,0]],
  "ena_read_count":obs[order[:,0]].astype("int64"),
  "rel_err_pct":(best*100).round(4),
  "margin_x":margin.round(2),
}).sort_values("margin_x").to_csv(f"{NEL}/proposed_crosswalk.tsv",sep="\t",index=False)
print(f"\n  wrote {NEL}/proposed_crosswalk.tsv  (sorted worst margin first)")
print("\n  head:")
print(pd.read_csv(f"{NEL}/proposed_crosswalk.tsv",sep="\t").head(8).to_string(index=False))
PY
