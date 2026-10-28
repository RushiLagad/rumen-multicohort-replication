#!/usr/bin/env bash
# =============================================================================
# 03_harmonize_metadata.sh
# One row per ANIMAL, identical columns across cohorts, checksummed.
#
# The hard part is not the table, it is the joins. S2 keys on SampleID
# ("238cr"); ENA keys on sample_alias ("470_Rumen"). Different numbering
# systems, no shared ID. So this script joins on READ COUNT, which is an
# incidental fingerprint neither group chose, and REFUSES to emit anything
# unless that join is 1:1 and unambiguous.
#
# It also reports design aliasing per cohort, because a phenotype perfectly
# nested in batch is unidentifiable and is better discovered now than in
# review. That is what sank the identifiability argument in PRJNA448333.
#
# Nothing here touches reads. Safe while downloads run.
# =============================================================================
set -uo pipefail

ROOT=/data/rumen_multicohort
mkdir -p "$ROOT/harmonized" "$ROOT/logs"

python3 -I - << 'PY'
import pandas as pd, numpy as np, glob, os, sys, hashlib, itertools

ROOT="/data/rumen_multicohort"
OUT=os.path.join(ROOT,"harmonized")
SCHEMA=["cohort","animal_id","breed","diet","batch","sex","sampling_method",
        "rfi_value","rfi_class","methane_value","run_accession","n_reads"]
blocks=[]
def rule(t): print("\n"+"="*70+"\n"+t+"\n"+"="*70)

def find_header(df, token):
    """Locate the row that actually holds the column names."""
    for i in range(min(6,len(df))):
        if df.iloc[i].astype(str).str.contains(token, case=False, na=False).any():
            return i
    return None

# =====================================================================
rule("NELORE  PRJNA987743  Bos indicus, n=52")
# =====================================================================
s2p=os.path.join(ROOT,"PRJNA987743/metadata/Table_2.XLSX")
enap=os.path.join(ROOT,"PRJNA987743/metadata/rumen_runs.tsv")

if not os.path.exists(s2p):
    print("  MISSING",s2p); sys.exit(1)
if not os.path.exists(enap):
    print("  MISSING",enap,"\n  run: bash scripts/01_fetch_nelore.sh"); sys.exit(1)

raw=pd.read_excel(s2p,header=None)
h=find_header(raw,"SampleID")
s2=pd.read_excel(s2p,header=h)
s2.columns=[str(c).strip() for c in s2.columns]
print(f"  S2: {s2.shape[0]} rows, header row {h}")

rum=s2[s2.Site.astype(str).str.strip().str.lower()=="rumen"].copy()
print(f"  rumen rows: {len(rum)}")

ena=pd.read_csv(enap,sep="\t",dtype=str).fillna("")
ena["read_count"]=pd.to_numeric(ena.read_count,errors="coerce")
print(f"  ENA rumen runs: {len(ena)}")

# ---- join on read count, testing both conventions -------------------
rum["_raw"]=pd.to_numeric(rum["#Raw_reads"],errors="coerce")
print("\n  --- read-count join test ---")
print(f"    S2  #Raw_reads  : {rum._raw.min():,.0f} .. {rum._raw.max():,.0f}")
print(f"    ENA read_count  : {ena.read_count.min():,.0f} .. {ena.read_count.max():,.0f}")

best=None
for name,factor in [("1x (spots = reads)",1.0),("2x (spots = pairs)",2.0)]:
    key=(rum._raw/factor).round().astype("Int64")
    m=key.isin(set(ena.read_count.dropna().astype("int64")))
    dup_s2=key.duplicated(keep=False).sum()
    dup_en=ena.read_count.duplicated(keep=False).sum()
    print(f"    {name:20s} matched {m.sum():3d}/{len(rum)}   "
          f"dup keys in S2: {dup_s2}  in ENA: {dup_en}")
    if best is None or m.sum()>best[1]:
        best=(factor,m.sum(),dup_s2,dup_en,name)

factor,nmatch,dup_s2,dup_en,label=best
if nmatch!=len(rum) or dup_s2 or dup_en:
    print(f"\n  JOIN REJECTED. best was {label}: {nmatch}/{len(rum)} matched, "
          f"{dup_s2} duplicate S2 keys, {dup_en} duplicate ENA keys.")
    print("  Not emitting a table on an unverified join. Fallback options,")
    print("  in order of how much I would trust them:")
    print("    1 join on fastq_bytes instead (also incidental, also unique)")
    print("    2 email the corresponding author for the SampleID to BioSample map")
    print("    3 match on (Diet, Weight) cells only, which is group-level and")
    print("      destroys the per-animal pairing the whole study depends on")
    print("\n  Diagnostics to send me:")
    print("    S2  head :", sorted(rum._raw.dropna().astype(int))[:6])
    print("    ENA head :", sorted(ena.read_count.dropna().astype(int))[:6])
    sys.exit(2)

print(f"\n  JOIN ACCEPTED via {label}, 1:1 on all {len(rum)} animals")
rum["_key"]=(rum._raw/factor).round().astype("Int64")
ena["_key"]=ena.read_count.astype("Int64")
j=rum.merge(ena[["_key","run_accession","sample_alias","sample_title"]],on="_key",how="left")
assert len(j)==len(rum) and j.run_accession.notna().all()

# ---- design diagnostics ---------------------------------------------
j["_day"]=pd.to_datetime(j.Slaughter_day,unit="us",errors="coerce").dt.date
print("\n  --- design structure ---")
print("    slaughter days:",sorted(set(j._day.dropna())))
ct=pd.crosstab(j.Diet,j._day)
print("\n    Diet x slaughter_day:")
print(ct.to_string().replace("\n","\n      "))
alias=(ct>0).sum(0).max()==1
print(f"\n    diet aliased with slaughter day: {'YES, PROBLEM' if alias else 'no, both diets on shared days'}")
print("\n    Diet x Weight:")
print(pd.crosstab(j.Diet,j.Weight).to_string().replace("\n","\n      "))
for c in ("RFI_GC","RME_GC","RFI_by_dieta","RME_by_dieta"):
    v=pd.to_numeric(j[c],errors="coerce")
    print(f"    {c:14s} n={v.notna().sum():3d}  mean {v.mean():+.3f}  sd {v.std():.3f}  "
          f"[{v.min():+.2f},{v.max():+.2f}]")

nel=pd.DataFrame({
  "cohort":"Nelore_PRJNA987743",
  "animal_id":j.SampleID.astype(str).str.strip(),
  "breed":"Nelore",
  "diet":j.Diet.astype(str).str.strip(),
  "batch":j._day.astype(str),
  "sex":"steer",
  "sampling_method":"slaughter",
  "rfi_value":pd.to_numeric(j.RFI_GC,errors="coerce"),
  "rfi_class":np.nan,
  "methane_value":pd.to_numeric(j.RME_GC,errors="coerce"),
  "run_accession":j.run_accession,
  "n_reads":j._raw.astype("Int64"),
})
blocks.append(nel); print(f"\n  emitted {len(nel)} animals")

# =====================================================================
rule("CHAROLAIS  PRJEB70895  binary replication arm")
# =====================================================================
mp=os.path.join(ROOT,"PRJEB70895/metadata/mapping.tab")
cm=pd.read_csv(mp,sep="\t",dtype=str).fillna("")
cm.columns=[c.strip() for c in cm.columns]
print(f"  mapping.tab: {len(cm)} rows, {cm['N Animal'].nunique()} animals")

BAD="20BAR-01326"
if (cm.id==BAD).any():
    r=cm[cm.id==BAD].iloc[0]
    print(f"\n  DROPPING {BAD}  name={r['Name']}  RFI_1={r['RFI_1']}  Diet={r['Diet']}")
    print("    name encodes N and M, fields say POSITIF and GS, and animal 8513")
    print("    is the only ID in both years while every other 2020 animal is 9xxx.")
    cm=cm[cm.id!=BAD]

ext=cm[cm.RFI_1.isin(["POSITIF","NEGATIF"])].copy()
mid=cm[cm.RFI_1=="0"]
print(f"\n  extreme-RFI samples: {len(ext)} from {ext['N Animal'].nunique()} animals")
print(f"  mid-RFI samples    : {len(mid)} from {mid['N Animal'].nunique()} animals  (monotonicity control)")
print("\n  Year x Diet x class:")
print(pd.crosstab([ext.Year,ext.Diet],ext.RFI_1).to_string().replace("\n","\n    "))
print("\n  paired sampling within animal:")
pr=ext.groupby("N Animal").Sampling.nunique()
print(f"    animals with both Sonde and Abattoir: {(pr>1).sum()}")

# ---- the 80-vs-64 reconciliation ------------------------------------
cena=os.path.join(ROOT,"metadata/PRJEB70895_runs.tsv")
if os.path.exists(cena):
    ce=pd.read_csv(cena,sep="\t",dtype=str).fillna("")
    print(f"\n  --- ENA reconciliation ---")
    print(f"    mapping.tab rows {len(cm)+1} (pre-drop)   ENA runs {len(ce)}   "
          f"ENA samples {ce.sample_accession.nunique()}")
    print(f"    ENA strategies: {sorted(set(ce.library_strategy))}")
    for strat,g in ce.groupby("library_strategy"):
        print(f"      {strat:12s} {len(g):3d} runs")
    hit=0
    for col in ("sample_alias","sample_title"):
        if col in ce:
            s=set(ce[col])
            n=sum(1 for v in cm.id if v in s)+sum(1 for v in cm.Name if v in s)
            print(f"    mapping ids found in ENA {col}: {n}")
            hit=max(hit,n)
    if hit==0:
        print("    NO overlap on id or Name. The 16S/mRNA deposit uses its own")
        print("    accession labels, so this arm needs a crosswalk before use.")
        print("    Check PRJEB70895 sample XML for a 'sample name' attribute.")
else:
    print(f"\n  ENA run table missing: {cena}")
    print("  run: bash scripts/00b_cohort_gate_fixed.sh")

cha=pd.DataFrame({
  "cohort":"Charolais_PRJEB70895",
  "animal_id":ext["N Animal"].astype(str)+"_"+ext.Sampling.astype(str),
  "breed":"Charolais",
  "diet":ext.Diet.astype(str),
  "batch":ext.Year.astype(str),
  "sex":"bull",
  "sampling_method":ext.Sampling.astype(str),
  "rfi_value":np.nan,
  "rfi_class":ext.RFI_1.map({"POSITIF":"high_RFI","NEGATIF":"low_RFI"}),
  "methane_value":np.nan,
  "run_accession":"",
  "n_reads":pd.NA,
})
blocks.append(cha); print(f"\n  emitted {len(cha)} samples")

# =====================================================================
rule("BASELINE  PRJNA448333  Li et al., n=35 with continuous RFI")
# =====================================================================
cands=sorted(set(itertools.chain.from_iterable(
  glob.glob(p) for p in [
    "/data/**/PRJNA448333*meta*", "/data/**/*448333*/**/*meta*",
    "/data/**/*li_et_al*", "/data/**/*rfi*.csv", "/data/**/*rfi*.tsv",
    "/data/**/harmonized*/*.csv", "/data/**/*phenotype*",
  ])))
cands=[c for c in cands if os.path.isfile(c) and os.path.getsize(c)<5_000_000]
print("  Your existing Li et al. phenotype file is not at a path I know.")
print("  Candidates on this disk:")
for c in cands[:25]: print("   ",c)
if not cands: print("    none found")
print("\n  Point me at the right one and I will wire it in. It needs animal ID,")
print("  continuous RFI, breed, and collection date. Without it the pooled")
print("  continuous analysis is 52 animals, not 87, and the power floor fails.")

# =====================================================================
rule("HARMONIZED TABLE")
# =====================================================================
h=pd.concat(blocks,ignore_index=True)[SCHEMA]
f=os.path.join(OUT,"animal_metadata.tsv")
h.to_tsv=None
h.to_csv(f,sep="\t",index=False)
md5=hashlib.md5(open(f,'rb').read()).hexdigest()
open(os.path.join(OUT,"animal_metadata.md5"),"w").write(f"{md5}  animal_metadata.tsv\n")

print(f"  wrote {f}")
print(f"  md5  {md5}")
print(f"  {len(h)} rows x {len(SCHEMA)} cols\n")
print(h.groupby("cohort").agg(
    n=("animal_id","size"),
    continuous_rfi=("rfi_value",lambda s:s.notna().sum()),
    binary_rfi=("rfi_class",lambda s:s.notna().sum()),
    methane=("methane_value",lambda s:s.notna().sum()),
).to_string())

ncont=h.rfi_value.notna().sum()
print(f"\n  continuous-RFI animals available now: {ncont}")
print(f"  power floor for r=0.3 at 80%: 85")
print("  VERDICT:", "clears" if ncont>=85 else
      f"SHORT BY {85-ncont}. PRJNA448333 must be wired in, or the pooled")
if ncont<85:
    print("           continuous analysis is underpowered and the paper needs")
    print("           to lead with the per-cohort replication framing instead.")
PY
