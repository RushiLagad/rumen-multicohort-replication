#!/usr/bin/env bash
# =============================================================================
# 03b_harmonize_metadata.sh
# One row per animal-sample, identical columns across cohorts, checksummed.
#
# Nelore join: S2 SampleID == ENA sample_alias, exact. Verified by a
# per-pair read-count check (S2 raw / 2 / ENA read_count must sit in a
# narrow band for every animal; a swapped label would fall outside it).
#
# RFI columns: rfi_value = RFI_GC (centered within slaughter batch, matches
# the published Table_7 associations); rfi_within_diet = RFI_by_dieta, the
# one to use for cross-cohort pooling because the diet variances differ ~3x.
# =============================================================================
set -uo pipefail
ROOT=/data/rumen_multicohort
mkdir -p "$ROOT/harmonized" "$ROOT/logs"

python3 -I - << 'PY'
import pandas as pd, numpy as np, os, sys, hashlib
ROOT="/data/rumen_multicohort"; OUT=f"{ROOT}/harmonized"
SCHEMA=["cohort","animal_id","breed","sex","diet","batch","weight_class",
        "sampling_method","rfi_value","rfi_within_diet","rfi_class",
        "methane_value","methane_within_diet","run_accession","sample_accession",
        "n_reads_published","n_reads_ena"]
blocks=[]
def rule(t): print("\n"+"="*70+"\n"+t+"\n"+"="*70)

# ====================================================== NELORE ========
rule("NELORE  PRJNA987743")
s2=pd.read_excel(f"{ROOT}/PRJNA987743/metadata/Table_2.XLSX",header=1)
s2.columns=[str(c).strip() for c in s2.columns]
s2["SampleID"]=s2.SampleID.astype(str).str.strip()
rum=s2[s2.Site.str.strip().str.lower()=="rumen"].copy()

ena=pd.read_csv(f"{ROOT}/PRJNA987743/metadata/rumen_runs.tsv",sep="\t",dtype=str).fillna("")
ena["sample_alias"]=ena.sample_alias.str.strip()
ena["read_count"]=pd.to_numeric(ena.read_count,errors="coerce")

j=rum.merge(ena[["sample_alias","run_accession","sample_accession","read_count"]],
            left_on="SampleID",right_on="sample_alias",how="left",validate="1:1")
miss=j.run_accession.isna().sum()
print(f"  S2 rumen rows {len(rum)}, ENA runs {len(ena)}, joined {len(j)-miss}, unmatched {miss}")
if miss: print("  UNMATCHED:",list(j[j.run_accession.isna()].SampleID)); sys.exit(2)

# per-pair sanity: raw/2/ena must be tight for every animal
j["_ratio"]=pd.to_numeric(j["#Raw_reads"])/2/j.read_count
lo,hi=j._ratio.quantile([.0,1.0])
print(f"  read-count ratio (S2 raw/2 over ENA): min {lo:.4f}  max {hi:.4f}  "
      f"spread {(hi-lo)*100:.2f}%")
bad=j[(j._ratio<0.95)|(j._ratio>1.10)]
if len(bad):
    print("  RATIO OUTLIERS, possible mislabel:"); print(bad[["SampleID","_ratio"]].to_string(index=False)); sys.exit(3)
print("  every pair inside [0.95, 1.10]: alias join confirmed")

j["_day"]=pd.to_datetime(j.Slaughter_day,unit="us",errors="coerce").dt.strftime("batch_%m%d")
print("\n  Diet x batch:"); print(pd.crosstab(j.Diet,j._day).to_string().replace("\n","\n    "))

blocks.append(pd.DataFrame({
  "cohort":"Nelore_PRJNA987743","animal_id":j.SampleID.str.replace("cr$","",regex=True),
  "breed":"Nelore","sex":"steer","diet":j.Diet.str.strip(),"batch":j._day,
  "weight_class":j.Weight.str.strip(),"sampling_method":"slaughter",
  "rfi_value":pd.to_numeric(j.RFI_GC),"rfi_within_diet":pd.to_numeric(j.RFI_by_dieta),
  "rfi_class":np.nan,
  "methane_value":pd.to_numeric(j.RME_GC),"methane_within_diet":pd.to_numeric(j.RME_by_dieta),
  "run_accession":j.run_accession,"sample_accession":j.sample_accession,
  "n_reads_published":pd.to_numeric(j["#Raw_reads"]).astype("Int64"),
  "n_reads_ena":j.read_count.astype("Int64"),
}))
print(f"\n  emitted {len(blocks[-1])} animals")

# =================================================== CHAROLAIS ========
rule("CHAROLAIS  PRJEB70895  binary arm, ENA crosswalk still unverified")
cm=pd.read_csv(f"{ROOT}/PRJEB70895/metadata/mapping.tab",sep="\t",dtype=str).fillna("")
cm.columns=[c.strip() for c in cm.columns]
cm=cm[cm.id!="20BAR-01326"]            # internally contradictory row, see log 03
ext=cm[cm.RFI_1.isin(["POSITIF","NEGATIF"])].copy()
print(f"  extreme-RFI samples {len(ext)} from {ext['N Animal'].nunique()} animals; "
      f"mid-RFI controls {(cm.RFI_1=='0').sum()}")

ce=pd.read_csv(f"{ROOT}/metadata/PRJEB70895_runs.tsv",sep="\t",dtype=str).fillna("")
print(f"  ENA: {len(ce)} runs, strategies {sorted(set(ce.library_strategy))}")
print("  ENA alias examples:",list(ce.sample_alias.unique()[:5]))
print("  ENA title examples:",list(ce.sample_title.unique()[:5]))
ov=lambda col: sum(v in set(ce[col]) for v in pd.concat([cm.id,cm.Name]))
print(f"  mapping ids found in ENA alias: {ov('sample_alias')}, title: {ov('sample_title')}")
print("  -> if both are 0, this arm needs the ENA sample XML before any reads are touched")

blocks.append(pd.DataFrame({
  "cohort":"Charolais_PRJEB70895","animal_id":ext["N Animal"],
  "breed":"Charolais","sex":"bull","diet":ext.Diet,"batch":ext.Year.str.replace("Year_","batch_"),
  "weight_class":np.nan,"sampling_method":ext.Sampling,
  "rfi_value":np.nan,"rfi_within_diet":np.nan,
  "rfi_class":ext.RFI_1.map({"POSITIF":"high_RFI","NEGATIF":"low_RFI"}),
  "methane_value":np.nan,"methane_within_diet":np.nan,
  "run_accession":"","sample_accession":"",
  "n_reads_published":pd.NA,"n_reads_ena":pd.NA,
}))
print(f"  emitted {len(blocks[-1])} samples (run accessions blank until crosswalk)")

# ======================================================= WRITE ========
rule("HARMONIZED")
h=pd.concat(blocks,ignore_index=True)[SCHEMA]
f=f"{OUT}/animal_metadata.tsv"; h.to_csv(f,sep="\t",index=False)
md5=hashlib.md5(open(f,"rb").read()).hexdigest()
open(f"{OUT}/animal_metadata.md5","w").write(f"{md5}  animal_metadata.tsv\n")
print(f"  {f}\n  md5 {md5}\n")
print(h.groupby("cohort").agg(rows=("animal_id","size"),animals=("animal_id","nunique"),
      continuous_rfi=("rfi_value",lambda s:s.notna().sum()),
      binary_rfi=("rfi_class",lambda s:s.notna().sum()),
      methane=("methane_value",lambda s:s.notna().sum()),
      with_runs=("run_accession",lambda s:(s!="").sum())).to_string())
print("\n  Li et al. PRJNA448333 not yet included. Give me the path to your")
print("  existing per-animal RFI table and it becomes a third block.")
PY
