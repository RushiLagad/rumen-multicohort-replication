#!/usr/bin/env bash
# =============================================================================
# 05_runs_table.sh
# Normalize: animal_metadata.tsv stays one row per animal; runs.tsv holds one
# row per (cohort, animal, layer, run). Needed because Charolais has 16S + mRNA
# per animal and Li et al. has DNA + total RNA + mRNA per animal.
# Also reconciles PRJEB70895 ENA aliases against mapping.tab in both directions.
# =============================================================================
set -uo pipefail
ROOT=/data/rumen_multicohort; H="$ROOT/harmonized"

python3 -I - << 'PY'
import pandas as pd, hashlib, re
ROOT="/data/rumen_multicohort"; H=f"{ROOT}/harmonized"
am=pd.read_csv(f"{H}/animal_metadata.tsv",sep="\t",dtype=str).fillna("")
runs=[]

# ---------------------------------------------------------------- Nelore
nel=am[am.cohort=="Nelore_PRJNA987743"]
for _,r in nel.iterrows():
    runs.append(dict(cohort=r.cohort,animal_id=r.animal_id,layer="DNA_shotgun",
                     run_accession=r.run_accession,sample_accession=r.sample_accession,
                     ena_alias=r.animal_id+"cr"))
print(f"Nelore: {len(nel)} animals -> {len(nel)} DNA runs")

# ------------------------------------------------------------- Charolais
cm=pd.read_csv(f"{ROOT}/PRJEB70895/metadata/mapping.tab",sep="\t",dtype=str).fillna("")
cm.columns=[c.strip() for c in cm.columns]
ce=pd.read_csv(f"{ROOT}/metadata/PRJEB70895_runs.tsv",sep="\t",dtype=str).fillna("")
ce["sample_alias"]=ce.sample_alias.str.strip()
print(f"\nCharolais ENA: {len(ce)} runs, {ce.sample_accession.nunique()} samples")
print(ce.groupby("library_strategy").size().to_string())

ena_al=set(ce.sample_alias); map_nm=set(cm.Name)
print(f"\n  ENA aliases: {len(ena_al)}   in mapping.tab: {len(ena_al&map_nm)}   "
      f"NOT in mapping.tab: {len(ena_al-map_nm)}")
for a in sorted(ena_al-map_nm): print("    unmapped ENA alias:",a)
print(f"  mapping Names with reads: {len(map_nm&ena_al)} of {len(map_nm)}")
print("  mapping rows WITHOUT reads, by Sampling x RFI_1:")
nr=cm[~cm.Name.isin(ena_al)]
print(pd.crosstab(nr.Sampling,nr.RFI_1).to_string().replace("\n","\n    "))

bad="8513-N-M-2-A20"
print(f"\n  contradictory row {bad} deposited on ENA: {bad in ena_al}")
if bad in ena_al:
    print("   ", ce[ce.sample_alias==bad][["run_accession","sample_accession","library_strategy"]]
          .to_string(index=False).replace("\n","\n    "))
    print("    -> its ENA sample XML may hold the true RFI class. Check before dropping.")

# one row per animal, abattoir + extreme + deposited only
cmx=cm[cm.RFI_1.isin(["POSITIF","NEGATIF"]) & cm.Name.isin(ena_al)].copy()
cmx=cmx[cmx.Name!=bad]
print(f"\n  usable Charolais animals: {cmx['N Animal'].nunique()}  "
      f"(extreme RFI, abattoir, reads on ENA, contradictory row excluded)")
print(pd.crosstab([cmx.Year,cmx.Diet],cmx.RFI_1).to_string().replace("\n","\n    "))

layer={"AMPLICON":"16S_amplicon","RNA-Seq":"mRNA"}
for _,r in cmx.iterrows():
    for _,e in ce[ce.sample_alias==r.Name].iterrows():
        runs.append(dict(cohort="Charolais_PRJEB70895",animal_id=r["N Animal"],
                         layer=layer.get(e.library_strategy,e.library_strategy),
                         run_accession=e.run_accession,sample_accession=e.sample_accession,
                         ena_alias=r.Name))

cha=pd.DataFrame({
  "cohort":"Charolais_PRJEB70895","animal_id":cmx["N Animal"],"breed":"Charolais","sex":"bull",
  "diet":cmx.Diet,"batch":cmx.Year.str.replace("Year_","batch_"),"weight_class":"",
  "sampling_method":"slaughter","rfi_value":"","rfi_within_diet":"",
  "rfi_class":cmx.RFI_1.map({"POSITIF":"high_RFI","NEGATIF":"low_RFI"}),
  "methane_value":"","methane_within_diet":"","run_accession":"","sample_accession":"",
  "n_reads_published":"","n_reads_ena":""})

# ----------------------------------------------------------------- write
am2=pd.concat([nel,cha[nel.columns]],ignore_index=True)
am2.to_csv(f"{H}/animal_metadata.tsv",sep="\t",index=False)
R=pd.DataFrame(runs); R.to_csv(f"{H}/runs.tsv",sep="\t",index=False)
with open(f"{H}/CHECKSUMS.md5","w") as fh:
    for f in ("animal_metadata.tsv","runs.tsv"):
        fh.write(f"{hashlib.md5(open(f'{H}/{f}','rb').read()).hexdigest()}  {f}\n")

print(f"\nanimal_metadata.tsv: {len(am2)} animals")
print(am2.groupby("cohort").size().to_string())
print(f"\nruns.tsv: {len(R)} runs")
print(R.groupby(["cohort","layer"]).size().to_string())
print(open(f"{H}/CHECKSUMS.md5").read())
PY
