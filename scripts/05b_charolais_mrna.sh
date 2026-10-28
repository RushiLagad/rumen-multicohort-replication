#!/usr/bin/env bash
# Add Charolais mRNA runs to runs.tsv (alias = bare animal id, or id-year for
# the duplicated 8513) and inspect the one contradictory BioSample.
set -uo pipefail
ROOT=/data/rumen_multicohort; H="$ROOT/harmonized"; M="$ROOT/PRJEB70895/metadata"

curl -sSL --retry 3 -o "$M/SAMEA114804964.xml" \
  "https://www.ebi.ac.uk/ena/browser/api/xml/SAMEA114804964"
echo "=== BioSample SAMEA114804964 (the 2020 '8513') ==="
grep -oE "<TAG>[^<]+</TAG>\s*<VALUE>[^<]*</VALUE>" "$M/SAMEA114804964.xml" \
  | sed -E 's/<TAG>([^<]+)<\/TAG>\s*<VALUE>([^<]*)<\/VALUE>/  \1 = \2/' || echo "  (no attributes)"
grep -oE 'alias="[^"]+"|<TITLE>[^<]*' "$M/SAMEA114804964.xml" | sed 's/^/  /'
echo

python3 -I - << 'PY'
import pandas as pd, hashlib
ROOT="/data/rumen_multicohort"; H=f"{ROOT}/harmonized"
am=pd.read_csv(f"{H}/animal_metadata.tsv",sep="\t",dtype=str).fillna("")
R=pd.read_csv(f"{H}/runs.tsv",sep="\t",dtype=str).fillna("")
ce=pd.read_csv(f"{ROOT}/metadata/PRJEB70895_runs.tsv",sep="\t",dtype=str).fillna("")
rna=ce[ce.library_strategy=="RNA-Seq"].copy()
rna["sample_alias"]=rna.sample_alias.str.strip()

cha=am[am.cohort=="Charolais_PRJEB70895"]
new=[]; missing=[]
for _,a in cha.iterrows():
    yr=a.batch.replace("batch_","")
    hits=rna[(rna.sample_alias==a.animal_id)|(rna.sample_alias==f"{a.animal_id}-{yr}")]
    if hits.empty: missing.append(a.animal_id); continue
    for _,e in hits.iterrows():
        new.append(dict(cohort=a.cohort,animal_id=a.animal_id,layer="mRNA",
                        run_accession=e.run_accession,sample_accession=e.sample_accession,
                        ena_alias=e.sample_alias))
print(f"Charolais animals: {len(cha)}   mRNA runs attached: {len(new)}   animals lacking mRNA: {missing or 'none'}")
used={n["sample_alias"] if "sample_alias" in n else n["ena_alias"] for n in new}
print("RNA-Seq aliases not attached to any animal:", sorted(set(rna.sample_alias)-used))

R=pd.concat([R[~((R.cohort=="Charolais_PRJEB70895")&(R.layer=="mRNA"))],pd.DataFrame(new)],ignore_index=True)
R.to_csv(f"{H}/runs.tsv",sep="\t",index=False)
with open(f"{H}/CHECKSUMS.md5","w") as fh:
    for f in ("animal_metadata.tsv","runs.tsv"):
        fh.write(f"{hashlib.md5(open(f'{H}/{f}','rb').read()).hexdigest()}  {f}\n")
print("\nruns.tsv:"); print(R.groupby(["cohort","layer"]).size().to_string())
print(open(f"{H}/CHECKSUMS.md5").read())
PY
