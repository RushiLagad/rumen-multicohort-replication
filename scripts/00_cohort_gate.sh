#!/usr/bin/env bash
# =============================================================================
# 00_cohort_gate.sh
# Multi-cohort rumen replication study, candidate screening
#
# Decides whether the study is viable BEFORE any sequence download.
# Three questions per accession:
#   1 does it exist and how many runs
#   2 what is the library strategy and total download size
#   3 do per-animal phenotype fields exist, or only group labels
#
# Question 3 is the one that kills projects. Group labels alone mean you
# cannot model RFI as continuous, and the whole design collapses.
#
# Nothing here downloads reads. Metadata only, a few hundred KB.
# =============================================================================

set -uo pipefail

ROOT=/data/rumen_multicohort
mkdir -p "$ROOT"/{metadata,logs}
cd "$ROOT/metadata"

# curl flags that work around the OpenSSL handshake failure on this machine
CURL="curl -sL --ssl-allow-beast --tlsv1.2 --max-time 120"

ACCESSIONS=(
  PRJNA448333    # baseline, Li et al. 2019, 48 animals, DNA + total RNA + mRNA
  PRJEB10338     # Auffret, Scotland
  PRJEB31266     # Auffret 2020, RFI / FCR / methane
  PRJNA987743    # Nelore, residual methane
  PRJNA1151032   # Angus bulls, RFI-EPD
  PRJEB70895     # Charolais, 16S + mRNA, two diets two years
)

FIELDS="run_accession,sample_accession,secondary_sample_accession,experiment_accession,\
study_accession,sample_alias,sample_title,sample_description,library_strategy,\
library_source,library_layout,instrument_platform,read_count,base_count,fastq_bytes,\
scientific_name,tax_id,collection_date,country,first_public"

echo "=== Cohort gate $(date) ==="
echo

SUMMARY="$ROOT/metadata/cohort_summary.tsv"
printf "accession\texists\truns\tstrategy\tsource\tplatform\tGB\tspecies\n" > "$SUMMARY"

for A in "${ACCESSIONS[@]}"; do
  echo "----------------------------------------------------------------"
  echo "[$A]"
  OUT="${A}_runs.tsv"

  $CURL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${A}&result=read_run&fields=${FIELDS}&format=tsv" -o "$OUT"

  if [ ! -s "$OUT" ] || [ "$(wc -l < "$OUT")" -lt 2 ]; then
    echo "  NOT FOUND or empty on ENA"
    printf "%s\tNO\t0\t-\t-\t-\t-\t-\n" "$A" >> "$SUMMARY"
    rm -f "$OUT"
    continue
  fi

  python3 - "$OUT" "$A" "$SUMMARY" << 'PY'
import sys, pandas as pd
f, acc, summ = sys.argv[1], sys.argv[2], sys.argv[3]
d = pd.read_csv(f, sep="\t", dtype=str).fillna("")

n = len(d)
strat = ",".join(sorted(set(d.get("library_strategy", pd.Series(dtype=str)))) - {""}) or "-"
src   = ",".join(sorted(set(d.get("library_source",   pd.Series(dtype=str)))) - {""}) or "-"
plat  = ",".join(sorted(set(d.get("instrument_platform", pd.Series(dtype=str)))) - {""}) or "-"
spec  = ",".join(sorted(set(d.get("scientific_name",  pd.Series(dtype=str)))) - {""}) or "-"

def gb(s):
    tot = 0
    for v in s:
        for part in str(v).split(";"):
            if part.strip().isdigit(): tot += int(part)
    return round(tot/1e9, 1)
size = gb(d["fastq_bytes"]) if "fastq_bytes" in d else 0.0

print(f"  runs: {n}   strategy: {strat}   source: {src}")
print(f"  platform: {plat}   species: {spec}")
print(f"  download size: {size} GB")

# unique biological samples vs runs
if "sample_accession" in d:
    print(f"  unique samples: {d.sample_accession.nunique()}  (runs per sample: "
          f"{round(n / max(d.sample_accession.nunique(), 1), 1)})")

# what the sample names look like, since these carry group labels in many studies
for col in ("sample_alias", "sample_title", "sample_description"):
    if col in d and d[col].str.len().max() > 0:
        ex = [x for x in d[col].unique() if x][:4]
        print(f"  {col}: {ex}")

with open(summ, "a") as fh:
    fh.write(f"{acc}\tYES\t{n}\t{strat}\t{src}\t{plat}\t{size}\t{spec}\n")
PY

  # ---- the phenotype check -----------------------------------------------
  # Sample attributes live in the BioSample XML, not the run report. Pull a few
  # and list every attribute tag so we can see whether RFI, DMI, ADG, FCR or
  # methane are recorded per animal.
  SAMP=$(awk -F'\t' 'NR>1 && $2!="" {print $2}' "$OUT" | sort -u | head -3 | tr '\n' ',' | sed 's/,$//')
  if [ -n "$SAMP" ]; then
    echo "  --- sample attributes (first 3 samples) ---"
    $CURL "https://www.ebi.ac.uk/ena/browser/api/xml/${SAMP}" -o "${A}_samples.xml"
    python3 - "${A}_samples.xml" << 'PY'
import sys, re
try:
    x = open(sys.argv[1]).read()
except Exception:
    print("    (no XML returned)"); raise SystemExit
tags = re.findall(r"<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>", x, re.S)
if not tags:
    print("    no attributes found")
else:
    seen = {}
    for t, v in tags:
        seen.setdefault(t.strip(), set()).add(v.strip()[:40])
    KEY = ("rfi","residual","feed","intake","dmi","adg","gain","fcr","conversion",
           "methane","ch4","efficiency","breed","diet","weight","sex","age")
    for t in sorted(seen):
        hit = "  <== PHENOTYPE" if any(k in t.lower() for k in KEY) else ""
        print(f"    {t}: {sorted(seen[t])[:3]}{hit}")
PY
  fi
  echo
done

echo "================================================================"
column -t -s$'\t' "$SUMMARY"
echo
echo "VERDICT GUIDE"
echo "  A cohort is usable only if per-animal RFI (or methane) is recorded"
echo "  somewhere: in the BioSample attributes above, in the sample alias, or"
echo "  in a supplementary table from the paper that can be joined by animal ID."
echo "  Group labels alone (high/low) mean continuous RFI modelling is out and"
echo "  the cohort drops to a binary contrast at best."
echo
echo "Next: for any cohort with no phenotype in BioSample, check the paper's"
echo "supplementary tables before committing to the download."
