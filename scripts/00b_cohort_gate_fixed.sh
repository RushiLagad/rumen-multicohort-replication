#!/usr/bin/env bash
# =============================================================================
# 00b_cohort_gate_fixed.sh
# Fixes two bugs in 00_cohort_gate.sh:
#   - sorted(set(x)) - {""} raised TypeError, killing the summary table
#   - keyword matcher used substring matching, so "package" and "lineage"
#     matched "age". Now uses explicit token matching.
#
# Adds: per-cohort run counts, download size, unique sample count, and a
# scan of sample aliases/titles for animal identifiers, since a joinable ID
# is what decides whether a supplementary phenotype table can be attached.
# =============================================================================

set -uo pipefail

ROOT=/data/rumen_multicohort
mkdir -p "$ROOT"/{metadata,logs}
cd "$ROOT/metadata"

CURL="curl -sL --ssl-allow-beast --tlsv1.2 --max-time 120"

ACCESSIONS=(PRJNA448333 PRJEB10338 PRJEB31266 PRJNA987743 PRJEB70895)

FIELDS="run_accession,sample_accession,experiment_accession,sample_alias,sample_title,\
library_strategy,library_source,library_layout,instrument_platform,read_count,\
base_count,fastq_bytes,scientific_name,collection_date,country"

SUMMARY="$ROOT/metadata/cohort_summary.tsv"
printf "accession\truns\tsamples\tstrategy\tsource\tplatform\tGB\tspecies\n" > "$SUMMARY"

echo "=== Cohort gate, fixed $(date) ==="
echo

for A in "${ACCESSIONS[@]}"; do
  echo "----------------------------------------------------------------"
  echo "[$A]"
  OUT="${A}_runs.tsv"
  $CURL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${A}&result=read_run&fields=${FIELDS}&format=tsv" -o "$OUT"

  if [ ! -s "$OUT" ] || [ "$(wc -l < "$OUT")" -lt 2 ]; then
    echo "  NOT FOUND on ENA"
    printf "%s\t0\t0\t-\t-\t-\t-\t-\n" "$A" >> "$SUMMARY"
    continue
  fi

  python3 - "$OUT" "$A" "$SUMMARY" << 'PY'
import sys, pandas as pd, re
f, acc, summ = sys.argv[1:4]
d = pd.read_csv(f, sep="\t", dtype=str).fillna("")

def uniq(col):
    if col not in d: return "-"
    vals = sorted(set(d[col]) - {""})          # <-- the fix
    return ",".join(vals)[:60] or "-"

n  = len(d)
ns = d.sample_accession.nunique() if "sample_accession" in d else 0

tot = 0
if "fastq_bytes" in d:
    for v in d.fastq_bytes:
        for p in str(v).split(";"):
            if p.strip().isdigit(): tot += int(p)
size = round(tot/1e9, 1)

print(f"  runs: {n}   unique samples: {ns}   runs/sample: {round(n/max(ns,1),1)}")
print(f"  strategy: {uniq('library_strategy')}   source: {uniq('library_source')}")
print(f"  platform: {uniq('instrument_platform')}   species: {uniq('scientific_name')}")
print(f"  download: {size} GB")
if "read_count" in d:
    rc = pd.to_numeric(d.read_count, errors="coerce").dropna()
    if len(rc): print(f"  reads per run: median {int(rc.median()):,}  min {int(rc.min()):,}")

# does anything look like a per-animal identifier
print("  --- animal ID candidates ---")
found = False
for col in ("sample_alias", "sample_title"):
    if col in d:
        ex = [x for x in d[col].unique() if x][:5]
        if ex:
            print(f"    {col}: {ex}")
            found = True
if not found:
    print("    none in alias or title")

with open(summ, "a") as fh:
    fh.write(f"{acc}\t{n}\t{ns}\t{uniq('library_strategy')}\t{uniq('library_source')}"
             f"\t{uniq('instrument_platform')}\t{size}\t{uniq('scientific_name')}\n")
PY
  echo
done

echo "================================================================"
column -t -s$'\t' "$SUMMARY"
echo
echo "Total download across cohorts:"
awk -F'\t' 'NR>1 {s+=$7} END {printf "  %.1f GB\n", s}' "$SUMMARY"
