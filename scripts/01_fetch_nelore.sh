#!/usr/bin/env bash
# =============================================================================
# 01_fetch_nelore.sh
# PRJNA987743, Nelore cohort for the multi-cohort rumen replication study
#
# 52 Bos indicus steers, two diets, RFI and residual methane emission.
# 104 metagenomes: 52 rumen + 52 faecal. We want RUMEN ONLY.
#
# Stage 1 classifies every run as rumen or faecal from its sample title and
# alias, and STOPS so the split can be checked before 140 GB is pulled.
# Stage 2 downloads, and only runs when you pass "go" as the first argument.
#
# Downloads are resumable: rerun the same command and curl -C - picks up
# where it stopped.
# =============================================================================

set -uo pipefail

ROOT=/data/rumen_multicohort/PRJNA987743
mkdir -p "$ROOT"/{metadata,fastq,logs}
cd "$ROOT/metadata"

CURL="curl -sL --ssl-allow-beast --tlsv1.2 --max-time 180"
GO="${1:-check}"
THREADS=8          # ENA throttles above this, 16 failed last time

FIELDS="run_accession,sample_accession,sample_alias,sample_title,fastq_ftp,\
fastq_bytes,read_count,library_strategy,instrument_platform"

echo "=== PRJNA987743 Nelore $(date) ==="
echo "mode: $GO    (pass 'go' to start the download)"
echo

if [ ! -s runs_full.tsv ]; then
  $CURL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA987743&result=read_run&fields=${FIELDS}&format=tsv" -o runs_full.tsv
fi
echo "run table: $(( $(wc -l < runs_full.tsv) - 1 )) runs"

python3 - << 'PY'
import pandas as pd, re
d = pd.read_csv("runs_full.tsv", sep="\t", dtype=str).fillna("")

def tissue(row):
    s = (row.sample_title + " " + row.sample_alias).lower()
    if "fece" in s or "faec" in s or re.search(r"fz\b", s): return "feces"
    if "rumen" in s or re.search(r"rz\b|_r\b", s):          return "rumen"
    return "unknown"

d["tissue"] = d.apply(tissue, axis=1)
d["animal"] = d.sample_alias.str.extract(r"^(\d+)")[0]

print("\n--- tissue classification ---")
print(d.tissue.value_counts().to_string())
print(f"\nunique animal IDs: {d.animal.nunique()}")

for t in d.tissue.unique():
    sub = d[d.tissue == t]
    gb = sum(int(p) for v in sub.fastq_bytes for p in str(v).split(";") if p.strip().isdigit())/1e9
    print(f"  {t:8s}  runs {len(sub):3d}  animals {sub.animal.nunique():3d}  {gb:6.1f} GB")
    print(f"            examples: {list(sub.sample_title.head(3))}")

if (d.tissue == "unknown").any():
    print("\n  UNCLASSIFIED, inspect these before downloading:")
    print(d[d.tissue=="unknown"][["run_accession","sample_alias","sample_title"]].head(10).to_string(index=False))

rum = d[d.tissue == "rumen"]
urls = []
for v in rum.fastq_ftp:
    urls += [u for u in str(v).split(";") if u.strip()]
open("rumen_urls.txt","w").write("\n".join(urls) + "\n")
rum.to_csv("rumen_runs.tsv", sep="\t", index=False)
print(f"\nwrote rumen_urls.txt  ({len(urls)} files)")
print(f"wrote rumen_runs.tsv  ({len(rum)} runs, {rum.animal.nunique()} animals)")
PY

if [ "$GO" != "go" ]; then
  echo
  echo "================================================================"
  echo "STOPPED before downloading. Check the tissue split above."
  echo "If rumen = 52 runs / 52 animals, start with:"
  echo "  bash $0 go"
  exit 0
fi

echo
echo "=== downloading rumen samples ==="
cd "$ROOT/fastq"
N=$(wc -l < ../metadata/rumen_urls.txt)
echo "$N files, $THREADS parallel, resumable"

cat ../metadata/rumen_urls.txt | parallel -j $THREADS --retries 3 \
  'curl -sL -C - --ssl-allow-beast --tlsv1.2 --max-time 1800 -O ftp://{}' \
  2>&1 | tail -5

echo
echo "downloaded: $(ls *.fastq.gz 2>/dev/null | wc -l) of $N"
du -sh .
echo
echo "integrity check"
ls *.fastq.gz 2>/dev/null | parallel -j 16 'gzip -t {} 2>/dev/null || echo "CORRUPT {}"' | grep CORRUPT || echo "  all files intact"
echo
echo "If any are missing, rerun 'bash $0 go' — curl resumes partial files."
