#!/usr/bin/env bash
# =============================================================================
# 02b_fetch_phenotypes.sh
# Per-animal phenotype tables for both new cohorts, then a full structure dump.
#
# Why this replaces 02_fetch_phenotypes.sh:
#   1. -A $UA was unquoted, so curl parsed the user-agent as extra URLs and
#      -o bound to the first garbage one. Every real payload went to stdout.
#      No user-agent is needed here at all, so it is gone.
#   2. pmc.ncbi.nlm.nih.gov/articles/instance/.../bin/*.XLSX now returns a
#      JavaScript proof-of-work page, not the file. Unfixable from curl.
#      Europe PMC serves the identical supplement set as one zip.
#
# Metadata only, well under 5 MB. Safe while reads are downloading.
# =============================================================================
set -uo pipefail

NEL=/data/rumen_multicohort/PRJNA987743/metadata
CHA=/data/rumen_multicohort/PRJEB70895/metadata
mkdir -p "$NEL" "$CHA" /data/rumen_multicohort/logs

ok () { [ -s "$1" ] && ! head -c 200 "$1" | grep -qi "<html"; }

# ---------------------------------------------------------------- Nelore ----
echo "=== Nelore supplements, PMC10754987 ==="
Z="$NEL/supplementary.zip"

if ! ok "$Z"; then
  echo "  source 1: Europe PMC"
  curl -sSL --retry 4 --retry-delay 5 -o "$Z" \
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10754987/supplementaryFiles"
fi

if ok "$Z" && file -b "$Z" | grep -qi zip; then
  echo "  got zip: $(stat -c%s "$Z") bytes"
  unzip -o -j -d "$NEL" "$Z" >/dev/null 2>&1
  echo "  extracted:"
  ls -1S "$NEL" | grep -iE '\.xlsx?$' | sed 's/^/    /'
else
  echo "  Europe PMC did not return a zip."
  echo "  source 2: scraping Frontiers for direct supplement links"
  curl -sSL --retry 3 -o "$NEL/frontiers.html" \
    "https://www.frontiersin.org/articles/10.3389/fmicb.2023.1282851/full"
  grep -oiE 'https?://[^"'"'"' ]+\.xlsx?' "$NEL/frontiers.html" | sort -u \
    > "$NEL/frontiers_xlsx_urls.txt"
  n=$(grep -c . "$NEL/frontiers_xlsx_urls.txt" || true)
  echo "  found $n candidate URLs"
  if [ "$n" -gt 0 ]; then
    sed 's/^/    /' "$NEL/frontiers_xlsx_urls.txt"
    i=1
    while read -r u; do
      curl -sSL --retry 3 -o "$NEL/frontiers_S${i}.xlsx" "$u" && \
        echo "    saved frontiers_S${i}.xlsx $(stat -c%s "$NEL/frontiers_S${i}.xlsx") bytes"
      i=$((i+1))
    done < "$NEL/frontiers_xlsx_urls.txt"
  else
    echo
    echo "  BOTH AUTOMATED SOURCES FAILED. Do this by hand, it is 30 seconds:"
    echo "    open https://pmc.ncbi.nlm.nih.gov/articles/PMC10754987/ in Firefox"
    echo "    scroll to Supplementary material, click the 2nd XLSX link (19.2 KB)"
    echo "    then:  mv ~/Downloads/Table_2.XLSX $NEL/Supp_Table_S2.xlsx"
    echo "  The browser solves the proof-of-work automatically. curl cannot."
  fi
fi

# ------------------------------------------------------------- Charolais ----
echo
echo "=== Charolais mapping.tab, doi:10.57745/ACNRVQ ==="
M="$CHA/mapping.tab"
if ! ok "$M"; then
  curl -sSL --retry 4 --retry-delay 5 -o "$M" \
    "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/ACNRVQ"
fi
if ok "$M"; then
  echo "  saved: $(stat -c%s "$M") bytes, $(( $(wc -l < "$M") - 1 )) rows"
  md5sum "$M" | sed 's/^/  md5 /'
else
  echo "  FAILED"
fi

# ------------------------------------------------------------- structure ----
echo
echo "================================================================"
echo "STRUCTURE"
echo "================================================================"

python3 -I - << 'PY'
import pandas as pd, glob, os
NEL="/data/rumen_multicohort/PRJNA987743/metadata"

files=[f for f in glob.glob(os.path.join(NEL,"*"))
       if f.lower().endswith((".xlsx",".xls"))]
if not files:
    print("\n  no Nelore spreadsheet on disk yet, see the manual step above")

for f in sorted(files, key=os.path.getsize):
    print("\n"+"="*64); print(os.path.basename(f), os.path.getsize(f),"bytes"); print("="*64)
    try: xl=pd.ExcelFile(f)
    except Exception as e: print("  UNREADABLE:",e); continue
    print("  sheets:", xl.sheet_names)
    for sh in xl.sheet_names:
        best=None
        for h in (0,1,2):
            try: d=xl.parse(sh,header=h,nrows=80)
            except Exception: continue
            named=sum(1 for c in d.columns if not str(c).startswith("Unnamed"))
            if best is None or named>best[0]: best=(named,h)
        if best is None: continue
        hdr=best[1]; full=xl.parse(sh,header=hdr)
        print(f"\n  --- '{sh}'  header_row={hdr}  {full.shape[0]} rows x {full.shape[1]} cols")
        for c in full.columns:
            s=full[c].dropna()
            if not len(s): print(f"      {str(c)[:38]:40s} EMPTY"); continue
            num=pd.to_numeric(s,errors="coerce").dropna()
            if len(num)>=0.6*len(s):
                print(f"      {str(c)[:38]:40s} NUM  n={len(num):3d}  "
                      f"{num.min():.3f} .. {num.max():.3f}  mean {num.mean():.3f}  sd {num.std():.3f}")
            else:
                print(f"      {str(c)[:38]:40s} TXT  n={len(s):3d}  uniq={s.nunique():3d}  "
                      f"{[str(x)[:16] for x in s.unique()[:4]]}")
PY

echo
echo "Paste the Nelore block. The Charolais side is already resolved."
