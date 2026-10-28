#!/usr/bin/env bash
# =============================================================================
# 02_fetch_phenotypes.sh
# Pulls the per-animal phenotype tables for both new cohorts and dumps every
# sheet name and column header, so the harmonization join keys are known
# facts rather than guesses.
#
#   Nelore  PRJNA987743 : Frontiers Microbiol 2023, PMC10754987
#                         Table_2.XLSX = Supplementary Table S2
#                         (nine XLSX supplements, S1-S9, in file order)
#                         Table_1.XLSX = S1, pulled too as a cross-check
#   Charolais PRJEB70895: Dataverse doi:10.57745/ACNRVQ, mapping.tab
#
# Metadata only. A few hundred KB. Safe to run while reads are downloading.
# =============================================================================
set -uo pipefail

NEL=/data/rumen_multicohort/PRJNA987743/metadata
CHA=/data/rumen_multicohort/PRJEB70895/metadata
mkdir -p "$NEL" "$CHA"

UA='Mozilla/5.0 (X11; Linux x86_64)'
CURL="curl -sSL --ssl-allow-beast --tlsv1.2 --retry 4 --retry-delay 5 -A $UA"

echo "=== Nelore supplementary tables ==="
for n in 1 2; do
  f="$NEL/Supp_Table_S${n}.xlsx"
  [ -s "$f" ] || $CURL -o "$f" \
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/10754987/bin/Table_${n}.XLSX"
  echo "  S${n}: $(stat -c%s "$f" 2>/dev/null || echo 0) bytes  $(file -b "$f" | cut -c1-40)"
done

echo
echo "=== Charolais mapping.tab ==="
f="$CHA/mapping.tab"
[ -s "$f" ] || $CURL -o "$f" \
  "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/ACNRVQ"
echo "  $(stat -c%s "$f" 2>/dev/null || echo 0) bytes  $(file -b "$f" | cut -c1-40)"

echo
echo "================================================================"
echo "STRUCTURE"
echo "================================================================"

python3 - << 'PY'
import pandas as pd, glob, os
pd.set_option("display.width", 200)

for f in sorted(glob.glob("/data/rumen_multicohort/PRJNA987743/metadata/Supp_Table_S*.xlsx")):
    print("\n" + "="*64)
    print(os.path.basename(f))
    print("="*64)
    try:
        xl = pd.ExcelFile(f)
    except Exception as e:
        print("  UNREADABLE:", e); continue
    print("  sheets:", xl.sheet_names)
    for sh in xl.sheet_names:
        # header row is sometimes below a title row, so try the first three
        best, hdr = None, 0
        for h in (0, 1, 2):
            d = xl.parse(sh, header=h, nrows=60)
            named = sum(1 for c in d.columns if not str(c).startswith("Unnamed"))
            if best is None or named > best[0]:
                best, hdr = (named, d), h
        d = best[1]
        full = xl.parse(sh, header=hdr)
        print(f"\n  --- sheet '{sh}'  header_row={hdr}  {full.shape[0]} rows x {full.shape[1]} cols")
        for c in full.columns:
            s = full[c].dropna()
            num = pd.to_numeric(s, errors="coerce").dropna()
            if len(num) >= max(3, 0.6*len(s)) and len(s):
                print(f"      {str(c)[:38]:40s} NUM   n={len(num):3d}  "
                      f"range {num.min():.3f} .. {num.max():.3f}  mean {num.mean():.3f}")
            else:
                ex = [str(x)[:16] for x in s.unique()[:4]]
                print(f"      {str(c)[:38]:40s} TXT   n={len(s):3d}  uniq={s.nunique():3d}  {ex}")

f = "/data/rumen_multicohort/PRJEB70895/metadata/mapping.tab"
if os.path.exists(f) and os.path.getsize(f) > 0:
    print("\n" + "="*64)
    print("Charolais mapping.tab")
    print("="*64)
    d = pd.read_csv(f, sep="\t", dtype=str).fillna("")
    print(f"  {d.shape[0]} rows x {d.shape[1]} cols")
    for c in d.columns:
        s = d[c][d[c] != ""]
        print(f"      {str(c)[:38]:40s} n={len(s):3d}  uniq={s.nunique():3d}  "
              f"{[str(x)[:14] for x in s.unique()[:4]]}")
PY

echo
echo "Paste the block above. That is what the harmonization script joins on."
