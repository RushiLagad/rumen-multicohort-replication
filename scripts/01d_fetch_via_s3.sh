#!/usr/bin/env bash
# PRJNA987743 via the NCBI SRA mirror on AWS Open Data.
# Measured: S3 111 GB/hr vs ENA 4 GB/hr in the same 25 seconds, same machine.
# Each worker: fetch .sra -> fasterq-dump -> compress -> delete .sra.
set -uo pipefail
ROOT=/data/rumen_multicohort/PRJNA987743
FQ="$ROOT/fastq"; SRA="$ROOT/sra"; TMP="$ROOT/sra_tmp"; LOG="$ROOT/logs"
mkdir -p "$FQ" "$SRA" "$TMP" "$LOG"
JOBS="${JOBS:-4}"; TPJ="${TPJ:-6}"

mapfile -t RUNS < <(grep -oE 'SRR[0-9]+' "$ROOT/metadata/rumen_urls.txt" | sort -u)
echo "=== PRJNA987743 via AWS Open Data ==="
echo "  runs: ${#RUNS[@]}   already complete: $(ls -1 "$FQ"/*_1.fastq.gz 2>/dev/null | wc -l)"
echo "  parallel: $JOBS x $TPJ   free: $(df -h /data | awk 'NR==2{print $4}')"
[ "${1:-run}" = "check" ] && exit 0

export FQ SRA TMP LOG TPJ ROOT
fetch_one () {
  acc="$1"; o1="$FQ/${acc}_1.fastq.gz"; o2="$FQ/${acc}_2.fastq.gz"
  if [ -s "$o1" ] && [ -s "$o2" ] && gzip -t "$o1" 2>/dev/null && gzip -t "$o2" 2>/dev/null; then
    echo "  skip $acc"; return 0; fi
  f="$SRA/${acc}.sra"
  curl -sS -C - --retry 6 --retry-delay 5 --retry-all-errors \
       -o "$f" "https://sra-pub-run-odp.s3.amazonaws.com/sra/${acc}/${acc}" \
    || { echo "  DOWNLOAD FAILED $acc"; return 1; }
  fasterq-dump --split-files --threads "$TPJ" --temp "$TMP" -O "$TMP" "$f" \
    > "$LOG/${acc}.fasterq.log" 2>&1 \
    || { echo "  FASTERQ FAILED $acc"; rm -f "$f" "$TMP/${acc}"*.fastq; return 1; }
  if command -v pigz >/dev/null; then
    pigz -p "$TPJ" -f "$TMP/${acc}_1.fastq" "$TMP/${acc}_2.fastq"
  else
    gzip -f "$TMP/${acc}_1.fastq" "$TMP/${acc}_2.fastq"
  fi
  mv -f "$TMP/${acc}_1.fastq.gz" "$o1" && mv -f "$TMP/${acc}_2.fastq.gz" "$o2" \
    || { echo "  MOVE FAILED $acc"; return 1; }
  rm -f "$f"
  n=$(( $(zcat "$o1" | wc -l) / 4 ))
  echo "  done $acc  $n pairs  $(( ($(stat -c%s "$o1")+$(stat -c%s "$o2"))/1000000 )) MB"
  echo -e "${acc}\t${n}" >> "$ROOT/metadata/s3_read_counts.tsv"
}
export -f fetch_one
printf '%s\n' "${RUNS[@]}" | parallel -j "$JOBS" --line-buffer fetch_one {}

echo; echo "=== complete: $(ls -1 "$FQ"/*_1.fastq.gz 2>/dev/null | wc -l) of ${#RUNS[@]} ==="
echo "  size: $(( $(du -sb "$FQ" | cut -f1)/1000000000 )) GB"
rmdir "$TMP" "$SRA" 2>/dev/null
[ -s "$ROOT/metadata/s3_read_counts.tsv" ] && python3 -I - \
  "$ROOT/metadata/s3_read_counts.tsv" "$ROOT/metadata/rumen_runs.tsv" << 'PYEOF'
import sys, csv
s3={r[0]:int(r[1]) for r in csv.reader(open(sys.argv[1]),delimiter="\t") if len(r)==2}
ena={}
for r in csv.DictReader(open(sys.argv[2]),delimiter="\t"):
    try: ena[r["run_accession"]]=int(float(r["read_count"]))
    except Exception: pass
bad=[(a,s3[a],ena[a]) for a in s3 if a in ena and abs(s3[a]-ena[a])/max(ena[a],1)>0.01]
print(f"  read-count check vs ENA: {len(s3)} compared, {len(bad)} differ by >1%")
for a,x,y in bad[:10]: print(f"    {a}: S3 {x:,} vs ENA {y:,}")
if not bad and s3: print("  all runs match ENA's recorded read counts")
PYEOF
