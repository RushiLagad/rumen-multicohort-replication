#!/usr/bin/env bash
# =============================================================================
# 07_preprocess.sh   fastp -> minimap2 host depletion, manifest-driven
#   bash 07_preprocess.sh nelore [check]
# Resumable. fastp intermediates deleted as soon as host removal succeeds.
# =============================================================================
set -uo pipefail
ROOT=/data/rumen_multicohort; REF=/data/refs/host
COHORT="${1:-}"; MODE="${2:-run}"
[ -n "$COHORT" ] || { echo "usage: $0 {nelore|charolais_mrna|li} [check]"; exit 1; }
case "$COHORT" in
  nelore)         KEY="Nelore_PRJNA987743";   FQ="$ROOT/PRJNA987743/fastq"; HOST="${HOST_REF_NELORE:-ARS-UCD2.0}";;
  charolais_mrna) KEY="Charolais_PRJEB70895"; FQ="$ROOT/PRJEB70895/fastq";  HOST="ARS-UCD2.0";;
  li)             KEY="Li_PRJNA448333";       FQ="/data/rumen448333/fastq"; HOST="ARS-UCD2.0";;
  *) echo "unknown cohort"; exit 1;;
esac
OUT="$ROOT/processed/$COHORT"; TMP="$OUT/tmp"; LOG="$OUT/logs"
mkdir -p "$OUT" "$TMP" "$LOG"
JOBS="${JOBS:-3}"; TPJ="${TPJ:-8}"

ACCF="$ROOT/processed/${COHORT}.acc"
if [ -s "$ACCF" ]; then mapfile -t RUNS < "$ACCF"
else mapfile -t RUNS < <(awk -F'\t' -v k="$KEY" '$1==k{print $4}' "$ROOT/harmonized/runs.tsv" | sort -u); fi
echo "=== $COHORT ==="
echo "  manifest : ${#RUNS[@]} runs"
echo "  done     : $(ls -1 "$OUT"/*_clean_1.fq.gz 2>/dev/null | wc -l)"
echo "  host ref : $HOST"
echo "  parallel : $JOBS x $TPJ threads"
[ "$MODE" = "check" ] && exit 0
[ -s "$REF/${HOST}.mmi" ] || { echo "MISSING $REF/${HOST}.mmi"; exit 1; }

export FQ OUT TMP LOG REF HOST TPJ
process_one () {
  r="$1"; o1="$OUT/${r}_clean_1.fq.gz"; o2="$OUT/${r}_clean_2.fq.gz"
  [ -s "$o1" ] && [ -s "$o2" ] && { echo "  skip $r"; return 0; }
  R1="$FQ/${r}_1.fastq.gz"; R2="$FQ/${r}_2.fastq.gz"
  gzip -t "$R1" 2>/dev/null && gzip -t "$R2" 2>/dev/null || { echo "  not ready $r"; return 1; }
  t1="$TMP/${r}_t1.fq.gz"; t2="$TMP/${r}_t2.fq.gz"
  fastp -i "$R1" -I "$R2" -o "$t1" -O "$t2" \
        --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
        --length_required 50 --thread 4 \
        --json "$LOG/${r}.fastp.json" --html "$LOG/${r}.fastp.html" \
        > "$LOG/${r}.fastp.log" 2>&1 || { echo "  FASTP FAILED $r"; rm -f "$t1" "$t2"; return 1; }
  minimap2 -ax sr -t "$TPJ" "$REF/${HOST}.mmi" "$t1" "$t2" 2> "$LOG/${r}.mm2.log" \
    | samtools view -@ 2 -bu -f 12 -F 256 - \
    | samtools collate -@ 2 -O -u - \
    | samtools fastq -@ 2 -1 "$o1" -2 "$o2" -0 /dev/null -s /dev/null -n - \
      > "$LOG/${r}.sam.log" 2>&1 \
    || { echo "  HOST REMOVAL FAILED $r"; rm -f "$o1" "$o2" "$t1" "$t2"; return 1; }
  rm -f "$t1" "$t2"
  raw=$(python3 -c "import json;print(json.load(open('$LOG/${r}.fastp.json'))['summary']['before_filtering']['total_reads'])" 2>/dev/null || echo 0)
  aft=$(python3 -c "import json;print(json.load(open('$LOG/${r}.fastp.json'))['summary']['after_filtering']['total_reads'])" 2>/dev/null || echo 0)
  kept=$(( $(zcat "$o1" | wc -l) / 4 * 2 ))
  echo -e "$r\t$raw\t$aft\t$kept" >> "$OUT/read_accounting.tsv"
  echo "  done $r  raw $raw -> trimmed $aft -> nonhost $kept"
}
export -f process_one
printf '%s\n' "${RUNS[@]}" | parallel -j "$JOBS" --line-buffer process_one {}

echo; echo "=== $COHORT: $(ls -1 "$OUT"/*_clean_1.fq.gz 2>/dev/null | wc -l) of ${#RUNS[@]} ==="
[ -s "$OUT/read_accounting.tsv" ] && python3 -I - "$OUT/read_accounting.tsv" << 'PY'
import sys,csv,statistics as st
seen={}
for r in csv.reader(open(sys.argv[1]),delimiter="\t"):
    if len(r)==4: seen[r[0]]=r
rows=list(seen.values())
kept=[int(r[3]) for r in rows]; aft=[int(r[2]) for r in rows]
host=[100*(1-k/max(a,1)) for a,k in zip(aft,kept)]
print(f"  samples {len(rows)}  non-host median {st.median(kept):,.0f}")
print(f"  host fraction median {st.median(host):.1f}%  range {min(host):.1f}-{max(host):.1f}%")
o=[r[0] for r,h in zip(rows,host) if h>60]
if o: print(f"  OUTLIERS to inspect: {o}")
PY
