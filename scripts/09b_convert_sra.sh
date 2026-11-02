#!/usr/bin/env bash
# SRA -> paired gzipped FASTQ. CPU and I/O heavy: run only when preprocessing
# is finished. Deletes each .sra once its pair is written.
#   bash 09b_convert_sra.sh {charolais|li_metagenome} [check]
set -uo pipefail
COHORT="${1:-}"
case "$COHORT" in
  charolais)     SRA=/data/rumen_multicohort/PRJEB70895/sra; FQ=/data/rumen_multicohort/PRJEB70895/fastq ;;
  li_metagenome) SRA=/data/rumen448333/sra;                   FQ=/data/rumen448333/fastq ;;
  *) echo "usage: $0 {charolais|li_metagenome} [check]"; exit 1 ;;
esac
TMP="${FQ}_tmp"; LOG="${FQ}_logs"
mkdir -p "$FQ" "$TMP" "$LOG"
JOBS="${JOBS:-4}"; TPJ="${TPJ:-6}"
mapfile -t ACC < <(ls -1 "$SRA"/*.sra 2>/dev/null | xargs -rn1 basename | sed 's/.sra$//')
echo "=== $COHORT: ${#ACC[@]} sra -> $FQ ==="
echo "  already converted: $(ls -1 "$FQ"/*_1.fastq.gz 2>/dev/null | wc -l)"
echo "  free: $(df -h /data | awk 'NR==2{print $4}')"
[ "${2:-run}" = "check" ] && exit 0

export SRA FQ TMP LOG TPJ
conv () {
  a="$1"; o1="$FQ/${a}_1.fastq.gz"; o2="$FQ/${a}_2.fastq.gz"
  [ -s "$o1" ] && { echo "  skip $a"; return 0; }
  fasterq-dump --split-files --threads "$TPJ" --temp "$TMP" -O "$TMP" "$SRA/${a}.sra" \
    > "$LOG/${a}.log" 2>&1 || { echo "  FAILED $a"; rm -f "$TMP/${a}"*.fastq; return 1; }
  # single-end runs produce only _1; handle both
  if [ -s "$TMP/${a}_2.fastq" ]; then
    pigz -p "$TPJ" -f "$TMP/${a}_1.fastq" "$TMP/${a}_2.fastq"
    mv -f "$TMP/${a}_1.fastq.gz" "$o1"; mv -f "$TMP/${a}_2.fastq.gz" "$o2"
    echo "  done $a paired"
  else
    pigz -p "$TPJ" -f "$TMP/${a}_1.fastq"; mv -f "$TMP/${a}_1.fastq.gz" "$o1"
    echo "  done $a SINGLE-END"
  fi
  rm -f "$SRA/${a}.sra"
}
export -f conv
printf '%s\n' "${ACC[@]}" | parallel -j "$JOBS" --line-buffer conv {}
echo "=== $(ls -1 "$FQ"/*_1.fastq.gz 2>/dev/null | wc -l) of ${#ACC[@]}, $(( $(du -sb "$FQ" | cut -f1)/1000000000 )) GB ==="
echo "  single-end runs: $(grep -c 'SINGLE-END' <<< "$(cat $LOG/../*.log 2>/dev/null)" 2>/dev/null || echo see_log)"
rmdir "$TMP" "$SRA" 2>/dev/null
