#!/usr/bin/env bash
# Download SRA files only, no conversion. Network-bound, negligible CPU, so it
# runs alongside preprocessing without contending for threads.
#   bash 09_prefetch_sra.sh {charolais|li_metagenome} [check]
set -uo pipefail
ROOT=/data/rumen_multicohort
COHORT="${1:-}"; JOBS="${JOBS:-6}"
case "$COHORT" in
  charolais)
    mapfile -t ACC < <(awk -F'\t' '$1=="Charolais_PRJEB70895"{print $4}' "$ROOT/harmonized/runs.tsv" | sort -u)
    OUT="$ROOT/PRJEB70895/sra" ;;
  li_metagenome)
    mapfile -t ACC < <(awk -F'\t' 'NR>1 && tolower($0) ~ /wgs/ {print $1}' /data/rumen448333/metadata/runs.tsv | sort -u)
    OUT="/data/rumen448333/sra" ;;
  *) echo "usage: $0 {charolais|li_metagenome} [check]"; exit 1 ;;
esac
mkdir -p "$OUT"
echo "=== $COHORT: ${#ACC[@]} accessions -> $OUT ==="
echo "  have: $(ls -1 "$OUT"/*.sra 2>/dev/null | wc -l)   free: $(df -h /data | awk 'NR==2{print $4}')"
[ "${2:-run}" = "check" ] && { printf '  %s\n' "${ACC[@]:0:5}"; exit 0; }
export OUT
get_one () {
  a="$1"; f="$OUT/${a}.sra"
  [ -s "$f" ] && { echo "  skip $a"; return 0; }
  curl -sS -C - --retry 6 --retry-delay 5 --retry-all-errors \
       -o "$f" "https://sra-pub-run-odp.s3.amazonaws.com/sra/${a}/${a}" \
    && echo "  got $a  $(( $(stat -c%s "$f")/1000000 )) MB" \
    || { echo "  FAILED $a"; rm -f "$f"; }
}
export -f get_one
printf '%s\n' "${ACC[@]}" | parallel -j "$JOBS" --line-buffer get_one {}
echo "=== $(ls -1 "$OUT"/*.sra 2>/dev/null | wc -l) of ${#ACC[@]}, $(( $(du -sb "$OUT" | cut -f1)/1000000000 )) GB ==="
