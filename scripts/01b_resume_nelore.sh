#!/usr/bin/env bash
set -uo pipefail
ROOT=/data/rumen_multicohort/PRJNA987743
FQ="$ROOT/fastq"; URLS="$ROOT/metadata/rumen_urls.txt"
mkdir -p "$FQ" "$ROOT/logs"; cd "$FQ"
MODE="${1:-resume}"
THREADS="${THREADS:-6}"
total=$(wc -l < "$URLS")

report () {
  echo "  files : $(ls -1 *.fastq.gz 2>/dev/null | wc -l) of $total"
  echo "  size  : $(( $(du -sb . | cut -f1) / 1000000000 )) GB of ~143 GB"
}

verify () {
  echo "=== integrity check ==="
  : > "$ROOT/metadata/incomplete.txt"
  ls -1 *.fastq.gz 2>/dev/null | parallel -j 16 \
    'gzip -t {} 2>/dev/null || echo {}' >> "$ROOT/metadata/incomplete.txt"
  bad=$(grep -c . "$ROOT/metadata/incomplete.txt" || true)
  echo "  incomplete: $bad  (partial, not corrupt - they resume next run)"
}

echo "=== PRJNA987743 resume $(date) ==="
report; echo
[ "$MODE" = "check" ] && { verify; exit 0; }

echo "=== resuming, $THREADS streams, no total-time limit ==="
cat "$URLS" | parallel -j "$THREADS" --retries 2 \
  'curl -sS -C - --ssl-allow-beast --tlsv1.2 \
        --speed-limit 10000 --speed-time 120 \
        --retry 5 --retry-delay 15 --retry-all-errors \
        -O ftp://{}' 2>&1 | tail -10
echo; report; echo; verify
echo
echo "Still incomplete? Run again. Every file resumes from its last byte."
