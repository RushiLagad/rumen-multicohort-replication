#!/usr/bin/env bash
# =============================================================================
# 01c_download_supervised.sh
# Fixes what killed the overnight run: --speed-limit 10000 --speed-time 120
# aborted transfers that merely slowed down when the Kraken2 database competed
# for bandwidth, and parallel --retries 2 then gave up and exited.
#
#   - speed floor relaxed to 1 KB/s over 10 minutes (stalled, not merely slow)
#   - supervisor loop reruns passes until nothing is incomplete
#   - PROTO=https by default; ENA serves the same files over both
# =============================================================================
set -uo pipefail
ROOT=/data/rumen_multicohort/PRJNA987743
FQ="$ROOT/fastq"; URLS="$ROOT/metadata/rumen_urls.txt"
mkdir -p "$FQ" "$ROOT/logs"; cd "$FQ"
PROTO="${PROTO:-https}"; THREADS="${THREADS:-8}"; MAXPASS="${MAXPASS:-40}"
total=$(wc -l < "$URLS")

incomplete () { ls -1 *.fastq.gz 2>/dev/null | parallel -j 16 'gzip -t {} 2>/dev/null || echo {}' | grep -c . ; }

for pass in $(seq 1 "$MAXPASS"); do
  bad=$(incomplete)
  gb=$(( $(du -sb . | cut -f1) / 1000000000 ))
  echo "[$(date +%H:%M:%S)] pass $pass  ${gb} GB  incomplete ${bad}/${total}"
  [ "$bad" -eq 0 ] && { echo "ALL COMPLETE"; break; }

  cat "$URLS" | parallel -j "$THREADS" \
    "curl -sS -C - --retry 10 --retry-delay 10 --retry-all-errors \
          --speed-limit 1000 --speed-time 600 -O ${PROTO}://{}" 2>&1 | tail -3

  sleep 10
done
echo "[$(date +%H:%M:%S)] final: $(( $(du -sb . | cut -f1)/1000000000 )) GB, incomplete $(incomplete)"
