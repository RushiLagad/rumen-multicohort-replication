#!/usr/bin/env bash
# Database-independent work only. Nothing here depends on the Kraken2 decision.
set -uo pipefail
cd /data/rumen_multicohort
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

say "waiting for profiling and charolais conversion to finish"
while pgrep -f 08_profile >/dev/null || pgrep -f 09b_convert >/dev/null; do sleep 120; done
say "prerequisites done"

say "=== preprocess charolais mRNA (31 runs) ==="
JOBS=3 TPJ=8 bash scripts/07_preprocess.sh charolais_mrna

say "=== convert li metagenomes (48 runs) ==="
JOBS=4 TPJ=6 bash scripts/09b_convert_sra.sh li_metagenome

say "=== preprocess li metagenomes (48 runs) ==="
JOBS=3 TPJ=8 bash scripts/07_preprocess.sh li

say "ALL OVERNIGHT WORK COMPLETE"
df -h /data | tail -1
