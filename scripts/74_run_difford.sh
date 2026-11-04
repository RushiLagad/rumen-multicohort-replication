#!/bin/bash
# Amendment 28 driver for Difford: bacterial and archaeal markers in parallel (16 threads each). LIMIT=n for a test.
set -euo pipefail; cd ${RUMEN_ROOT:-/data/rumen_multicohort}; L=logs/16s/difford; mkdir -p $L
grep -q "^\*\*Amendment 30 (" plan/ANALYSIS_PLAN.md && git diff --quiet plan/ANALYSIS_PLAN.md || { echo "commit Amendment 28 first"; exit 1; }
[ -s work/charolais16s/silva138_train.fa ] || { echo "SILVA training file missing"; exit 1; }
t() { echo "$(date +%F_%T) $*"; }
run() { m=$1; ins=$2; t "$m trim"; bash scripts/72_trim_difford.sh $m > $L/${m}_trim.log 2>&1; tail -n 3 $L/${m}_trim.log
  t "$m dada2"; NT=16 mamba run -p /home/rlagad/mambaforge/envs/silva_taxonomy Rscript scripts/73_dada2_difford.R $m $ins > $L/${m}_dada2.log 2>&1
  grep -E "runs with|truncLen|merge rate|FALLBACK|off-target|^table|forward reads" $L/${m}_dada2.log; t "$m DONE"; }
run bact 471 & p1=$!; run arch 509 & p2=$!
ok=1; wait $p1 || { ok=0; t "bact FAILED"; }; wait $p2 || { ok=0; t "arch FAILED"; }
[ $ok = 1 ] && { touch $L/DONE; t ALL DONE; } || { t "NOT DONE, see logs"; exit 1; }
