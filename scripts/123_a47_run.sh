#!/bin/bash
# Amendment 47: late analysis, run after all primary, late and Amendment 46 results were seen. Writes only work/a46_inputs/wallace_prokaryotes_rfi_nrc.rds (never committed)
# and results/late_a47/; never overwrites. Usage (inside systemd-run, see the launch block): 123_a47_run.sh
#   0. dump of the Wallace RFI NRC table (81 with A47=1 TRAIT=rfi_nrc) for 47.4
#   1. 120 (47.1: generator check, Checks 1 and 2, selection, intervals on the real data)
#   2. 121 (47.2 labels and 47.4 four-cohort pooling), after 120
#   3. 122 (47.3 null-scan calibration, 200 outer replicates per table)
# Env: C (cores, default 16).
set -u; cd /data/rumen_multicohort; C=${C:-16}
MB=/home/rlagad/.local/share/mamba/envs; RR=$MB/rumen-r/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
L=logs/late_a47; mkdir -p $L
S="scripts/lib_a46.R scripts/lib_a47.R scripts/lib_replication.R scripts/81_estimate_methane.R scripts/120_a47_intervals.R scripts/121_a47_labels.R scripts/122_a47_null_scan.R scripts/123_a47_run.sh"
grep -q '^\*\*Amendment 47 (' plan/ANALYSIS_PLAN.md && git diff --quiet HEAD -- plan/ANALYSIS_PLAN.md && git ls-files --error-unmatch $S >/dev/null 2>&1 && git diff --quiet HEAD -- $S \
  || { echo "Amendment 47 and its scripts must be committed (and unmodified) before any run"; exit 1; }
echo "$(date '+%F %T') start | Amendment 47 commit $(git log -1 --format=%h -- plan/ANALYSIS_PLAN.md) | scripts commit $(git log -1 --format=%h -- scripts/lib_a47.R) | C $C"
A47=1 TRAIT=rfi_nrc $RR scripts/81_estimate_methane.R wallace_prokaryotes real > $L/dump_wallace_rfi_nrc.log 2>&1
[ -s work/a46_inputs/wallace_prokaryotes_rfi_nrc.rds ] || { echo "$(date '+%F %T') Wallace RFI NRC dump FAILED (see $L/dump_wallace_rfi_nrc.log)"; exit 1; }
echo "$(date '+%F %T') dump done: $(grep -o 'animals [0-9]* | features [0-9]* | df [0-9]*' $L/dump_wallace_rfi_nrc.log)"
A47_CORES=$C $RR scripts/120_a47_intervals.R > $L/120_intervals.log 2>&1 && echo "$(date '+%F %T') 120 done: $(sed -n 2p results/late_a47/intervals/chosen_method.txt)" \
  || { echo "$(date '+%F %T') 120 FAILED (see $L/120_intervals.log)"; exit 1; }
$RR scripts/121_a47_labels.R > $L/121_labels.log 2>&1 && echo "$(date '+%F %T') 121 done" || echo "$(date '+%F %T') 121 FAILED (see $L/121_labels.log)"
A47_CORES=$C $RR scripts/122_a47_null_scan.R > $L/122_null_scan.log 2>&1 && echo "$(date '+%F %T') 122 done" || echo "$(date '+%F %T') 122 FAILED (see $L/122_null_scan.log)"
touch $L/DONE; echo "$(date '+%F %T') end"
