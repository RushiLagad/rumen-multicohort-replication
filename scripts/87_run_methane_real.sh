#!/bin/bash
# Amendment 37 item 5: real methane fits only after the 200-permutation nulls exist; then pooling.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript
[ -f logs/methane_dryrun/NULL200_DONE ] || { echo "nulls not finished"; exit 1; }
mkdir -p logs/methane
for c in difford_arch difford_bact wallace_prokaryotes wallace_bacteria; do echo "$(date +%F_%T) $c real"; $RR scripts/81_estimate_methane.R $c real > logs/methane/$c.log 2>&1 || { echo "$c FAILED"; exit 1; }; done
echo "$(date +%F_%T) wallace sensitivity: no SE1, CH4 per DMI"
NO_SE1=1 $RR scripts/81_estimate_methane.R wallace_prokaryotes real > logs/methane/wallace_prokaryotes_noSE1.log 2>&1 || echo "noSE1 FAILED"
TRAIT=ch4_dmi $RR scripts/81_estimate_methane.R wallace_prokaryotes real > logs/methane/wallace_prokaryotes_ch4_dmi.log 2>&1 || echo "ch4_dmi FAILED"
$RR scripts/83_pool_methane.R real > logs/methane/pool.log 2>&1 && echo "$(date +%F_%T) POOLED" || echo "pooling FAILED"
