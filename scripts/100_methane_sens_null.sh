#!/bin/bash
# Calibration nulls for the two pre-registered Wallace methane sensitivities (no SE1; CH4 per kg DMI), 200 within-farm permutations each. Resumable.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; mkdir -p logs/methane_sens
job() { set -- $1; v=$1; s=$2; o=results/methane_dryrun/wallace_prokaryotes_${v}_seed${s}_targets.tsv; [ -s $o ] && return
  if [ $v = noSE1 ]; then NO_SE1=1 SEED=$s $RR scripts/81_estimate_methane.R wallace_prokaryotes shuffle > logs/methane_sens/${v}_${s}.log 2>&1 || echo "$v $s FAILED"
  else TRAIT=ch4_dmi SEED=$s $RR scripts/81_estimate_methane.R wallace_prokaryotes shuffle > logs/methane_sens/${v}_${s}.log 2>&1 || echo "$v $s FAILED"; fi; }
export -f job; export RR
for s in $(seq 1001 1200); do echo "noSE1 $s"; echo "ch4_dmi $s"; done | xargs -P 12 -I{} bash -c 'job "{}"'
touch logs/methane_sens/DONE; echo "$(date +%F_%T) done: noSE1 $(ls results/methane_dryrun/wallace_prokaryotes_noSE1_seed*_targets.tsv | wc -l), ch4_dmi $(ls results/methane_dryrun/wallace_prokaryotes_ch4_dmi_seed*_targets.tsv | wc -l)"
