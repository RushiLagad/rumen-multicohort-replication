#!/bin/bash
# Amendment 41: Wallace depth-adjusted real fit plus 200 within-farm permutations. Resumable.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; mkdir -p logs/depth_sens
[ -s results/methane/wallace_prokaryotes_depth_targets.tsv ] || DEPTH=1 $RR scripts/81_estimate_methane.R wallace_prokaryotes real > logs/depth_sens/real.log 2>&1 || echo "REAL FAILED"
job() { s=$1; o=results/methane_dryrun/wallace_prokaryotes_depth_seed${s}_targets.tsv
  [ -s $o ] || DEPTH=1 SEED=$s $RR scripts/81_estimate_methane.R wallace_prokaryotes shuffle > logs/depth_sens/seed${s}.log 2>&1 || echo "seed $s FAILED"; }
export -f job; export RR
seq 1001 1200 | xargs -P 12 -I{} bash -c 'job {}'
touch logs/depth_sens/DONE; echo "$(date +%F_%T) depth nulls done: $(ls results/methane_dryrun/wallace_prokaryotes_depth_seed*_targets.tsv | wc -l)"
