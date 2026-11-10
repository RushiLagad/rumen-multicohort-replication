#!/bin/bash
# Amendment 43: 200 within-farm permutations for the Wallace RFI NRC sensitivity fit. Resumable.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; mkdir -p logs/wallace_rfi
job() { s=$1; o=results/methane_dryrun/wallace_prokaryotes_rfi_nrc_seed${s}_targets.tsv
  [ -s $o ] || TRAIT=rfi_nrc SEED=$s $RR scripts/81_estimate_methane.R wallace_prokaryotes shuffle > logs/wallace_rfi/seed${s}.log 2>&1 || echo "seed $s FAILED"; }
export -f job; export RR
seq 1001 1200 | xargs -P 12 -I{} bash -c 'job {}'
touch logs/wallace_rfi/DONE; echo "$(date +%F_%T) Wallace RFI nulls done: $(ls results/methane_dryrun/wallace_prokaryotes_rfi_nrc_seed*_targets.tsv | wc -l)"
