#!/bin/bash
# Amendment 37: K within-farm/herd permutations per table; three tables run in parallel. Usage: K=200 86_null_shuffles.sh
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; K=${K:-200}
one() { c=$1; for s in $(seq 1001 $((1000 + K))); do
  [ -s results/methane_dryrun/${c}_seed${s}_targets.tsv ] || SEED=$s $RR scripts/81_estimate_methane.R $c shuffle > logs/methane_dryrun/${c}_seed${s}.log 2>&1 || echo "$c seed $s FAILED"
done; echo "$(date +%F_%T) $c: $K shuffles done"; }
export -f one; export K RR
printf "%s\n" wallace_prokaryotes wallace_bacteria difford_bact difford_arch | xargs -P 3 -I{} bash -c 'one {}'
touch logs/methane_dryrun/NULL200_DONE; echo "$(date +%F_%T) all done"
