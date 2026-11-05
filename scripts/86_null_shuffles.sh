#!/bin/bash
# Amendment 37 (pending): empirical null for ANCOM-BC2 W per target, K within-farm/herd permutations per cohort.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; K=${K:-30}
for c in wallace_prokaryotes difford_bact difford_arch; do for s in $(seq 1001 $((1000 + K))); do
  [ -s results/methane_dryrun/${c}_seed${s}_targets.tsv ] || SEED=$s $RR scripts/81_estimate_methane.R $c shuffle > logs/methane_dryrun/${c}_seed${s}.log 2>&1 || echo "$c seed $s FAILED"
done; echo "$(date +%F_%T) $c: $K shuffles done"; done; touch logs/methane_dryrun/NULL_DONE
