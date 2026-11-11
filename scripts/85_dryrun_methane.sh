#!/bin/bash
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript
for c in difford_arch difford_bact wallace_prokaryotes; do echo "$(date +%F_%T) $c start"; $RR scripts/81_estimate_methane.R $c shuffle > logs/methane_dryrun/$c.log 2>&1 && echo "$(date +%F_%T) $c done" || { echo "$(date +%F_%T) $c FAILED"; exit 1; }; done
touch logs/methane_dryrun/DONE
