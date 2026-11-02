#!/bin/bash
# Run all target estimates (8 jobs in parallel). Dry run first: scripts/21_run_all.sh dry
cd /data/rumen_multicohort && mkdir -p logs/replication
M=${1:-fit}
run1() { c=$1; r=$2; mamba run -n rumen-r Rscript scripts/21_estimate_targets.R $c $r $M > logs/replication/${c}_${r}_$M.log 2>&1
  echo "$(date +%T) $c $r $M exit=$?"; }
export -f run1; export M
if [ "$M" = dry ]; then run1 nelore family; for p in "nelore genus" "nelore phylum" "nelore species" "charolais_mrna genus" "charolais_mrna family" "charolais_mrna phylum" "charolais_mrna species"; do run1 $p; done
else for c in nelore charolais_mrna; do for r in genus family phylum species; do echo "$c $r"; done; done | xargs -P 8 -L1 bash -c 'run1 $0 $1'; fi
