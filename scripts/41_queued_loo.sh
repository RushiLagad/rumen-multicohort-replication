#!/bin/bash
# Queued behind the 16S job. Amendment 18 (exploratory): leave-one-animal-out for Charolais mRNA.
set -euo pipefail; cd /data/rumen_multicohort; L=logs/loo; mkdir -p $L results/loo
grep -q "^\*\*Amendment 18 (" plan/ANALYSIS_PLAN.md && git diff --quiet plan/ANALYSIS_PLAN.md || { echo "commit Amendment 18 first"; exit 1; }
t() { echo "$(date +%F_%T) $*"; }; t waiting for the 16S job
while pgrep -f 40_overnight_16s.sh > /dev/null; do sleep 60; done
awk -F'\t' '$1=="Charolais_PRJEB70895" && $3=="mRNA"{print $2}' harmonized/runs.tsv | sort -u > $L/animals.txt; [ "$(wc -l < $L/animals.txt)" = 30 ]
t fits; for rk in genus family phylum; do while read a; do echo "$a $rk"; done < $L/animals.txt; done | \
  xargs -P 8 -L1 bash -c 'DROP=$0 mamba run -n rumen-r Rscript scripts/35_estimate_loo.R charolais_mrna $1 > logs/loo/drop_$0_$1.log 2>&1' || echo "some fits failed; see logs/loo"
t report; mamba run -n rumen-r Rscript scripts/34_report_loo.R 2>&1 | tee results/loo/summary.txt
touch $L/DONE; t DONE
