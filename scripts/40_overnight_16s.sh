#!/bin/bash
# Overnight driver, Charolais 16S layer (Amendment 17). Stops at the first failure; commits nothing.
set -euo pipefail; cd /data/rumen_multicohort; L=logs/charolais16s; W=work/charolais16s; mkdir -p $L $W results/charolais16s
grep -q "^\*\*Amendment 17 (" plan/ANALYSIS_PLAN.md && git diff --quiet plan/ANALYSIS_PLAN.md || { echo "commit Amendment 17 first"; exit 1; }
grep -q charolais_16s scripts/21_estimate_targets.R || { echo "estimator not patched"; exit 1; }
t() { echo "$(date +%F_%T) $*"; }
[ -f $W/.trim_ok ] || { t trim; bash scripts/30_trim_charolais16s.sh > $L/trim.log 2>&1; touch $W/.trim_ok; }
[ -f $W/.silva_ok ] || { t silva; bash scripts/31_silva_trainset.sh > $L/silva.log 2>&1; touch $W/.silva_ok; }
t dada2; mamba run -p /home/rlagad/mambaforge/envs/silva_taxonomy Rscript scripts/32_dada2_charolais16s.R > $L/dada2.log 2>&1; tail -n 3 $L/dada2.log
t estimate; p=""
for rk in genus family; do PROFILE=profiles_16s mamba run -n rumen-r Rscript scripts/21_estimate_targets.R charolais_16s $rk > $L/est_$rk.log 2>&1 & p="$p $!"; done
for i in $p; do wait $i; done
t report; mamba run -n rumen-r Rscript scripts/33_report_charolais16s.R 2>&1 | tee results/charolais16s/summary.txt
touch $L/DONE; t DONE
