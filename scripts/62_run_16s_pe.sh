#!/bin/bash
# Amendment 21 driver for one paired 16S cohort. Usage: 62_run_16s_pe.sh <BioProject> <cohort> <fwd primer> <rev primer> <expected insert bp>
set -euo pipefail; cd ${RUMEN_ROOT:-/data/rumen_multicohort}; P=$1; C=$2; L=logs/16s/$P; mkdir -p $L
grep -q "^\*\*Amendment 21 (" plan/ANALYSIS_PLAN.md && git diff --quiet plan/ANALYSIS_PLAN.md || { echo "commit Amendment 21 first"; exit 1; }
[ -s work/charolais16s/silva138_train.fa ] || { echo "SILVA training file missing"; exit 1; }
t() { echo "$(date +%F_%T) $*"; }
t trim; bash scripts/60_trim_16s_pe.sh $P "$3" "$4" > $L/trim.log 2>&1; tail -n 3 $L/trim.log
t dada2; mamba run -p /home/rlagad/mambaforge/envs/silva_taxonomy Rscript scripts/61_dada2_16s_pe.R $P $C $5 > $L/dada2.log 2>&1; grep -E "truncLen|merge rate|FALLBACK|^mode" $L/dada2.log
touch $L/DONE; t DONE
