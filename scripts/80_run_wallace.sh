#!/bin/bash
# Amendments 28 and 33: Wallace markers one after another; stops at the first failure.
cd /data/rumen_multicohort; RS=/home/rlagad/mambaforge/envs/silva_taxonomy/bin/Rscript; mkdir -p logs/16s/wallace
for m in Prokaryotes Archaea Bacteria; do
  echo "$(date +%F_%T) $m start"
  $RS scripts/79_dada2_wallace.R $m > logs/16s/wallace/${m}_dada2.log 2>&1 || { echo "$(date +%F_%T) $m FAILED"; exit 1; }
  touch logs/16s/wallace/${m}.DONE; echo "$(date +%F_%T) $m done"
done
