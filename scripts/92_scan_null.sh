#!/bin/bash
# Amendment 38 item 3: 20 within-design permutations with the full ANCOM-BC2 table kept, RFI genus tables (methane tables already have 200 in results/methane_dryrun/). Resumable.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; mkdir -p results/scan_null logs/scan_null
job() { set -- $1; k=$1; s=$2
  case $k in
    nelore_genus|charolais_mrna_genus) c=${k%_*}; o=results/scan_null/${k}_seed${s}_ancombc2.csv
      [ -s $o ] || SHUFFLE_SEED=$s KEEP_CSV=1 $RR scripts/21_estimate_targets.R $c genus > logs/scan_null/${k}_seed${s}.log 2>&1 || echo "$k $s FAILED";;
    alberta_genus) o=results/scan_null/alberta_genus_seed${s}_ancombc2.csv
      [ -s $o ] || SHUFFLE_SEED=$s KEEP_CSV=1 $RR scripts/65_estimate_alberta.R genus > logs/scan_null/${k}_seed${s}.log 2>&1 || echo "$k $s FAILED";;
  esac; }
export -f job; export RR
for s in $(seq 2001 2200); do for k in nelore_genus charolais_mrna_genus alberta_genus; do echo "$k $s"; done; done | xargs -P 12 -I{} bash -c 'job "{}"'
touch logs/scan_null/DONE; echo "$(date +%F_%T) scan nulls done: $(ls results/scan_null/*_ancombc2.csv | wc -l) tables"
