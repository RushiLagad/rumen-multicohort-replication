#!/bin/bash
# Amendment 37 item 4: 200 within-design permutations for every RFI table; 12 fits in parallel. Resumable.
cd /data/rumen_multicohort; RR=/home/rlagad/.local/share/mamba/envs/rumen-r/bin/Rscript; mkdir -p results/rfi_null logs/rfi_null
job() { set -- $1; k=$1; s=$2
  case $k in
    nelore_*|charolais_mrna_*) c=${k%_*}; r=${k##*_}; o=results/rfi_null/${k}_seed${s}_targets.tsv
      [ -s $o ] || SHUFFLE_SEED=$s $RR scripts/21_estimate_targets.R $c $r > logs/rfi_null/${k}_seed${s}.log 2>&1 || echo "$k $s FAILED";;
    alberta_*_all) r=$(echo $k | cut -d_ -f2); o=results/rfi_null/alberta_${r}_all_seed${s}_targets.tsv
      [ -s $o ] || SHUFFLE_SEED=$s ALBERTA_ALL=1 $RR scripts/65_estimate_alberta.R $r > logs/rfi_null/${k}_seed${s}.log 2>&1 || echo "$k $s FAILED";;
    alberta_*) r=${k#alberta_}; o=results/rfi_null/alberta_${r}_seed${s}_targets.tsv
      [ -s $o ] || SHUFFLE_SEED=$s $RR scripts/65_estimate_alberta.R $r > logs/rfi_null/${k}_seed${s}.log 2>&1 || echo "$k $s FAILED";;
  esac; }
export -f job; export RR
# slow tables first so they do not trail at the end
{ for s in $(seq 1001 1200); do for k in nelore_species nelore_genus; do echo "$k $s"; done; done
  for s in $(seq 1001 1200); do for k in nelore_family nelore_phylum charolais_mrna_species charolais_mrna_genus charolais_mrna_family charolais_mrna_phylum alberta_genus alberta_family alberta_genus_all alberta_family_all; do echo "$k $s"; done; done
} | xargs -P 12 -I{} bash -c 'job "{}"'
touch logs/rfi_null/DONE; echo "$(date +%F_%T) RFI nulls done"
