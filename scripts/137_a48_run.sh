#!/bin/bash
# Amendment 48: late analysis, run after all primary, late, Amendment 46 and Amendment 47 results were seen. Writes only work/a46_inputs/wallace_prokaryotes_ch4_dmi.rds
# (never committed) and results/late_a48/; never overwrites. Usage (inside systemd-run, see the launch block): 137_a48_run.sh
#   0. dump of the Wallace CH4 per kg DMI table (81 with A48=1 TRAIT=ch4_dmi) for 48.1 and 48.9
#   1. priority 1: 128 (48.1 seven claims), 129 (48.2 ANCOM-BC2 filter), 130 (48.3 harmonized dairy)
#   2. priority 2: 131 (48.4 farms), 132 (48.5 CLR grid), 133 (48.6 Charolais influence, after 132), 134 (48.7 power)
#   3. 135 (48.8 mapping), 136 (48.9 clade ratio)
# A failed step is logged and the next one runs (none depends on another except 133 on 132); the run ends with DONE.
set -u; cd /data/rumen_multicohort
MB=/home/rlagad/.local/share/mamba/envs; RR=$MB/rumen-r/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
L=logs/late_a48; mkdir -p $L
S="scripts/lib_a46.R scripts/lib_a47.R scripts/lib_a48.R scripts/81_estimate_methane.R scripts/128_a48_seven.R scripts/129_a48_ancombc2.R scripts/130_a48_dairy.R scripts/131_a48_farms.R scripts/132_a48_clr_grid.R scripts/133_a48_charolais.R scripts/134_a48_power.R scripts/135_a48_mapping.R scripts/136_a48_clade_ratio.R scripts/137_a48_run.sh"
grep -q '^\*\*Amendment 48 (' plan/ANALYSIS_PLAN.md && git diff --quiet HEAD -- plan/ANALYSIS_PLAN.md && git ls-files --error-unmatch $S >/dev/null 2>&1 && git diff --quiet HEAD -- $S \
  || { echo "Amendment 48 and its scripts must be committed (and unmodified) before any run"; exit 1; }
echo "$(date '+%F %T') start | Amendment 48 commit $(git log -1 --format=%h -- plan/ANALYSIS_PLAN.md) | scripts commit $(git log -1 --format=%h -- scripts/lib_a48.R)"
A48=1 TRAIT=ch4_dmi $RR scripts/81_estimate_methane.R wallace_prokaryotes real > $L/dump_wallace_ch4_dmi.log 2>&1
[ -s work/a46_inputs/wallace_prokaryotes_ch4_dmi.rds ] || { echo "$(date '+%F %T') Wallace CH4 per kg DMI dump FAILED (see $L/dump_wallace_ch4_dmi.log)"; exit 1; }
echo "$(date '+%F %T') dump done: $(grep -o 'animals [0-9]* | features [0-9]* | df [0-9]*' $L/dump_wallace_ch4_dmi.log)"
run() { local s=$1 n=$2; $RR scripts/$s > $L/$n.log 2>&1 && echo "$(date '+%F %T') $n done" || echo "$(date '+%F %T') $n FAILED (see $L/$n.log)"; }
run 128_a48_seven.R 128_seven
run 129_a48_ancombc2.R 129_ancombc2
run 130_a48_dairy.R 130_dairy
run 131_a48_farms.R 131_farms
run 132_a48_clr_grid.R 132_clr_grid
run 133_a48_charolais.R 133_charolais
run 134_a48_power.R 134_power
run 135_a48_mapping.R 135_mapping
run 136_a48_clade_ratio.R 136_clade_ratio
touch $L/DONE; echo "$(date '+%F %T') end"
