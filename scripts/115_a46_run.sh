#!/bin/bash
# Amendment 46: late analysis, run after all primary results were seen. Writes only to work/a46_inputs/ (dumps, never committed) and results/late_a46/; never overwrites.
# Usage (inside systemd-run, see the launch block): 115_a46_run.sh
#   1. dumps: 21/65/81 with A46=1 for the 12 primary target tables, difford_arch (M1) and Wallace depth-adjusted (scan only)
#   2. 46.3 ANCOM-BC2: two Freedman-Lane fits per methane table (seeds 46101, 46102), timed; the other 198 per table only if the projected time on P workers is under 24 h
#   3. while those run: 112 (46.1, 46.3 CLR), 113 (46.2), 114 (46.4), one after another on C cores; then the ANCOM-BC2 summary (112 ancombc2)
# Env: P (ANCOM-BC2 workers, default 10), C (cores for 113 and 114, default 4).
set -u; cd /data/rumen_multicohort; P=${P:-10}; C=${C:-4}
MB=/home/rlagad/.local/share/mamba/envs; RR=$MB/rumen-r/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RR
L=logs/late_a46; mkdir -p $L/dump $L/ancombc2_fl; export L
grep -q '^\*\*Amendment 46 (' plan/ANALYSIS_PLAN.md && git diff --quiet HEAD -- plan/ANALYSIS_PLAN.md && git ls-files --error-unmatch scripts/lib_a46.R scripts/115_a46_run.sh >/dev/null 2>&1 \
  && git diff --quiet HEAD -- scripts/lib_a46.R scripts/lib_clr.R scripts/21_estimate_targets.R scripts/65_estimate_alberta.R scripts/81_estimate_methane.R scripts/112_a46_effect_size.R scripts/113_a46_coverage.R scripts/114_a46_scan.R scripts/115_a46_run.sh \
  || { echo "Amendment 46 and its scripts must be committed (and unmodified) before any run"; exit 1; }
echo "$(date '+%F %T') start | Amendment 46 commit $(git log -1 --format=%h -- plan/ANALYSIS_PLAN.md) | scripts commit $(git log -1 --format=%h -- scripts/lib_a46.R) | P $P | C $C"
dump() { set -- $1; t=$1; a=""
  case $t in
    nelore_*|charolais_mrna_*) cmd="scripts/21_estimate_targets.R ${t%_*} ${t##*_}";;
    alberta_*_all) r=${t#alberta_}; cmd="scripts/65_estimate_alberta.R ${r%_all}"; a="ALBERTA_ALL=1";;
    wallace_prokaryotes_depth) cmd="scripts/81_estimate_methane.R wallace_prokaryotes real"; a="DEPTH=1";;
    difford_bact|difford_arch|wallace_prokaryotes) cmd="scripts/81_estimate_methane.R $t real";;
    *) echo "unknown table $t"; return 1;;
  esac
  env $a A46=1 $RR $cmd > $L/dump/$t.log 2>&1 || echo "$(date +%T) dump $t FAILED (see $L/dump/$t.log)"; }
flfit() { set -- $1; t=$1; s=$2; a=$(date +%s)
  A46_FL_SEED=$s $RR scripts/81_estimate_methane.R $t real > $L/ancombc2_fl/${t}_$s.log 2>&1 || echo "$(date +%T) ANCOM-BC2 $t $s FAILED"
  echo -e "$t\t$s\t$(( $(date +%s) - a ))" >> results/late_a46/ancombc2_fl/timing.tsv; }
export -f dump flfit
TABLES="nelore_species nelore_genus nelore_family nelore_phylum charolais_mrna_species charolais_mrna_genus charolais_mrna_family charolais_mrna_phylum alberta_genus_all alberta_family_all difford_bact difford_arch wallace_prokaryotes wallace_prokaryotes_depth"
for t in $TABLES; do echo $t; done | xargs -P 4 -I{} bash -c 'dump "{}"'
n=$(ls work/a46_inputs/*.rds 2>/dev/null | wc -l); echo "$(date '+%F %T') dumps: $n of 14"; [ "$n" -eq 14 ] || { echo "dumps incomplete, stopping"; exit 1; }
mkdir -p results/late_a46/ancombc2_fl
for t in difford_bact wallace_prokaryotes; do for s in 46101 46102; do echo "$t $s"; done; done | xargs -P 4 -I{} bash -c 'flfit "{}"'
K=$($RR -e 'x <- read.delim("results/late_a46/ancombc2_fl/timing.tsv", header=FALSE, col.names=c("table", "seed", "sec")); per <- tapply(x$sec, x$table, max)
  P <- as.numeric(Sys.getenv("P", "10")); h <- sum(per * 198) / P / 3600; run <- length(per) == 2 && h < 24
  writeLines(c("late analysis, run after all primary results were seen (Amendment 46.3)", sprintf("seconds per ANCOM-BC2 fit (slowest of two): %s", paste(names(per), per, collapse="; ")),
    sprintf("projected wall time for the remaining 2 x 198 Freedman-Lane fits on %g workers: %.1f h", P, h), sprintf("rule: 200 permutations only if under 24 h -> %s", if (run) "run" else "not run")),
    "results/late_a46/ancombc2_fl/permutation_decision.txt"); cat(if (run) 1 else 0)' 2>>$L/ancombc2_fl/decision.err)
cat results/late_a46/ancombc2_fl/permutation_decision.txt
if [ "$K" = 1 ]; then { for s in $(seq 46103 46300); do for t in difford_bact wallace_prokaryotes; do echo "$t $s"; done; done | xargs -P $P -I{} bash -c 'flfit "{}"'; echo "$(date '+%F %T') ANCOM-BC2 Freedman-Lane fits done"; } & BG=$!; fi
$RR scripts/112_a46_effect_size.R > $L/112_effect_size.log 2>&1 && echo "$(date '+%F %T') 112 done" || echo "112 FAILED (see $L/112_effect_size.log)"
A46_CORES=$C $RR scripts/113_a46_coverage.R > $L/113_coverage.log 2>&1 && echo "$(date '+%F %T') 113 done" || echo "113 FAILED (see $L/113_coverage.log)"
A46_CORES=$C $RR scripts/114_a46_scan.R > $L/114_scan.log 2>&1 && echo "$(date '+%F %T') 114 done" || echo "114 FAILED (see $L/114_scan.log)"
[ -n "${BG:-}" ] && wait $BG
echo "$(date '+%F %T') ANCOM-BC2 Freedman-Lane files: $(ls results/late_a46/ancombc2_fl/*_fl*_targets.tsv 2>/dev/null | wc -l)"
$RR scripts/112_a46_effect_size.R ancombc2 > $L/112_ancombc2.log 2>&1 && echo "$(date '+%F %T') ANCOM-BC2 summary done" || echo "112 ancombc2 FAILED (see $L/112_ancombc2.log)"
touch $L/DONE; echo "$(date '+%F %T') end"
