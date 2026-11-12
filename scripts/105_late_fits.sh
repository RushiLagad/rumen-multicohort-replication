#!/bin/bash
# Amendment 44: late analyses, run after all primary results were seen. Writes only to results/late_*; resumable (existing outputs are skipped, never overwritten).
# Usage (inside systemd-run, see Amendment 44): 105_late_fits.sh <linda|aldex2|hedges|diversity>
#   linda:     real fit + 200 within-design permutations for the 12 primary tables, then pooling (106).
#   aldex2:    real fit + 2 permutations per table (timing), then 200 permutations if the projected run is under 12 h, else 50; then pooling (106).
#   hedges:    Charolais mRNA Hedges' g, four ranks.   diversity: Li et al. diversity indices in the three beef cohorts, then pooling (107).
#   Amendment 45: poscontrol (108, no new fits); mapping (Charolais mRNA and 16S target CLR, then 109); spikein (Wallace, Alberta 529, Nelore genus; 20 permuted bases x 3 effect sizes; then 110).
set -u; cd /data/rumen_multicohort; M=${1:?method}; P=${P:-10}
MB=/home/rlagad/.local/share/mamba/envs
RRR=$MB/rumen-r/bin/Rscript; if [ "$M" = linda ]; then RR=$MB/rumen-late/bin/Rscript; else RR=$RRR; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RR M P
L=logs/late_$M; mkdir -p $L; export L
grep -q '^\*\*Amendment 44 (' plan/ANALYSIS_PLAN.md && { case $M in poscontrol|mapping|spikein) grep -q '^\*\*Amendment 45 (' plan/ANALYSIS_PLAN.md;; esac; } && git diff --quiet HEAD -- plan/ANALYSIS_PLAN.md && git ls-files --error-unmatch scripts/lib_late.R >/dev/null 2>&1 \
  && git diff --quiet HEAD -- scripts/lib_late.R scripts/lib_replication.R scripts/21_estimate_targets.R scripts/65_estimate_alberta.R scripts/81_estimate_methane.R scripts/105_late_fits.sh scripts/106_late_pool.R scripts/107_late_hedges_diversity.R scripts/108_late_positive_controls.R scripts/109_late_mapping.R scripts/110_late_spikein.R || { echo "Amendment 44 and the scripts must be committed (and unmodified) before any late run"; exit 1; }
echo "$(date '+%F %T') start $M | Amendment 44 commit $(git log -1 --format=%h -- plan/ANALYSIS_PLAN.md) | Rscript $RR | P $P"
unit() { set -- $1; t=$1; s=$2; a=""
  case $t in
    nelore_*|charolais_mrna_*) cmd="scripts/21_estimate_targets.R ${t%_*} ${t##*_}";;
    charolais_16s_*) cmd="scripts/21_estimate_targets.R charolais_16s ${t##*_}"; a="PROFILE=profiles_16s";;
    alberta529_*) cmd="scripts/65_estimate_alberta.R ${t#alberta529_}";;   # Alberta without the 35 ID matches: the table whose 200-permutation null is on file (scan_null)
    alberta_*) cmd="scripts/65_estimate_alberta.R ${t#alberta_}"; a="ALBERTA_ALL=1";;
    difford_bact|wallace_prokaryotes) cmd="scripts/81_estimate_methane.R $t real";;
    *) echo "unknown table $t"; return 1;;
  esac
  env $a LATE=$M LATE_SEEDS=$s $RR $cmd > $L/${t}_${s/:/-}.log 2>&1 || echo "$(date +%T) $t $s FAILED (see $L/${t}_${s/:/-}.log)"; }
export -f unit
TABLES="nelore_species nelore_genus charolais_mrna_species charolais_mrna_genus wallace_prokaryotes difford_bact nelore_family nelore_phylum charolais_mrna_family charolais_mrna_phylum alberta_genus alberta_family"
chunks() { a=$1; b=$2; shift 2; tl=${*:-$TABLES}; while [ $a -le $b ]; do e=$((a + 24)); [ $e -gt $b ] && e=$b; for t in $tl; do echo "$t $a:$e"; done; a=$((e + 1)); done; }
# linda() calls vcov() on one multivariate lm, a (features x coefficients)^2 matrix: about 15 to 30 GB for the species tables, so they run on P_SPECIES workers after the rest.
SPECIES="nelore_species charolais_mrna_species"; NONSP=$(for t in $TABLES; do case $t in *_species) ;; *) echo -n "$t ";; esac; done)
case $M in
  linda)
    { for t in $NONSP; do echo "$t real"; done; chunks 1001 1200 $NONSP; } | xargs -P $P -I{} bash -c 'unit "{}"'
    echo "$(date '+%F %T') non-species tables done; species tables on ${P_SPECIES:-2} workers"
    { for t in $SPECIES; do echo "$t real"; done; chunks 1001 1200 $SPECIES; } | xargs -P ${P_SPECIES:-2} -I{} bash -c 'unit "{}"'
    K=200;;
  aldex2)
    { for t in $TABLES; do echo "$t real"; echo "$t 1001:1002"; done; } | xargs -P $P -I{} bash -c 'unit "{}"'
    K=$($RR -e 'x <- read.delim("results/late_aldex2/timing.tsv", header=FALSE, col.names=c("table","unit","sec","features","animals"))
      p <- x[grepl("^perm_x", x$unit), ]; p$per <- p$sec / as.integer(sub("perm_x", "", p$unit)); per <- tapply(p$per, p$table, max)
      if (length(per) < 12) per <- c(per, setNames(rep(max(per), 12 - length(per)), paste0("untimed", seq_len(12 - length(per))))); h <- sum(per * 198) / as.numeric(Sys.getenv("P", "10")) / 3600; hl <- max(per * 25) / 3600
      K <- if (h + hl < 12) 200 else 50
      writeLines(c(sprintf("projected wall time for 198 further permutations per table on %s workers: %.2f h (+ %.2f h for the slowest 25-seed chunk)", Sys.getenv("P", "10"), h, hl),
        sprintf("rule (Amendment 44 item 4): 200 if under 12 h, else 50 -> K = %d", K), sprintf("seconds per permutation: %s", paste(names(per), round(per, 1), collapse="; "))), "results/late_aldex2/permutation_count.txt")
      cat(K)')
    cat results/late_aldex2/permutation_count.txt
    chunks 1003 $((1000 + K)) | xargs -P $P -I{} bash -c 'unit "{}"';;
  hedges)
    for r in species genus family phylum; do unit "charolais_mrna_$r real"; done;;
  diversity)
    for t in nelore_genus charolais_mrna_genus alberta_genus; do unit "$t real"; done
    $RRR scripts/107_late_hedges_diversity.R > $L/pool.log 2>&1 && echo "$(date '+%F %T') diversity pooled" || echo "107 FAILED";;
  poscontrol)
    $RRR scripts/108_late_positive_controls.R > $L/positive_controls.log 2>&1 && echo "$(date '+%F %T') positive controls done" || echo "108 FAILED";;
  mapping)
    for t in charolais_mrna_genus charolais_mrna_family charolais_mrna_phylum charolais_mrna_species charolais_16s_genus charolais_16s_family charolais_16s_phylum; do unit "$t real"; done
    $RRR scripts/109_late_mapping.R > $L/mapping.log 2>&1 && echo "$(date '+%F %T') mapping summary done" || echo "109 FAILED";;
  spikein)
    for s in 4001 4006 4011 4016; do for t in nelore_genus wallace_prokaryotes alberta529_genus; do echo "$t $s:$((s + 4))"; done; done | xargs -P $P -I{} bash -c 'unit "{}"'
    echo "$(date '+%F %T') spike-in fits: $(find results/late_spikein -name 'b*_seed*.tsv' | wc -l) of 180"
    $RRR scripts/110_late_spikein.R > $L/spikein_summary.log 2>&1 && echo "$(date '+%F %T') spike-in summary done" || echo "110 FAILED";;
  *) echo "unknown method $M"; exit 1;;
esac
if [ "$M" = linda ] || [ "$M" = aldex2 ]; then
  echo "$(date '+%F %T') fits done: real $(ls results/late_$M/*_targets.tsv 2>/dev/null | wc -l) of 12; null files $(find results/late_$M/null -name 'seed*.tsv' 2>/dev/null | wc -l) (expected $((12 * K)))"
  $RRR scripts/106_late_pool.R $M > $L/pool.log 2>&1 && echo "$(date '+%F %T') $M pooled" || echo "106 FAILED (see $L/pool.log)"
fi
touch $L/DONE; echo "$(date '+%F %T') end $M"
