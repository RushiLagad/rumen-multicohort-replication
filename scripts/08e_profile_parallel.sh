#!/usr/bin/env bash
# Parallel version of 08d: xargs hands each sample to exactly one worker. Resumable (bracken_G = done). Usage: P=6 bash 08e_profile_parallel.sh
set -uo pipefail
DB=${DB:-/data/refs/k2_gtdb207_full}; ROOT=${ROOT:-/data/rumen_multicohort}; OUT=${OUT:-$ROOT/profiles_gtdb_full}; P=${P:-5}
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
export DB ROOT OUT
for c in nelore charolais_mrna li; do
  mkdir -p $OUT/$c
  f1=$(ls $ROOT/processed/$c/*.f*q.gz | head -1)
  L=$(zcat -f $f1 | head -8000 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
  [ -s $DB/database${L}mers.kmer_distrib ] || { say "FATAL: no Bracken distribution for $c L=$L"; exit 1; }
  mapfile -t R1 < <(ls $ROOT/processed/$c/* | grep -E '_1\.(fastq|fq)(\.gz)?$')
  [ ${#R1[@]} -eq 0 ] && mapfile -t R1 < <(ls $ROOT/processed/$c/*.f*q.gz)
  say "$c: read length $L, ${#R1[@]} samples, $P at a time"
  printf '%s\n' "${R1[@]}" | xargs -P $P -I{} bash $ROOT/scripts/08e_one_sample.sh $c {} $L
  n=$(ls $OUT/$c/*.bracken_G 2>/dev/null | wc -l); [ $n -eq ${#R1[@]} ] || { say "FATAL: $c has $n/${#R1[@]} samples"; exit 1; }
  python3 $ROOT/scripts/08c_matrices.py $OUT/$c $c
  if [ $c = nelore ]; then
    say "--- adequacy test (Amendment 1): median classified fraction, Nelore, threshold 60% ---"
python3 - <<PY
import re, glob, statistics
v = [float(re.search(r'sequences classified \(([\d.]+)%\)', open(f).read()).group(1)) for f in sorted(glob.glob("$OUT/nelore/*.kraken.log"))]
print(f"n={len(v)} median={statistics.median(v):.1f}% range {min(v):.1f}-{max(v):.1f}% ->", "PASS" if statistics.median(v) >= 60 else "FAIL")
PY
  fi
done
say "PROFILING COMPLETE"
