#!/usr/bin/env bash
# Kraken2 + Bracken on all cohorts against k2_gtdb207. Resumable. Builds Bracken distributions per read length.
set -uo pipefail
DB=/data/refs/k2_gtdb207; T=24; ROOT=/data/rumen_multicohort; OUT=$ROOT/profiles_gtdb
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
mkdir -p $OUT
for c in nelore charolais_mrna li; do
  mkdir -p $OUT/$c
  f1=$(ls $ROOT/processed/$c/*.f*q.gz | head -1)
  L=$(zcat -f $f1 | head -8000 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
  if [ ! -s $DB/database${L}mers.kmer_distrib ]; then
    say "$c: building Bracken distribution for read length $L"
    kmer2read_distr --seqid2taxid $DB/seqid2taxid.map --taxonomy $DB/taxonomy --kraken $DB/database.kraken \
      --output $DB/database${L}mers.kraken -k 35 -l $L -t $T > $DB/kmer2read_${L}.log 2>&1 || { say "FATAL kmer2read_distr L=$L"; exit 1; }
    generate_kmer_distribution.py -i $DB/database${L}mers.kraken -o $DB/database${L}mers.kmer_distrib \
      > $DB/kmer_distrib_${L}.log 2>&1 || { say "FATAL generate_kmer_distribution L=$L"; exit 1; }
  fi
  say "$c: read length $L, distribution ready"
  mapfile -t R1 < <(ls $ROOT/processed/$c/* | grep -E '_1\.(fastq|fq)(\.gz)?$')
  [ ${#R1[@]} -eq 0 ] && mapfile -t R1 < <(ls $ROOT/processed/$c/*.f*q.gz)
  for r1 in "${R1[@]}"; do
    r2=${r1/_1./_2.}; s=$(basename $r1 | sed -E 's/(_1)?\.(fastq|fq)(\.gz)?$//')
    o=$OUT/$c/$s
    [ -s $o.bracken_G ] && continue
    gz=""; [[ $r1 == *.gz ]] && gz="--gzip-compressed"
    if [ "$r2" != "$r1" ] && [ -s "$r2" ]; then
      kraken2 --db $DB --threads $T --memory-mapping $gz --paired --report $o.kreport --output /dev/null $r1 $r2 2> $o.kraken.log
    else
      kraken2 --db $DB --threads $T --memory-mapping $gz --report $o.kreport --output /dev/null $r1 2> $o.kraken.log
    fi
    [ -s $o.kreport ] || { say "FATAL kraken2 failed on $c/$s"; tail -2 $o.kraken.log; exit 1; }
    bracken -d $DB -i $o.kreport -o $o.bracken_S -w $o.bracken_S.report -r $L -l S -t 10 > $o.bracken.log 2>&1
    bracken -d $DB -i $o.kreport -o $o.bracken_G -w $o.bracken_G.report -r $L -l G -t 10 >> $o.bracken.log 2>&1
    [ -s $o.bracken_G ] || { say "FATAL bracken failed on $c/$s"; tail -2 $o.bracken.log; exit 1; }
    pct=$(grep -oE 'classified \([0-9.]+%' $o.kraken.log | head -1 | cut -d'(' -f2)
    say "$c/$s: $pct classified · $(ls $OUT/$c/*.bracken_G | wc -l)/${#R1[@]} done"
  done
  python3 $ROOT/scripts/08c_matrices.py $OUT/$c $c
done
say "--- adequacy test (Amendment 1): median classified fraction, Nelore, threshold 60% ---"
python3 - <<PY
import re, glob, statistics
v = [float(re.search(r'sequences classified \(([\d.]+)%\)', open(f).read()).group(1)) for f in sorted(glob.glob("$OUT/nelore/*.kraken.log"))]
print(f"n={len(v)} median={statistics.median(v):.1f}% range {min(v):.1f}-{max(v):.1f}% ->", "PASS" if statistics.median(v) >= 60 else "FAIL")
PY
say "PROFILING COMPLETE"
