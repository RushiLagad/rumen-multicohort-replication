#!/usr/bin/env bash
# One sample: kraken2 + Bracken S and G. Usage: 08e_one_sample.sh <cohort> <R1> <read_len>
set -uo pipefail
c=$1; r1=$2; L=$3
DB=${DB:-/data/refs/k2_gtdb207_full}; ROOT=${ROOT:-/data/rumen_multicohort}; OUT=${OUT:-$ROOT/profiles_gtdb_full}; T=${T:-8}
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
r2=${r1/_1./_2.}; s=$(basename $r1 | sed -E 's/(_1)?\.(fastq|fq)(\.gz)?$//'); o=$OUT/$c/$s
[ -s $o.bracken_G ] && exit 0
rm -f $o.kreport $o.bracken_S $o.bracken_G
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
say "$c/$s: $pct classified · $(ls $OUT/$c/*.bracken_G | wc -l) done"
