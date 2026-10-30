#!/usr/bin/env bash
# Verify /data/refs/k2_gtdb207 after 10_build_gtdb_db.sh. Count gate always; heavy checks only when no kraken2 is running.
set -uo pipefail
DB=/data/refs/k2_gtdb207; LOG=/data/refs/logs/kraken2_build.log; fail=0
expected=$(wc -l < $DB/seqid2taxid.map)
got=$(grep -oP 'Completed processing of \K[0-9]+' $LOG | tail -1)
bp=$(grep -oP 'sequences, \K[0-9]+(?= bp)' $LOG | tail -1)
echo "sequences processed: ${got:-NONE} / $expected expected · bp: ${bp:-NONE} (expect ~2.0e11)"
[ "${got:-0}" -eq "$expected" ] || { echo "FAIL: sequence count"; fail=1; }
for f in hash.k2d opts.k2d taxo.k2d; do [ -s $DB/$f ] && echo "ok $f $(du -h $DB/$f | cut -f1)" || { echo "FAIL: $f missing"; fail=1; }; done
[ "$(stat -c %s $DB/hash.k2d 2>/dev/null || echo 0)" -gt 80000000000 ] || { echo "FAIL: hash.k2d too small"; fail=1; }
if pgrep -x kraken2 >/dev/null || pgrep -f bracken-build >/dev/null; then
  echo "kraken2/bracken running; skipping inspect + smoke test, rerun later"
else
  echo "--- kraken2-inspect ---"; kraken2-inspect --db $DB --threads 8 --memory-mapping > /tmp/k2_inspect.txt 2>/dev/null
  echo "taxa: $(wc -l < /tmp/k2_inspect.txt) · genera: $(awk '$4=="G"' /tmp/k2_inspect.txt | wc -l)"
  for g in Prevotella Ruminococcus Fibrobacter Butyrivibrio Methanobrevibacter Succiniclasticum Selenomonas Treponema; do
    awk -v g=$g '$4=="G" && $0 ~ g {c++; if(c<=3) printf "  %s (%s minimizers)\n", $6, $3} END{if(!c) print "  "g": ABSENT"}' /tmp/k2_inspect.txt; done
  echo "--- smoke: 100k pairs, one Nelore sample (not the pre-registered test) ---"
  r1=$(ls /data/rumen_multicohort/processed/nelore/* 2>/dev/null | grep -E '_1\.(fastq|fq)(\.gz)?$' | head -1); r2=${r1/_1./_2.}
  if [ -s "$r1" ] && [ -s "$r2" ]; then zcat -f $r1 | head -400000 > /tmp/smoke_1.fq; zcat -f $r2 | head -400000 > /tmp/smoke_2.fq
    kraken2 --db $DB --threads 24 --memory-mapping --paired /tmp/smoke_1.fq /tmp/smoke_2.fq --report /tmp/smoke.report --output /dev/null 2>&1 | grep -E 'classified'
  else echo "nelore fastq not found for smoke test"; fi
fi
[ $fail -eq 0 ] && echo "VERIFY PASS" || { echo "VERIFY FAIL"; exit 1; }
