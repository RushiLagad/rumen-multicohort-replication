#!/bin/bash
# Amendment 19 support: fetch runs that failed on ENA from NCBI SRA instead. Usage: 51_ncbi_fallback.sh <BioProject>. File layout from ENA's file list; checked by read count against ENA.
set -uo pipefail; cd ${RUMEN_ROOT:-/data/rumen_multicohort}; export P=$1; mkdir -p work/sra logs/download
grep '^FAIL' logs/download/$P.log | sed 's#.*/##; s#_[12]\.fastq\.gz$##; s#\.fastq\.gz$##' | sort -u > logs/download/$P.ncbi_runs
echo "$(date +%F_%T) $P: $(wc -l < logs/download/$P.ncbi_runs) runs to fetch from NCBI"
one() { r=$1; nf=$(awk -F'\t' -v r=$r '$1==r{print split($5, u, ";")}' metadata/candidates/$P.files.tsv); want=$(awk -F'\t' -v r=$r '$1==r{print $4}' metadata/candidates/$P.files.tsv)
  t=work/sra/$r
  if ! ls $t/*.fastq > /dev/null 2>&1; then rm -rf $t; mkdir -p $t
    prefetch -q -O $t $r > $t/log 2>&1 && fasterq-dump -q --split-files -e 2 -O $t $t/$r >> $t/log 2>&1 || { echo "FAIL $r ncbi"; return 1; }; fi
  if [ "$nf" = 2 ]; then fs="${r}_1 ${r}_2"; else fs=$r; [ -s $t/$r.fastq ] || mv $t/${r}_1.fastq $t/$r.fastq 2>/dev/null; fi
  for f in $fs; do [ -s $t/$f.fastq ] || { echo "FAIL $r missing $f"; return 1; }; done
  got=$(( $(wc -l < $t/$(echo $fs | cut -d' ' -f1).fastq) / 4 ))
  [ "$got" = "$want" ] || { echo "FAIL $r reads $got want $want"; return 1; }
  for f in $fs; do gzip -c $t/$f.fastq > raw/$P/$f.fastq.gz.tmp && mv raw/$P/$f.fastq.gz.tmp raw/$P/$f.fastq.gz; done
  rm -rf $t; echo "ok $r $got reads"; }
export -f one; xargs -P ${JOBS:-6} -I{} bash -c 'one {}' < logs/download/$P.ncbi_runs > logs/download/$P.ncbi.log 2>&1
awk -v p=$P 'BEGIN{OFS="\t"} /^ok/{print $2, "NCBI SRA", $3}' logs/download/$P.ncbi.log >> metadata/candidates/$P.ncbi_sourced.tsv
echo "$(date +%F_%T) $P: ok $(grep -c '^ok' logs/download/$P.ncbi.log) FAIL $(grep -c '^FAIL' logs/download/$P.ncbi.log)"
