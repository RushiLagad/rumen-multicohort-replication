#!/bin/bash
# Amendment 21 step 1: remove both primers from a paired 16S cohort. Usage: 60_trim_16s_pe.sh <BioProject> <fwd primer> <rev primer>
set -euo pipefail; cd ${RUMEN_ROOT:-/data/rumen_multicohort}; export P=$1 F=$2 R=$3
O=work/16s/$P; mkdir -p $O/trimmed logs/16s/$P results/16s/$P
miss=$(cut -f3 logs/download/$P.list | while read f; do [ -s "$f" ] || echo "$f"; done | wc -l)
[ "$miss" = 0 ] || { echo "$P: $miss listed files not downloaded yet; rerun the downloader first"; exit 1; }
cut -f3 logs/download/$P.list | grep '_1\.fastq\.gz$' | sed 's#.*/##; s#_1\.fastq\.gz$##' | sort > $O/runs.txt
echo "$P: $(wc -l < $O/runs.txt) paired runs"
one() { r=$1; mamba run -n cutadapt cutadapt -g "$F" -G "$R" -e 0.12 -O 12 --discard-untrimmed -m 100 -j 1 \
  -o work/16s/$P/trimmed/${r}_1.fastq.gz -p work/16s/$P/trimmed/${r}_2.fastq.gz \
  raw/$P/${r}_1.fastq.gz raw/$P/${r}_2.fastq.gz > logs/16s/$P/cutadapt_$r.log 2>&1; }
export -f one; xargs -P 12 -I{} bash -c 'one {}' < $O/runs.txt
for r in $(cat $O/runs.txt); do awk -v r=$r -F: '/Total read pairs processed/{gsub(/[ ,]/,"",$2);a=$2} /Pairs written/{split($2,b," ");gsub(/,/,"",b[1]);w=b[1]} END{print r"\t"a"\t"w"\t"(a>0?w/a:0)}' logs/16s/$P/cutadapt_$r.log; done > results/16s/$P/trim_summary.tsv
awk -F'\t' '{a+=$2; w+=$3} END{printf "pairs in %d, kept %d (%.3f)\n", a, w, w/a}' results/16s/$P/trim_summary.tsv; sort -k4,4n results/16s/$P/trim_summary.tsv | sed -n '1p;$p'
