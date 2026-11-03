#!/bin/bash
# Amendment 17 step 1: remove 341F/806R and the spacer before them. No phenotype data used.
set -euo pipefail; cd /data/rumen_multicohort
O=work/charolais16s; mkdir -p $O/trimmed logs/charolais16s results/charolais16s
mamba run -n cutadapt cutadapt --version > /dev/null 2>&1 || mamba create -y -q -n cutadapt -c conda-forge -c bioconda cutadapt > logs/charolais16s/install.log 2>&1
awk -F'\t' '$1=="Charolais_PRJEB70895" && $3=="16S_amplicon"{print $4}' harmonized/runs.tsv | sort > $O/runs.txt
[ "$(wc -l < $O/runs.txt)" = 30 ]
one() { r=$1; mamba run -n cutadapt cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACHVGGGTWTCTAAT -e 0.12 -O 12 --discard-untrimmed -m 100 \
  -o work/charolais16s/trimmed/${r}_1.fastq.gz -p work/charolais16s/trimmed/${r}_2.fastq.gz \
  PRJEB70895/fastq/${r}_1.fastq.gz PRJEB70895/fastq/${r}_2.fastq.gz > logs/charolais16s/cutadapt_$r.log 2>&1; }
export -f one; xargs -P 6 -I{} bash -c 'one {}' < $O/runs.txt
for r in $(cat $O/runs.txt); do awk -v r=$r -F: '/Total read pairs processed/{gsub(/[ ,]/,"",$2);a=$2} /Pairs written/{split($2,b," ");gsub(/,/,"",b[1]);w=b[1]} END{print r"\t"a"\t"w"\t"w/a}' logs/charolais16s/cutadapt_$r.log; done > results/charolais16s/trim_summary.tsv
sort -k4,4n results/charolais16s/trim_summary.tsv | sed -n '1p;15p;30p'; ls $O/trimmed | wc -l
