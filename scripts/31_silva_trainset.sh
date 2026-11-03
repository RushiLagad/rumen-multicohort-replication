#!/bin/bash
# Amendment 17: convert the local mothur-format SILVA 138.1 nr reference to a DADA2 training file (6 ranks).
set -euo pipefail; D=/data/references/silva; O=/data/rumen_multicohort/work/charolais16s; mkdir -p $O
awk -F'\t' 'FNR==NR{t[$1]=$2;next} /^>/{if(h!=""&&s!="")print h"\n"s; id=substr($1,2); h=(id in t)?">"t[id]:""; s=""; next} h!=""{x=toupper($0); gsub(/[.-]/,"",x); gsub(/U/,"T",x); s=s x} END{if(h!=""&&s!="")print h"\n"s}' $D/silva.nr_v138_1.tax $D/silva.nr_v138_1.align > $O/silva138_train.fa
grep -c '^>' $O/silva138_train.fa; grep '^>' $O/silva138_train.fa | cut -d';' -f1 | sort | uniq -c
grep -m2 -A1 '^>' $O/silva138_train.fa | cut -c1-90
awk '!/^>/{n++; l+=length($0)} END{print "mean length", l/n}' $O/silva138_train.fa
grep -c '_ge;$' $O/silva138_train.fa
