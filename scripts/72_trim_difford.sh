#!/bin/bash
# Amendment 28 items 1-3: primer removal for one Difford marker. Usage: 72_trim_difford.sh <bact|arch>   (LIMIT=n for a test)
set -euo pipefail; cd ${RUMEN_ROOT:-/data/rumen_multicohort}; M=$1; P=PRJEB28065
case $M in bact) F=AGAGTTTGATCCTGGCTCAG; R=ATTACCGCGGCTGCTGG;; arch) F=CAGCMGCCGCGGTAA; R=GGCCATGCACCWCCTCTC;; *) echo "marker must be bact or arch"; exit 1;; esac
export P M F R; O=work/16s/difford_$M; mkdir -p $O/trimmed logs/16s/difford_$M results/16s/difford_$M
awk -F'\t' -v m=$M 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["in_s4"]=="True" && $c["marker"]==m {print $c["run"]}' harmonized/difford_runs.tsv | sort > $O/runs.txt
[ -n "${LIMIT:-}" ] && { head -n $LIMIT $O/runs.txt > $O/r.tmp; mv $O/r.tmp $O/runs.txt; }
echo "difford_$M: $(wc -l < $O/runs.txt) paired runs"
one() { r=$1; /home/rlagad/.local/share/mamba/envs/cutadapt/bin/cutadapt -g "$F" -G "$R" -e 0.12 -O 12 --discard-untrimmed -m 100 -j 1 \
  -o work/16s/difford_$M/trimmed/${r}_1.fastq.gz -p work/16s/difford_$M/trimmed/${r}_2.fastq.gz \
  raw/$P/${r}_1.fastq.gz raw/$P/${r}_2.fastq.gz > logs/16s/difford_$M/cutadapt_$r.log 2>&1; }
export -f one; xargs -P 12 -I{} bash -c 'one {}' < $O/runs.txt
mkdir -p $O/trimmed_sub; xargs -P 12 -I{} python3 -I scripts/75_subsample_pairs.py $O/trimmed/{}_1.fastq.gz $O/trimmed/{}_2.fastq.gz $O/trimmed_sub/{}_1.fastq.gz $O/trimmed_sub/{}_2.fastq.gz 50000 {} < $O/runs.txt > results/16s/difford_$M/subsample.tsv
for r in $(cat $O/runs.txt); do awk -v r=$r -F: '/Total read pairs processed/{gsub(/[ ,]/,"",$2);a=$2} /Pairs written/{split($2,b," ");gsub(/,/,"",b[1]);w=b[1]} END{print r"\t"a"\t"w"\t"(a>0?w/a:0)}' logs/16s/difford_$M/cutadapt_$r.log; done > results/16s/difford_$M/trim_summary.tsv
awk -F'\t' '{a+=$2; w+=$3} END{printf "pairs in %d, kept %d (%.3f)\n", a, w, w/a}' results/16s/difford_$M/trim_summary.tsv; sort -k4,4n results/16s/difford_$M/trim_summary.tsv | sed -n '1p;$p'
