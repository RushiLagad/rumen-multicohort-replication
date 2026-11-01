#!/usr/bin/env bash
# Bracken distributions for k2_gtdb207, bypassing bracken-build's 62,443-argument kraken2 command line (ARG_MAX).
# The same three steps bracken-build runs: classify library -> kmer2read_distr -> generate_kmer_distribution.
set -uo pipefail
DB=/data/refs/k2_gtdb207; T=24; K=35; L=150
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
cd $DB
if [ -s database.kraken ] && [ ! -f database.kraken.incomplete ]; then say "database.kraken present, skipping classification"; else
  rm -f database.kraken database.kraken.tmp; touch database.kraken.incomplete
  say "classifying library (7,745,746 contigs, 202 Gbp) against the DB"
  kraken2 --db $DB --threads $T <(find $DB/library/added -name '*.fna' -print0 | xargs -0 cat) > database.kraken 2> database.kraken.log || { say "FATAL kraken2 classification failed"; tail -3 database.kraken.log; exit 1; }
  rm -f database.kraken.incomplete; say "classified: $(wc -l < database.kraken) lines · $(grep -E 'classified' database.kraken.log | tr '\n' ' ')"
fi
n=$(wc -l < database.kraken); [ "$n" -eq 7745746 ] || { say "FATAL database.kraken has $n lines, expected 7745746"; exit 1; }
say "kmer2read_distr k=$K l=$L (the slow step)"
kmer2read_distr --seqid2taxid $DB/seqid2taxid.map --taxonomy $DB/taxonomy --kraken $DB/database.kraken --output $DB/database${L}mers.kraken -k $K -l $L -t $T > kmer2read.log 2>&1 || { say "FATAL kmer2read_distr failed"; tail -5 kmer2read.log; exit 1; }
say "generate_kmer_distribution"
generate_kmer_distribution.py -i $DB/database${L}mers.kraken -o $DB/database${L}mers.kmer_distrib > kmer_distrib.log 2>&1 || { say "FATAL generate_kmer_distribution failed"; tail -5 kmer_distrib.log; exit 1; }
say "BRACKEN COMPLETE: $(wc -l < database${L}mers.kmer_distrib) rows"; ls -lh database${L}mers.*
