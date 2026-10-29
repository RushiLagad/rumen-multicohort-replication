#!/usr/bin/env bash
# Downloads only. No tools required. Indexing happens tomorrow via 06_setup_references.sh
set -uo pipefail
REF=/data/refs
say () { echo "[$(date +%H:%M:%S)] $*"; }

get () {
  local out="$1" url="$2"
  [ -s "$out" ] && { say "$(basename "$out") present, $(( $(stat -c%s "$out")/1000000 )) MB"; return; }
  say "fetching $(basename "$out")"
  curl -sSL -C - --retry 8 --retry-delay 15 --retry-all-errors \
       --speed-limit 10000 --speed-time 180 -o "$out" "$url" \
    && say "done $(basename "$out"): $(( $(stat -c%s "$out")/1000000 )) MB" \
    || say "FAILED $(basename "$out")"
}

get "$REF/host/ARS-UCD2.0.fa.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz"
get "$REF/host/UOA_Brahman_1.fa.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/369/695/GCF_003369695.1_UOA_Brahman_1/GCF_003369695.1_UOA_Brahman_1_genomic.fna.gz"

RAM=$(free -g | awk 'NR==2{print $2}')
say "RAM: ${RAM} GB"
if   [ "$RAM" -ge 100 ]; then K2=k2_standard_20240605.tar.gz
elif [ "$RAM" -ge 40 ];  then K2=k2_standard_16gb_20240605.tar.gz
else                          K2=k2_standard_8gb_20240605.tar.gz; fi
say "kraken2 build: $K2"
get "$REF/kraken2/$K2" "https://genome-idx.s3.amazonaws.com/kraken/$K2"
[ -s "$REF/kraken2/$K2" ] && { say "extracting"; tar -xzf "$REF/kraken2/$K2" -C "$REF/kraken2" && rm -f "$REF/kraken2/$K2" && say "extracted"; }
say "ALL DOWNLOADS DONE"
