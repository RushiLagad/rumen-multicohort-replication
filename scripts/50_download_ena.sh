#!/bin/bash
# Amendment 19: download reads for new cohorts from ENA, resumable, MD5-checked. Phenotype files are not touched.
set -uo pipefail; R=${RUMEN_ROOT:-/data/rumen_multicohort}; cd $R; mkdir -p raw logs/download
one() { url=$1; md5=$2; out=$3; mkdir -p "$(dirname $out)"
  if [ -s "$out" ] && [ "$(md5sum < $out | cut -c1-32)" = "$md5" ]; then echo "skip $out"; return 0; fi
  for try in 1 2 3; do
    wget -q -c -T 120 -t 5 -O "$out.part" "https://$url" 2>/dev/null || true
    if [ -s "$out.part" ] && [ "$(md5sum < $out.part | cut -c1-32)" = "$md5" ]; then mv "$out.part" "$out"; echo "ok $out"; return 0; fi
    rm -f "$out.part"
  done; echo "FAIL $out"; return 1; }
export -f one
for p in "$@"; do
  awk -F'\t' -v p=$p 'NR>1{n=split($5,u,";"); split($7,m,";"); for(i=1;i<=n;i++){k=split(u[i],f,"/"); print u[i]"\t"m[i]"\traw/"p"/"f[k]}}' metadata/candidates/$p.files.tsv > logs/download/$p.list
  echo "$(date +%F_%T) $p: $(wc -l < logs/download/$p.list) files"
  tr '\t' ' ' < logs/download/$p.list | xargs -P ${JOBS:-8} -L1 bash -c 'one "$0" "$1" "$2"' > logs/download/$p.log 2>&1
  echo "$(date +%F_%T) $p: ok $(grep -c '^ok' logs/download/$p.log) skip $(grep -c '^skip' logs/download/$p.log) FAIL $(grep -c '^FAIL' logs/download/$p.log)"
done
touch logs/download/DONE
