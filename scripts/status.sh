#!/usr/bin/env bash
cd /data/rumen_multicohort
echo "=== preprocessing (nelore) ==="
pgrep -f 07_preprocess >/dev/null && echo "  running, minimap2 x$(pgrep -c minimap2 || echo 0)" || echo "  NOT RUNNING"
echo "  done: $(ls -1 processed/nelore/*_clean_1.fq.gz 2>/dev/null | wc -l) of 52"
tail -2 processed/nelore_run.log 2>/dev/null | sed 's/^/  /'
echo "=== prefetch ==="
pgrep -f 09_prefetch_sra >/dev/null && echo "  running" || echo "  not running"
echo "  charolais: $(ls -1 PRJEB70895/sra/*.sra 2>/dev/null | wc -l) of 61"
echo "  li wgs   : $(ls -1 /data/rumen448333/sra/*.sra 2>/dev/null | wc -l) of 48"
tail -2 logs/prefetch_all.log 2>/dev/null | sed 's/^/  /'
echo "=== disk ==="
df -h /data | tail -1 | sed 's/^/  /'
