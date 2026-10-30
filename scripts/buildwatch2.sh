#!/usr/bin/env bash
# Progress of a running build_db = bytes read from the library pipe / total .fna bytes. 30 s sample for rate and ETA.
pid=$(pgrep -x build_db | head -1)
if [ -z "$pid" ]; then echo "no build_db running"; tail -3 /data/refs/logs/kraken2_build.log | cut -c1-120; tail -2 /data/refs/logs/build.log; pgrep -af 'bracken|kraken2 ' | grep -v grep | cut -c1-100; exit 0; fi
tot=$(find /data/refs/k2_gtdb207/library/added -name '*.fna' -printf '%s\n' | awk '{s+=$1} END{print s}')
r1=$(awk '/^rchar/{print $2}' /proc/$pid/io); sleep 30; r2=$(awk '/^rchar/{print $2}' /proc/$pid/io)
awk -v a=$r2 -v b=$tot -v d=$((r2-r1)) -v e="$(ps -o etime= -p $pid | tr -d ' ')" 'BEGIN{rate=d/30; printf "build_db up %s · %.1f%% of %.0f GB read · %.0f MB/s · ETA %.1f h\n", e, 100*a/b, b/1e9, rate/1e6, (rate>0)?(b-a)/rate/3600:0}'
