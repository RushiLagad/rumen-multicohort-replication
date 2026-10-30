#!/usr/bin/env bash
# Kraken2 database from Struo2's GTDB r207 library.
# Pre-built hash.k2d is 287 GB vs 122 GB RAM, so we build capped at 85 GB.
# Checkpointed: rerunning skips completed stages.
set -uo pipefail
SRC=/data/refs/gtdb207
DB=/data/refs/k2_gtdb207
THREADS="${THREADS:-24}"
MAXDB="${MAXDB:-85000000000}"
mkdir -p "$DB" /data/refs/logs
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

if [ ! -s "$DB/taxonomy/nodes.dmp" ]; then
  say "fetching taxdump.tar.gz"
  rm -rf "$SRC/taxdump"; mkdir -p "$SRC/taxdump" "$DB/taxonomy"
  curl -sSL --retry 8 -o "$SRC/taxdump/taxdump.tar.gz" \
    "http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release207/taxdump/taxdump.tar.gz"
  tar -xzf "$SRC/taxdump/taxdump.tar.gz" -C "$DB/taxonomy"
  if [ ! -s "$DB/taxonomy/nodes.dmp" ]; then
    n=$(find "$DB/taxonomy" -name nodes.dmp | head -1)
    [ -n "$n" ] && mv "$(dirname "$n")"/*.dmp "$DB/taxonomy/"
  fi
  say "nodes.dmp $(wc -l < "$DB/taxonomy/nodes.dmp") lines, names.dmp $(wc -l < "$DB/taxonomy/names.dmp") lines"
else
  say "taxonomy present"
fi
[ -s "$DB/taxonomy/nodes.dmp" ] || { say "FATAL: no nodes.dmp"; exit 1; }

if [ ! -f "$DB/.library_complete" ]; then
  say "extracting library.tar.gz, 62 GB compressed, 20-40 min"
  rm -rf "$DB/library"
  tar -xzf "$SRC/library.tar.gz" -C "$DB" && touch "$DB/.library_complete"
  [ -f "$DB/.library_complete" ] || { say "FATAL: extraction failed"; exit 1; }
  say "library: $(find "$DB/library" -name "*.fna" -printf . | wc -c) fna, $(du -sh "$DB/library" | cut -f1)"
else
  say "library present: $(du -sh "$DB/library" | cut -f1)"
fi
cp -n "$SRC/seqid2taxid.map" "$DB/seqid2taxid.map" 2>/dev/null || true
df -h /data | tail -1

if [ -s "$DB/hash.k2d" ]; then
  say "hash.k2d present, skipping build"
else
  n=$(find "$DB/library" -name "*.fna" -printf . 2>/dev/null | wc -c)
  [ "$n" -lt 1000 ] && { say "FATAL: only $n fna files in library, refusing to build"; exit 1; }
  say "building from $n genomes, cap $(( MAXDB/1000000000 )) GB, $THREADS threads, several hours"
  kraken2-build --build --db "$DB" --threads "$THREADS" --max-db-size "$MAXDB" > /data/refs/logs/kraken2_build.log 2>&1 || { say "FATAL kraken2-build failed, see kraken2_build.log"; exit 1; }
fi
[ -s "$DB/hash.k2d" ] || { say "FATAL: no hash.k2d produced"; exit 1; }
say "db files: $(du -sh "$DB"/*.k2d | tr '\n' ' ')"

if [ ! -s "$DB/database150mers.kmer_distrib" ]; then
  say "bracken distributions, read length 150"
  bracken-build -d "$DB" -t "$THREADS" -k 35 -l 150 > /data/refs/logs/bracken_build.log 2>&1 || { say "FATAL bracken-build failed, see bracken_build.log"; exit 1; }
fi
say "BUILD COMPLETE"
du -sh "$DB"; df -h /data | tail -1
