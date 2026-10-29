#!/usr/bin/env bash
# =============================================================================
# 08_profile.sh   Kraken2 + Bracken -> genus and species count matrices
#   bash 08_profile.sh nelore [check]
# Genus is the PRIMARY feature space (Section 4). Species is sensitivity only.
# Runs the Amendment 1 adequacy test: median classified fraction must be >=60%.
# Samples run sequentially inside one process so the 78 GB database loads once.
# =============================================================================
set -uo pipefail
ROOT=/data/rumen_multicohort; DB=/data/refs/kraken2
COHORT="${1:-}"; MODE="${2:-run}"
[ -n "$COHORT" ] || { echo "usage: $0 {nelore|charolais_mrna|li} [check]"; exit 1; }
IN="$ROOT/processed/$COHORT"; OUT="$ROOT/profiles/$COHORT"
mkdir -p "$OUT"/{kraken,bracken,logs}
T="${THREADS:-24}"; READLEN="${READLEN:-150}"

mapfile -t S < <(ls -1 "$IN"/*_clean_1.fq.gz 2>/dev/null | xargs -rn1 basename | sed 's/_clean_1.fq.gz//')
echo "=== profiling $COHORT ==="
echo "  ready : ${#S[@]}   done : $(ls -1 "$OUT/bracken"/*.G.bracken 2>/dev/null | wc -l)"
[ "$MODE" = "check" ] && exit 0
[ -s "$DB/hash.k2d" ] || { echo "no kraken2 db"; exit 1; }

for s in "${S[@]}"; do
  rep="$OUT/kraken/${s}.report"
  if [ ! -s "$rep" ]; then
    echo "[$(date +%H:%M:%S)] kraken2 $s"
    kraken2 --db "$DB" --threads "$T" --paired --gzip-compressed \
      --report "$rep" --output /dev/null \
      "$IN/${s}_clean_1.fq.gz" "$IN/${s}_clean_2.fq.gz" \
      > "$OUT/logs/${s}.kraken2.log" 2>&1 || { echo "  FAILED $s"; rm -f "$rep"; continue; }
  fi
  for lvl in G S; do
    b="$OUT/bracken/${s}.${lvl}.bracken"
    [ -s "$b" ] || bracken -d "$DB" -i "$rep" -o "$b" -r "$READLEN" -l "$lvl" -t 10 \
      >> "$OUT/logs/${s}.bracken.log" 2>&1 || echo "  bracken $lvl failed $s"
  done
done

echo; echo "=== count matrices + Amendment 1 adequacy test ==="
python3 -I - "$OUT" "$COHORT" << 'PY'
import sys, os, glob, csv, collections, statistics as st
out, cohort = sys.argv[1], sys.argv[2]

frac=[]
for f in sorted(glob.glob(f"{out}/kraken/*.report")):
    lines=[l.split("\t") for l in open(f)]
    unc=next((float(l[0]) for l in lines if len(l)>5 and l[5].strip()=="unclassified"), 0.0)
    frac.append(100-unc)
if frac:
    m=st.median(frac)
    print(f"  classified-read fraction: median {m:.1f}%  range {min(frac):.1f}-{max(frac):.1f}%")
    print(f"  ADEQUACY TEST (Amendment 1, threshold 60%): {'PASS' if m>=60 else 'FAIL'}")
    if m<60:
        print("  -> database judged inadequate. Amendment 1 requires a rumen-inclusive")
        print("     build BEFORE any association testing. Do not proceed to analysis.")

for lvl,label in (("G","genus"),("S","species")):
    files=sorted(glob.glob(f"{out}/bracken/*.{lvl}.bracken"))
    if not files: continue
    counts=collections.defaultdict(dict); samples=[]
    for f in files:
        s=os.path.basename(f).split(f".{lvl}.")[0]; samples.append(s)
        for row in csv.DictReader(open(f),delimiter="\t"):
            counts[row["name"]][s]=int(row["new_est_reads"])
    taxa=sorted(counts); m=f"{out}/{cohort}_{label}_counts.tsv"
    with open(m,"w") as fh:
        fh.write("taxon\t"+"\t".join(samples)+"\n")
        for t in taxa:
            fh.write(t+"\t"+"\t".join(str(counts[t].get(s,0)) for s in samples)+"\n")
    tot=sum(sum(counts[t].values()) for t in taxa)
    print(f"\n  {label}: {len(taxa)} taxa x {len(samples)} samples -> {os.path.basename(m)}")
    for g in ("Prevotella","Bacteroides","Fibrobacter","Faecalibacterium",
              "Methanobrevibacter","Ruminococcus","Butyrivibrio"):
        r=sum(sum(counts[t].values()) for t in taxa if t.startswith(g))
        print(f"    {g:20s} {100*r/max(tot,1):6.2f}% of classified reads")
PY
echo
echo "Sanity: in rumen content Prevotella should far exceed Bacteroides and"
echo "Fibrobacter should exceed Faecalibacterium. If reversed, the database is"
echo "forcing reads onto human-gut relatives and Amendment 1 applies."
