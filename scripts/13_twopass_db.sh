#!/usr/bin/env bash
# Amendment 9 two-pass: full-resolution GTDB r207 DB from genomes seen in the capped screen, then re-profile
set -uo pipefail
OLD=${OLD:-/data/refs/k2_gtdb207}; NEW=${NEW:-/data/refs/k2_gtdb207_full}; SRC=${SRC:-/data/refs/gtdb207}
ROOT=${ROOT:-/data/rumen_multicohort}; LOG=${LOG:-/data/refs/logs/kraken2_build_full.log}; T=24
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
die(){ say "FATAL: $*"; exit 1; }
for x in kraken2 kraken2-build kmer2read_distr generate_kmer_distribution.py pigz python3; do command -v $x >/dev/null || die "missing $x"; done
mkdir -p $NEW
if [ ! -s $NEW/selected_genomes.tsv ]; then
python3 - <<PY || die "selection"
import os
files=[os.path.join(d,f) for d,_,fs in os.walk("$ROOT/profiles_gtdb") for f in fs if "report" in f and "bracken" not in f]
print("kraken reports found:",len(files),"(expect 131)")
if len(files)!=131: raise SystemExit("report count != 131")
sp=set()
for f in files:
    for l in open(f):
        p=l.rstrip("\n").split("\t")
        if len(p)>=6 and p[3].strip()=="S" and int(p[1])>=100: sp.add(p[4].strip())
par={}
for l in open("$OLD/taxonomy/nodes.dmp"):
    p=l.split("\t|\t"); par[p[0].strip()]=p[1].strip()
lv=set(l.rstrip("\n").split("\t")[1] for l in open("$OLD/seqid2taxid.map"))
sel=sorted((g for g in lv if g in sp or par.get(g) in sp),key=int)
print(f"species with >=100 clade reads: {len(sp)}; genomes under them: {len(sel)} of {len(lv)}")
if not sel: raise SystemExit("empty selection")
open("$NEW/selected_genomes.tsv","w").write("\n".join(sel)+"\n")
PY
fi
n=$(wc -l < $NEW/selected_genomes.tsv); say "selection: $n genomes"
if [ ! -f $NEW/.library_pruned ]; then
  m=$( (set +o pipefail; pigz -dc $SRC/library.tar.gz | tar -t | head -3 | tr '\n' ' ') 2>/dev/null )
  h=$( (set +o pipefail; pigz -dc $SRC/library.tar.gz | tar -x -O --wildcards '*.fna' | head -1) 2>/dev/null )
  say "tar members: $m"; say "first header: ${h:0:100}"
  t=${h#>}; t=${t%%[[:space:]]*}; t2=${t#kraken:taxid|}; t2=${t2#*|}
  awk -v a="$t" -v b="$t2" '$1==a||$1==b{f=1;exit} END{exit !f}' $OLD/seqid2taxid.map || die "first header not in seqid2taxid.map"
  rm -rf $NEW/library; mkdir -p $NEW/library; say "extracting (30-40 min)"
  pigz -dc $SRC/library.tar.gz | tar -x -C $NEW --wildcards '*.fna' || die "extract"
  say "extracted $(find $NEW/library -name '*.fna' | wc -l) fna; pruning"
python3 - <<PY || die "prune gate"
import os,sys
sel=set(l.strip() for l in open("$NEW/selected_genomes.tsv"))
m=dict(l.rstrip("\n").split("\t")[:2] for l in open("$OLD/seqid2taxid.map"))
kept=drop=0
for d,_,fs in os.walk("$NEW/library"):
    for f in fs:
        if not f.endswith(".fna"): continue
        p=os.path.join(d,f)
        t=open(p).readline()[1:].split()[0]
        x=m.get(t) or m.get(t.split("|",2)[-1])
        if x in sel: kept+=1
        else: os.remove(p); drop+=1
print(f"kept {kept} genomes, removed {drop}; selected {len(sel)}")
if kept!=len(sel): sys.exit("kept != selected")
PY
  touch $NEW/.library_pruned
fi
awk -F'\t' 'NR==FNR{s[$1]=1;next} ($2 in s)' $NEW/selected_genomes.tsv $OLD/seqid2taxid.map > $NEW/seqid2taxid.map
rm -rf $NEW/taxonomy; cp -r $OLD/taxonomy $NEW/
exp=$(wc -l < $NEW/seqid2taxid.map); say "map: $exp contigs, taxonomy copied"
if [ ! -s $NEW/hash.k2d ]; then
  say "kraken2-build, full minimizers, 96 GB RAM guard"
  kraken2-build --build --db $NEW --threads $T --max-db-size 96000000000 > $LOG 2>&1 || die "kraken2-build, see $LOG"
fi
grep -iE 'estimated|lower maximum' $LOG | head -3 | while read -r l; do say "build: $l"; done
got=$(grep -oP 'Completed processing of \K[0-9]+' $LOG | tail -1)
[ "${got:-0}" -eq "$exp" ] || die "build processed ${got:-0} of $exp contigs"
sz=$(stat -c %s $NEW/hash.k2d); say "VERIFY PASS: $got/$exp contigs, hash.k2d $((sz/1000000000)) GB"
[ "$sz" -ge 95900000000 ] && say "NOTE: table CAPPED at 96 GB guard -> write Amendment 11" || say "table uncapped (full minimizer set)"
cd $NEW
if [ ! -f .dbkraken_ok ]; then
  say "classifying library for Bracken"
  kraken2 --db $NEW --threads $T <(find $NEW/library -name '*.fna' -print0 | xargs -0 cat) > database.kraken 2> database.kraken.log || die "library classification"
  [ "$(wc -l < database.kraken)" -eq "$exp" ] || die "database.kraken lines != $exp"
  touch .dbkraken_ok
fi
for L in 136 150 100; do
  [ -s database${L}mers.kmer_distrib ] && continue
  kmer2read_distr --seqid2taxid seqid2taxid.map --taxonomy taxonomy --kraken database.kraken --output database${L}mers.kraken -k 35 -l $L -t $T > kmer2read_$L.log 2>&1 || die "kmer2read $L"
  generate_kmer_distribution.py -i database${L}mers.kraken -o database${L}mers.kmer_distrib > distrib_$L.log 2>&1 || die "distrib $L"
  say "bracken distribution $L bp done"
done
say "profiling 131 samples on full-resolution DB"
bash $ROOT/scripts/08d_profile_gtdb_full.sh
say "TWO-PASS COMPLETE"
