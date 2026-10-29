#!/usr/bin/env bash
# =============================================================================
# 07b_reference_choice.sh   ANALYSIS_PLAN Section 5 + Amendment 2
# Compare host depletion against ARS-UCD2.0 (Bos taurus) and UOA_Brahman_1
# (Bos indicus) on 8 Nelore samples. Adopt Bos indicus if the retained-read
# difference exceeds 2%. Subset rather than all 52 per Amendment 2.
# =============================================================================
set -uo pipefail
ROOT=/data/rumen_multicohort; REF=/data/refs/host
FQ="$ROOT/PRJNA987743/fastq"; WORK="$ROOT/PRJNA987743/refchoice"
mkdir -p "$WORK"; T="${THREADS:-24}"; N="${N_SAMPLES:-8}"
: > "$WORK/summary.tsv"

mapfile -t RUNS < <(ls -1 "$FQ"/*_1.fastq.gz 2>/dev/null | xargs -rn1 basename \
                    | sed 's/_1.fastq.gz//' | sort | head -"$N")
echo "=== reference choice on ${#RUNS[@]} samples ==="
[ "${#RUNS[@]}" -lt 3 ] && { echo "need at least 3 complete samples, have ${#RUNS[@]}"; exit 1; }

printf "%-14s %10s %12s %12s %9s\n" run pairs taurus_kept brahman_kept diff_pct
for r in "${RUNS[@]}"; do
  R1="$FQ/${r}_1.fastq.gz"; R2="$FQ/${r}_2.fastq.gz"
  gzip -t "$R1" 2>/dev/null && gzip -t "$R2" 2>/dev/null || { echo "  $r incomplete, skip"; continue; }
  S1="$WORK/${r}_s1.fq"; S2="$WORK/${r}_s2.fq"
  zcat "$R1" | head -1600000 > "$S1"; zcat "$R2" | head -1600000 > "$S2"
  tot=$(( $(wc -l < "$S1") / 4 ))
  for ref in ARS-UCD2.0 UOA_Brahman_1; do
    minimap2 -ax sr -t "$T" "$REF/${ref}.mmi" "$S1" "$S2" 2>/dev/null \
      | samtools view -c -f 12 -F 256 - > "$WORK/${r}.${ref}.n"
  done
  bt=$(( $(cat "$WORK/${r}.ARS-UCD2.0.n")/2 )); bi=$(( $(cat "$WORK/${r}.UOA_Brahman_1.n")/2 ))
  d=$(python3 -c "print(f'{abs($bi-$bt)/max($tot,1)*100:.3f}')")
  printf "%-14s %10d %12d %12d %9s\n" "$r" "$tot" "$bt" "$bi" "$d"
  echo -e "$r\t$tot\t$bt\t$bi\t$d" >> "$WORK/summary.tsv"
  rm -f "$S1" "$S2"
done

python3 -I - "$WORK/summary.tsv" << 'PY'
import sys, csv, statistics as st
rows=[r for r in csv.reader(open(sys.argv[1]),delimiter="\t") if len(r)==5]
if not rows: raise SystemExit("no results")
d=[float(r[4]) for r in rows]
hb=[100*int(r[2])/int(r[1]) for r in rows]
print(f"\n  retained-read difference: mean {st.mean(d):.3f}%  max {max(d):.3f}%  n={len(d)}")
print(f"  host fraction (ARS-UCD2.0): median {100-st.median(hb):.1f}%")
if st.mean(d) > 2.0:
    print("\n  DECISION: >2%. Use UOA_Brahman_1 for Nelore.")
    print("  export HOST_REF_NELORE=UOA_Brahman_1  before running 07_preprocess.sh")
else:
    print("\n  DECISION: <2%. Use ARS-UCD2.0 for all cohorts.")
    print("  One reference across cohorts removes a between-cohort technical factor.")
print("  Record in ANALYSIS_PLAN.md Section 13.")
PY
