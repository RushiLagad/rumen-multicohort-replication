#!/usr/bin/env bash
# 147: clean-environment reproduction check (draft 20 item 3, requested by Dr. Goswami). Reporting only: nothing in the working repository is read or changed.
# Usage: bash scripts/147_reproduction_check.sh [WORK]      WORK (default /data/repro_v104) must not exist. Run under systemd-run with PATH and SSH_AUTH_SOCK passed (see docs/REPRODUCTION_CHECK.md).
#  1. fresh clone of the GitHub repository (branch master, or REF) into WORK/repo; the four data archives downloaded from the published Zenodo record
#     (ZENODO_RECORD, default 24327001 = v1.0.3, whose data archives are the ones v1.0.4 carries), MD5-checked against the record and unpacked at the clone root;
#  2. two new conda environments built only from the clone's env/rumen-r.yml and env/reporting.yml;
#  3. regenerated from results/ and compared with the committed files: number manifests (126, 127, 143, 145, 146), Figures 1 to 5 (140 + rsvg-convert, 139, 103),
#     Additional files 1 to 4 (142);
#  4. end to end from the processed inputs: Amendment 46 model-input dumps (21, 65, 81), per-cohort CLR estimates with Freedman-Lane p (112), the 48.1 seven-claim
#     table and Table 4 patterns (128, written to WORK/out/late_a48 through A48_OUT); compared with the archived outputs;
#  5. permutation-null SD, B and b/B of every rescaled W recomputed from the stored null fits (147_compare.py perm);
#  6. WORK/REPRODUCTION_CHECK.md written from the record (copy it to docs/ and commit).
set -uo pipefail
WORK=${1:-/data/repro_v104}; [ ! -e "$WORK" ] || { echo "$WORK exists; give a new folder"; exit 1; }
URL=${REPO_URL:-git@github.com:RushiLagad/rumen-multicohort-replication.git}; REC=${ZENODO_RECORD:-24327001}
mkdir -p "$WORK"/{archived,dl,logs,out}; exec > >(tee -a "$WORK/logs/147.log") 2>&1
R="$WORK/repo"; CMP="$WORK/out/compare.tsv"; printf 'step\trc\tseconds\n' > "$WORK/out/steps.tsv"; : > "$WORK/out/commands.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LANG=C.UTF-8 LC_ALL=C.UTF-8
step() { local n=$1; shift; local t0=$(date +%s); echo "$*" >> "$WORK/out/commands.txt"; echo "$(date '+%F %T') $n"
  ( cd "$R" 2>/dev/null || cd "$WORK"; bash -c "$*" ) > "$WORK/logs/$n.log" 2>&1; local rc=$?
  printf '%s\t%s\t%s\n' "$n" "$rc" "$(( $(date +%s) - t0 ))" >> "$WORK/out/steps.tsv"; [ $rc -eq 0 ] || echo "   $n exit $rc (see $WORK/logs/$n.log)"; return $rc; }
aside() { for p in "$@"; do [ -e "$R/$p" ] || { echo "   nothing to set aside at $p"; continue; }; mkdir -p "$WORK/archived/$(dirname "$p")"; mv "$R/$p" "$WORK/archived/$p"; done; }
cmpf() { "$PY" "$R/scripts/147_compare.py" pair "$1" "$2" "$3" "$4" "$CMP"; }          # new archived label group
cmpa() { cmpf "$R/$1" "$WORK/archived/$1" "$1" "$2"; }                                  # regenerated in place vs set aside

# ---- 1. clone
CM=""; for t in mamba micromamba conda; do command -v $t >/dev/null && { CM=$(command -v $t); break; }; done; [ -n "$CM" ] || { echo "no mamba, micromamba or conda"; exit 1; }
step clone "git clone -q --branch master '$URL' '$R' && cd '$R' && { [ -z '${REF:-}' ] || git checkout -q '${REF:-}'; }" || exit 1
COMMIT=$(git -C "$R" rev-parse --short=8 HEAD); echo "   clone at $COMMIT"
for f in env/rumen-r.yml env/reporting.yml scripts/147_compare.py; do [ -s "$R/$f" ] || { echo "missing $f in the clone"; exit 1; }; done

# ---- 2. environments, only from the clone's env files
mkenv() { case "$CM" in *micromamba) echo "'$CM' create -y -q -p '$1' -f '$2'";; *) echo "'$CM' env create -y -q -p '$1' -f '$2'";; esac; }
step env_rumen_r "$(mkenv "$WORK/envs/rumen-r" "$R/env/rumen-r.yml")" || exit 1
step env_reporting "$(mkenv "$WORK/envs/reporting" "$R/env/reporting.yml")" || exit 1
RR="$WORK/envs/rumen-r/bin/Rscript"; PY="$WORK/envs/reporting/bin/python"; export PATH="$WORK/envs/reporting/bin:$PATH"
for x in "$RR" "$PY" "$WORK/envs/reporting/bin/pandoc" "$WORK/envs/reporting/bin/rsvg-convert" "$WORK/envs/reporting/bin/pdftotext"; do [ -x "$x" ] || { echo "missing $x"; exit 1; }; done

# ---- 1b. data archives from the published record
step zenodo_download "'$PY' - '$REC' '$WORK/dl' <<'EOF'
import sys, json, hashlib, urllib.request
rec, d = sys.argv[1], sys.argv[2]
r = json.load(urllib.request.urlopen(urllib.request.Request(f'https://zenodo.org/api/records/{rec}', headers={'Accept': 'application/json'})))
fs = [f for f in r['files'] if f['key'].endswith('.tar.gz') and not f['key'].startswith('rumen-multicohort-replication-')]
assert len(fs) == 4, [f['key'] for f in fs]
for f in fs:
    p = f'{d}/{f[\"key\"]}'; urllib.request.urlretrieve(f['links']['self'], p); h = hashlib.md5(open(p, 'rb').read()).hexdigest()
    ok = 'md5:' + h == f['checksum']; print(f['key'], h, 'matches the record' if ok else 'DOES NOT MATCH ' + f['checksum']); assert ok
EOF" || exit 1
ARCH=$(awk '/matches the record/{printf "%s%s", s, $1; s=", "}' "$WORK/logs/zenodo_download.log")
step unpack "for a in '$WORK'/dl/*.tar.gz; do tar -xzf \"\$a\" -C '$R' || exit 1; done; ls -d profiles_gtdb_full/nelore results/scan_null results/late_a46/ancombc2_fl" || exit 1

# ---- 3. reporting outputs regenerated from results/ and compared
aside results/manifest/table_numbers.tsv results/manifest/draft17_numbers.tsv results/manifest/draft18_numbers.tsv results/manifest/draft19_numbers.tsv \
      results/manifest/draft20_numbers.tsv results/manifest/draft21_numbers.tsv results/late_p_compare results/figures_draft18 results/figures_manuscript results/additional_files_draft18
mkdir -p "$R/results/figures_draft18"; cp "$WORK/archived/results/figures_draft18/README.md" "$R/results/figures_draft18/"
step manifest_126 "RUMEN_ROOT='$R' '$PY' scripts/126_manifest_tables.py"
step manifest_127 "RUMEN_ROOT='$R' '$PY' scripts/127_manifest_draft17.py"
step manifest_143 "'$PY' scripts/143_manifest_draft18.py '$R'"
step manifest_145 "'$PY' scripts/145_manifest_draft19.py '$R'"
step p_compare_146 "'$PY' scripts/146_p_compare.py '$R'"
step figure1_140 "'$PY' scripts/140_figure1_compact.py '$R' && cd results/figures_draft18 && rsvg-convert -d 300 -p 300 -b white -f png -o fig1_study_design.png fig1_study_design.svg && rsvg-convert -b white -f pdf -o fig1_study_design.pdf fig1_study_design.svg"
step figures_139 "RUMEN_ROOT='$R' '$RR' scripts/139_figures_draft18.R"
step figures_103 "RUMEN_ROOT='$R' '$RR' scripts/103_figures_2to4.R"
step additional_files_142 "'$PY' scripts/142_additional_files_draft18.py '$R'"
G="Tables 1 to 5: number manifests (each script stops if a printed value differs from its result file)"
for f in table_numbers draft17_numbers draft18_numbers draft19_numbers draft20_numbers draft21_numbers; do cmpa results/manifest/$f.tsv "$G"; done
for f in seven_claims_p all_targets_p summary; do cmpa results/late_p_compare/$f.tsv "$G"; done
G="Figures 1 to 5 and Table 4 data"
for f in fig1_study_design.svg fig1_study_design.png fig1_study_design.pdf fig2_forest_pooled.png fig2_forest_pooled.pdf fig4_clade_ratio.png fig4_clade_ratio.pdf \
         fig5_dairy_effects.png fig5_dairy_effects.pdf af4_seven_claims.png af4_seven_claims.pdf table4_seven_claims.tsv; do cmpa results/figures_draft18/$f "$G"; done
for f in fig3_calibration.png fig3_calibration.pdf; do cmpa results/figures_manuscript/$f "$G"; done
G="Additional files 1 to 4"
for f in Additional_file_1_targets.xlsx Additional_file_2_estimates.xlsx Additional_file_3_analysis_plan.docx Additional_file_4_seven_claims.pdf; do cmpa results/additional_files_draft18/$f "$G"; done

# ---- 4. end to end from the processed inputs (count tables and harmonized phenotypes)
dump() { local t=$1 a="" c; case $t in
  nelore_*|charolais_mrna_*) c="scripts/21_estimate_targets.R ${t%_*} ${t##*_}";;
  alberta_*_all) local r=${t#alberta_}; c="scripts/65_estimate_alberta.R ${r%_all}"; a="ALBERTA_ALL=1";;
  wallace_prokaryotes_depth) c="scripts/81_estimate_methane.R wallace_prokaryotes real"; a="DEPTH=1";;
  *) c="scripts/81_estimate_methane.R $t real";; esac
  echo "env $a A46=1 RUMEN_ROOT='$R' '$RR' $c"; }
T="nelore_species nelore_genus nelore_family nelore_phylum charolais_mrna_species charolais_mrna_genus charolais_mrna_family charolais_mrna_phylum alberta_genus_all alberta_family_all difford_bact difford_arch wallace_prokaryotes wallace_prokaryotes_depth"
D=""; for t in $T; do D="$D $(dump $t) &&"; done
step dumps_a46 "${D} env A48=1 TRAIT=ch4_dmi RUMEN_ROOT='$R' '$RR' scripts/81_estimate_methane.R wallace_prokaryotes real && ls work/a46_inputs/*.rds | wc -l"
aside results/late_a46/effect_size
step clr_estimates_112 "RUMEN_ROOT='$R' '$RR' scripts/112_a46_effect_size.R"
step seven_claims_128 "A48_OUT='$WORK/out/late_a48' RUMEN_ROOT='$R' '$RR' scripts/128_a48_seven.R"
G="End to end from the processed inputs: per-cohort CLR estimates (112) and the 48.1 seven claims with Table 4 patterns (128)"
for f in cohort_estimates pooled exclusion_counts a37_perm_p_recomputed; do cmpa results/late_a46/effect_size/$f.tsv "$G"; done
for f in seven_claims.tsv seven_pooled.tsv concordance.tsv secondary_matched.tsv counts_by_classification.tsv s1_yield.tsv fig_seven_claims.png fig_seven_claims.pdf; do
  cmpf "$WORK/out/late_a48/seven/$f" "$R/results/late_a48/seven/$f" "results/late_a48/seven/$f" "$G"; done

# ---- 5. permutation nulls behind every rescaled W
step permutation_nulls "'$PY' scripts/147_compare.py perm '$R' '$WORK/out/perm.tsv'"

# ---- 6. record and report
"$PY" - "$WORK" "$URL" "$REC" "$COMMIT" "$ARCH" "$CM" "$RR" <<'EOF'
import sys, json, subprocess, platform, datetime, os
work, url, rec, commit, arch, cm, rr = sys.argv[1:]
sh = lambda c: subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip()
osn = sh(". /etc/os-release && echo $PRETTY_NAME")
mach = {"Machine": f"{platform.node()} ({osn}, Linux {platform.release()})", "CPU": f"{sh('lscpu | sed -n s/^Model\\ name:\\ *//p')} ({os.cpu_count()} threads)",
        "Memory": sh("free -g | awk '/^Mem:/{print $2\" GB\"}'"), "Threads per process": "1 (OMP, OpenBLAS and MKL set to 1, as in the original runs)"}
rv = sh(f"'{rr}' -e 'cat(R.version.string, \"|\", paste(sapply(c(\"ANCOMBC\", \"metafor\", \"ggplot2\", \"patchwork\", \"dplyr\", \"zCompositions\"), function(p) paste(p, packageVersion(p))), collapse=\"; \"))'")
pv = sh("python -c \"import sys,pandas,numpy,openpyxl,PIL; print('Python', sys.version.split()[0], '| pandas', pandas.__version__, '| numpy', numpy.__version__, '| openpyxl', openpyxl.__version__, '| Pillow', PIL.__version__)\"")
soft = {"Environment manager": f"{os.path.basename(cm)} {sh(cm + ' --version').splitlines()[0] if sh(cm + ' --version') else ''}",
        "R (env/rumen-r.yml)": rv, "Python (env/reporting.yml)": pv, "pandoc": sh("pandoc --version | head -1"), "rsvg-convert": sh("rsvg-convert --version"),
        "pdftotext": sh("pdftotext -v 2>&1 | head -1"), "git": sh("git --version")}
cmds = [f"systemd-run --user --unit=repro-v104 -p MemoryMax=32G -p MemorySwapMax=0 --setenv=PATH=... --setenv=SSH_AUTH_SOCK=... --working-directory=/data/rumen_multicohort bash scripts/147_reproduction_check.sh {work}", "# commands run by the script, in order (from the clone root):"] + open(f"{work}/out/commands.txt").read().splitlines()
json.dump({"date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M %Z").strip(), "url": url, "record": rec, "commit": commit, "archives": arch or "see log", "work": work,
           "machine": mach, "software": soft, "commands": cmds}, open(f"{work}/out/run.json", "w"), indent=1)
EOF
"$PY" "$R/scripts/147_compare.py" report "$WORK"
echo "$(date '+%F %T') done: $WORK/REPRODUCTION_CHECK.md"
