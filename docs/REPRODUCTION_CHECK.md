# Reproduction check (clean environment)

Draft 20, item 3, requested by Dr. Goswami. Run 2026-11-25 15:55 by scripts/147_reproduction_check.sh; reporting only (no new analysis, and nothing in the working repository was changed).

**Result: PASS: every regenerated output matches the committed file.**

## What was done

1. Fresh clone of git@github.com:RushiLagad/rumen-multicohort-replication.git (branch master) at commit `6f9654e6` into a new folder (/data/repro_v104/repo).
2. The four data archives downloaded from the published Zenodo record 24327001 (https://zenodo.org/records/24327001), MD5 checked against the record (read_qc_and_classification.tar.gz, scan_null.tar.gz, profiles_gtdb_full_counts.tar.gz, late_a46_ancombc2_fl.tar.gz), and unpacked at the clone root.
3. Two new conda environments built only from the clone's environment files: `env/rumen-r.yml` (R analysis and figures) and `env/reporting.yml` (Python reporting, pandoc, rsvg-convert).
4. From the committed results, regenerated and compared with the committed files: the number manifests behind Tables 1 to 5 (scripts 126, 127, 143, 145, 146), Figures 1 to 5 (140 and rsvg-convert, 139, 103) and Additional files 1 to 4 (142).
5. End to end from the processed inputs (genus and species count tables, harmonized per-animal phenotypes): the Amendment 46 model-input dumps (21, 65, 81), the per-cohort CLR partial correlations with 9,999 Freedman-Lane permutations (112), and the 48.1 seven-claim table with its Table 4 patterns (128), compared with the archived outputs.
6. The permutation-null SD, number of permutations and b/B behind every rescaled W in the manuscript, recomputed from the stored null fits.

## Machine and software

| Item | Value |
|---|---|
| Machine | linux-mini-gini (Ubuntu 26.04.1 LTS, Linux 6.17.0-1017-oem) |
| CPU | AMD Ryzen Threadripper PRO 7955WX 16-Cores (32 threads) |
| Memory | 122 GB |
| Threads per process | 1 (OMP, OpenBLAS and MKL set to 1, as in the original runs) |
| Environment manager | mamba 2.1.1 |
| R (env/rumen-r.yml) | R version 4.3.3 (2024-02-29); ANCOMBC 2.4.0; metafor 4.8.0; ggplot2 3.5.2; patchwork 1.3.2; dplyr 1.1.4; zCompositions 1.5.0.5 |
| Python (env/reporting.yml) | Python 3.12.11; pandas 3.0.3; numpy 2.4.6; openpyxl 3.1.5; Pillow 13.0.0 |
| pandoc | pandoc 3.7.0.2 |
| rsvg-convert | rsvg-convert version 2.58.0 |
| pdftotext | pdftotext version 25.07.0 |
| git | git version 2.53.0 |

## Commands

```
systemd-run --user --unit=repro-v104 -p MemoryMax=32G -p MemorySwapMax=0 --setenv=PATH=... --setenv=SSH_AUTH_SOCK=... --working-directory=/data/rumen_multicohort bash scripts/147_reproduction_check.sh /data/repro_v104
# commands run by the script, in order (from the clone root):
git clone -q --branch master 'git@github.com:RushiLagad/rumen-multicohort-replication.git' '/data/repro_v104/repo' && cd '/data/repro_v104/repo' && { [ -z '' ] || git checkout -q ''; }
'/home/rlagad/mambaforge/bin/mamba' env create -y -q -p '/data/repro_v104/envs/rumen-r' -f '/data/repro_v104/repo/env/rumen-r.yml'
'/home/rlagad/mambaforge/bin/mamba' env create -y -q -p '/data/repro_v104/envs/reporting' -f '/data/repro_v104/repo/env/reporting.yml'
'/data/repro_v104/envs/reporting/bin/python' - '24327001' '/data/repro_v104/dl' <<'EOF'
import sys, json, hashlib, urllib.request
rec, d = sys.argv[1], sys.argv[2]
r = json.load(urllib.request.urlopen(urllib.request.Request(f'https://zenodo.org/api/records/{rec}', headers={'Accept': 'application/json'})))
fs = [f for f in r['files'] if f['key'].endswith('.tar.gz') and not f['key'].startswith('rumen-multicohort-replication-')]
assert len(fs) == 4, [f['key'] for f in fs]
for f in fs:
    p = f'{d}/{f["key"]}'; urllib.request.urlretrieve(f['links']['self'], p); h = hashlib.md5(open(p, 'rb').read()).hexdigest()
    ok = 'md5:' + h == f['checksum']; print(f['key'], h, 'matches the record' if ok else 'DOES NOT MATCH ' + f['checksum']); assert ok
EOF
for a in '/data/repro_v104'/dl/*.tar.gz; do tar -xzf "$a" -C '/data/repro_v104/repo' || exit 1; done; ls -d profiles_gtdb_full/nelore results/scan_null results/late_a46/ancombc2_fl
RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/reporting/bin/python' scripts/126_manifest_tables.py
RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/reporting/bin/python' scripts/127_manifest_draft17.py
'/data/repro_v104/envs/reporting/bin/python' scripts/143_manifest_draft18.py '/data/repro_v104/repo'
'/data/repro_v104/envs/reporting/bin/python' scripts/145_manifest_draft19.py '/data/repro_v104/repo'
'/data/repro_v104/envs/reporting/bin/python' scripts/146_p_compare.py '/data/repro_v104/repo'
'/data/repro_v104/envs/reporting/bin/python' scripts/148_sample_map.py '/data/repro_v104/repo'
'/data/repro_v104/envs/reporting/bin/python' scripts/140_figure1_compact.py '/data/repro_v104/repo' && cd results/figures_draft18 && rsvg-convert -d 300 -p 300 -b white -f png -o fig1_study_design.png fig1_study_design.svg && rsvg-convert -b white -f pdf -o fig1_study_design.pdf fig1_study_design.svg
RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/139_figures_draft18.R
RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/103_figures_2to4.R
'/data/repro_v104/envs/reporting/bin/python' scripts/142_additional_files_draft18.py '/data/repro_v104/repo'
 env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R nelore species && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R nelore genus && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R nelore family && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R nelore phylum && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R charolais_mrna species && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R charolais_mrna genus && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R charolais_mrna family && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/21_estimate_targets.R charolais_mrna phylum && env ALBERTA_ALL=1 A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/65_estimate_alberta.R genus && env ALBERTA_ALL=1 A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/65_estimate_alberta.R family && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/81_estimate_methane.R difford_bact real && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/81_estimate_methane.R difford_arch real && env  A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/81_estimate_methane.R wallace_prokaryotes real && env DEPTH=1 A46=1 RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/81_estimate_methane.R wallace_prokaryotes real && env A48=1 TRAIT=ch4_dmi RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/81_estimate_methane.R wallace_prokaryotes real && ls work/a46_inputs/*.rds | wc -l
RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/112_a46_effect_size.R
A48_OUT='/data/repro_v104/out/late_a48' RUMEN_ROOT='/data/repro_v104/repo' '/data/repro_v104/envs/rumen-r/bin/Rscript' scripts/128_a48_seven.R
'/data/repro_v104/envs/reporting/bin/python' scripts/147_compare.py perm '/data/repro_v104/repo' '/data/repro_v104/out/perm.tsv'
```

## Steps

| Step | Exit code | Seconds |
|---|---|---|
| clone | 0 | 7 |
| env_rumen_r | 0 | 21 |
| env_reporting | 0 | 6 |
| zenodo_download | 0 | 328 |
| unpack | 0 | 2 |
| manifest_126 | 0 | 6 |
| manifest_127 | 0 | 2 |
| manifest_143 | 0 | 0 |
| manifest_145 | 0 | 0 |
| p_compare_146 | 0 | 0 |
| sample_map_148 | 0 | 11 |
| figure1_140 | 0 | 0 |
| figures_139 | 0 | 3 |
| figures_103 | 0 | 24 |
| additional_files_142 | 0 | 16 |
| dumps_a46 | 0 | 65 |
| clr_estimates_112 | 0 | 14 |
| seven_claims_128 | 0 | 1 |
| permutation_nulls | 0 | 2 |

## Comparisons

Status: identical = byte-identical; identical content = same pixels (PNG), same extracted text (PDF, DOCX; these files embed creation times), same cell values (XLSX) or same table values (TSV); numerically identical = numbers agree to 1e-9. The `commit` or `repo_commit` column of the manifests records the clone's HEAD and is not compared.

### Tables 1 to 5: number manifests (each script stops if a printed value differs from its result file)

| File | Status | Detail |
|---|---|---|
| `results/manifest/table_numbers.tsv` | identical content | all cells equal (repo_commit column ignored) |
| `results/manifest/draft17_numbers.tsv` | identical content | all cells equal (repo_commit column ignored) |
| `results/manifest/draft18_numbers.tsv` | identical content | all cells equal (repo_commit column ignored) |
| `results/manifest/draft19_numbers.tsv` | identical content | all cells equal (repo_commit column ignored) |
| `results/manifest/draft20_numbers.tsv` | identical content | all cells equal (commit column ignored) |
| `results/manifest/draft21_numbers.tsv` | identical content | all cells equal (commit column ignored) |
| `results/late_p_compare/seven_claims_p.tsv` | identical | byte-identical |
| `results/late_p_compare/all_targets_p.tsv` | identical | byte-identical |
| `results/late_p_compare/summary.tsv` | identical | byte-identical |

### Sample map: deposit run to BioSample, animal and count-table column (148)

| File | Status | Detail |
|---|---|---|
| `results/sample_map/sample_map.tsv` | identical | byte-identical |
| `results/sample_map/summary.tsv` | identical | byte-identical |

### Figures 1 to 5 and Table 4 data

| File | Status | Detail |
|---|---|---|
| `results/figures_draft18/fig1_study_design.svg` | identical | byte-identical |
| `results/figures_draft18/fig1_study_design.png` | identical | byte-identical |
| `results/figures_draft18/fig1_study_design.pdf` | identical content | extracted text equal (450 words) |
| `results/figures_draft18/fig2_forest_pooled.png` | identical | byte-identical |
| `results/figures_draft18/fig2_forest_pooled.pdf` | identical content | extracted text equal (292 words) |
| `results/figures_draft18/fig4_clade_ratio.png` | identical | byte-identical |
| `results/figures_draft18/fig4_clade_ratio.pdf` | identical content | extracted text equal (69 words) |
| `results/figures_draft18/fig5_dairy_effects.png` | identical | byte-identical |
| `results/figures_draft18/fig5_dairy_effects.pdf` | identical content | extracted text equal (95 words) |
| `results/figures_draft18/af4_seven_claims.png` | identical | byte-identical |
| `results/figures_draft18/af4_seven_claims.pdf` | identical content | extracted text equal (169 words) |
| `results/figures_draft18/table4_seven_claims.tsv` | identical | byte-identical |
| `results/figures_manuscript/fig3_calibration.png` | identical content | all 2007 x 1080 pixels equal |
| `results/figures_manuscript/fig3_calibration.pdf` | identical content | extracted text equal (124 words) |

### Additional files 1 to 4

| File | Status | Detail |
|---|---|---|
| `results/additional_files_draft18/Additional_file_1_targets.xlsx` | identical content | 7 sheets, 6646 cells equal |
| `results/additional_files_draft18/Additional_file_2_estimates.xlsx` | identical content | 58 sheets, 545954 cells equal |
| `results/additional_files_draft18/Additional_file_3_analysis_plan.docx` | identical content | extracted text equal (43067 words) |
| `results/additional_files_draft18/Additional_file_4_seven_claims.pdf` | identical content | extracted text equal (169 words) |

### End to end from the processed inputs: per-cohort CLR estimates (112) and the 48.1 seven claims with Table 4 patterns (128)

| File | Status | Detail |
|---|---|---|
| `results/late_a46/effect_size/cohort_estimates.tsv` | identical | byte-identical |
| `results/late_a46/effect_size/pooled.tsv` | identical | byte-identical |
| `results/late_a46/effect_size/exclusion_counts.tsv` | identical | byte-identical |
| `results/late_a46/effect_size/a37_perm_p_recomputed.tsv` | identical | byte-identical |
| `results/late_a48/seven/seven_claims.tsv` | identical | byte-identical |
| `results/late_a48/seven/seven_pooled.tsv` | identical | byte-identical |
| `results/late_a48/seven/concordance.tsv` | identical | byte-identical |
| `results/late_a48/seven/secondary_matched.tsv` | identical | byte-identical |
| `results/late_a48/seven/counts_by_classification.tsv` | identical | byte-identical |
| `results/late_a48/seven/s1_yield.tsv` | identical | byte-identical |
| `results/late_a48/seven/fig_seven_claims.png` | identical | byte-identical |
| `results/late_a48/seven/fig_seven_claims.pdf` | identical content | extracted text equal (150 words) |

### Permutation nulls behind the rescaled W

62 of 62 target-cohort estimates reproduced (SD, B, b/B and rescaled W to 1e-9) from 2600 stored null fits in `results/rfi_null/` and `results/methane_dryrun/`.

## Not covered

Read processing and taxonomic profiling (quality control, host removal, Kraken2 and Bracken against the GTDB r207 database, DADA2 against SILVA 138.1) were not rerun: they start from about 540 GB of reads and the Kraken2 database, which are not deposited. Their per-sample reports and logs are in the deposit (`read_qc_and_classification.tar.gz`), and the count tables they produced are the processed inputs used in step 5. ANCOM-BC2 model fits and their 200-permutation nulls were not refitted; step 6 recomputes every rescaled W and permutation p from the stored null fits instead.

Environments are built new from the environment files; if the environment manager links packages from a local package cache, the check does not show that every pinned package can still be downloaded.

