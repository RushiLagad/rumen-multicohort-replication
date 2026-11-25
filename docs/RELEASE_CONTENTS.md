# Release contents (v1.0.4-submission)

What the submitted version contains, where each part lives, and what is not deposited and why (draft 20, item 4). "Repository" means the git tag, which is also the repository archive on Zenodo (`rumen-multicohort-replication-v1.0.4-submission.tar.gz`). "Zenodo" means one of the four data archives in the same record; their paths are relative to the repository root. `docs/REPRODUCTION_CHECK.md` records a rebuild from a fresh clone and the four Zenodo archives in new environments.

## 1. Analysis scripts

| Item | Where |
|---|---|
| Every script, numbered in the order written (downloads, read processing, profiling, estimates, permutation nulls, pooling, late analyses, figures, tables, Additional files, number manifests, reproduction check) and the shared libraries (`lib_*.R`) | Repository, `scripts/`; `scripts/README.md` maps numbers to stages |
| The rules each script implements, with the date each was written relative to the results | Repository, `plan/ANALYSIS_PLAN.md`; `docs/DEVIATIONS.md`; `docs/ANALYSIS_CLASSES.md` |

## 2. Inputs

**Processed inputs (deposited).** Every model in the paper starts from these files.

| Input | Where |
|---|---|
| GTDB r207 genus and species count tables, Nelore metagenomes (52 animals) and Charolais metatranscriptomes (30 animals) | Zenodo, `profiles_gtdb_full_counts.tar.gz` (`profiles_gtdb_full/`) |
| SILVA 138.1 genus count tables, Alberta, Difford (bacterial and archaeal) and Wallace 16S, and the Charolais 16S layer; methanogen clade counts for M1 | Repository, `profiles_16s/` |
| Genus lineages used to build family and phylum tables | Repository, `results/replication/genus_lineage.tsv` (GTDB) and `results/16s/*/genus_lineage_silva.tsv` (SILVA) |
| Per-animal design and phenotype tables, run tables and the Alberta overlap list | Repository, `harmonized/` |
| Published targets and their mapping keys | Repository, `plan/target_keys.tsv` |
| Sample map: every run in the five deposits, with its BioSample, the animal ID used in the phenotype tables, the count table and column it enters, and the reason when it was not analysed | Repository, `results/sample_map/sample_map.tsv` (script 148); `summary.tsv` gives the identifier each count table uses |
| Read-processing record: fastp reports, host-removal logs, read accounting, Kraken 2 reports and Bracken outputs for the two shotgun cohorts | Zenodo, `read_qc_and_classification.tar.gz` |

**Raw inputs (public; not deposited).** Exact access:

| Input | Accession or source | How to obtain |
|---|---|---|
| Nelore shotgun metagenomes | PRJNA987743 | `scripts/01_fetch_nelore.sh`, or the faster AWS SRA mirror route `scripts/01d_fetch_via_s3.sh`; run lists in `PRJNA987743/` and `metadata/` |
| Charolais metatranscriptomes and 16S | PRJEB70895 | `scripts/09_prefetch_sra.sh charolais` and `scripts/09b_convert_sra.sh`; runs in `harmonized/runs.tsv` (`scripts/05b_charolais_mrna.sh`) |
| Alberta, Difford and Wallace 16S | PRJNA393057, PRJEB28065, PRJNA517480 | `scripts/50_download_ena.sh` (ENA, MD5-checked), with `scripts/51_ncbi_fallback.sh` for runs ENA fails to serve |
| Phenotype source files (published supplements and the INRAE Charolais sample sheet) | Conteville et al. 2023; INRAE doi:10.57745/ACNRVQ; Li et al. 2019; Difford et al. 2018; Wallace et al. 2019 | `docs/PHENOTYPES.md` gives each file, its MD5, the columns read and the join; `scripts/02b_fetch_phenotypes.sh` downloads the first two. The derived per-animal tables are in `harmonized/` |
| Host reference ARS-UCD2.0, Kraken 2 Standard index | Public downloads | `scripts/06a_fetch_refs.sh` |
| Kraken 2 GTDB r207 databases (capped and full resolution) | Struo2 GTDB release 207 library and taxdump | Built by `scripts/10_build_gtdb_db.sh`, `10b_renumber_taxids.py`, `10c_bracken_build.sh` and `13_twopass_db.sh`; details in `docs/DATABASES.md` |
| SILVA 138.1 reference for DADA2 | SILVA 138.1 nr | `scripts/31_silva_trainset.sh` |

**Sample identifiers in the count tables.** Nelore and Charolais: the SRA or ENA run accession with the suffix `_clean` (the two runs of Charolais animal 8574 are separate columns, summed into one animal profile before modelling). Alberta: the SRA run accession. Difford: the animal ID of the authors' phenotype table (Difford et al. 2018, S4), with all runs of that animal's BioSample summed. Wallace: the deposit's sample alias (BREED_FARM_Sn), with the two sequencing rounds summed. `results/sample_map/sample_map.tsv` links each of these to the deposit run and BioSample accessions.

## 3. Software environments

| Environment | Used for | File |
|---|---|---|
| rumen-r (R 4.3.3, ANCOM-BC2 2.4.0, metafor 4.8.0) | Every estimate, permutation null, pooling, late analysis and R figure | `env/rumen-r.yml`, `env/rumen-r_sessionInfo.txt` |
| reporting (Python 3.12.11, pandas 3.0.3, openpyxl 3.1.5, pandoc 3.7.0.2, rsvg-convert 2.58.0) | Number manifests, Additional files, Figure 1 render | `env/reporting.yml` |
| rumen-late (R 4.5.3, MicrobiomeStat) | LinDA (Amendment 44) | `env/rumen-late.yml`, `env/rumen-late_sessionInfo.txt` |
| silva_taxonomy (R 4.5.3, DADA2 1.38.0) and cutadapt | 16S read processing | `env/silva_taxonomy.yml`, `env/silva_taxonomy_sessionInfo.txt`, `env/cutadapt.yml` |
| Kraken 2 2.17.2, Bracken 3.0.1 | Shotgun profiling | Versions in `docs/DATABASES.md` |

## 4. Permutation outputs

| Reported quantity | Stored permutation output | Where |
|---|---|---|
| Rescaled ANCOM-BC2 W, permutation-null SD and empirical permutation p of every target estimate (Amendment 37; Tables 2 to 4; Additional file 2) | 200 within-design permutation fits per table (seeds 1001 to 1200; one estimate has 198 usable fits) | Repository, `results/rfi_null/` (feed efficiency) and `results/methane_dryrun/` (methane) |
| Untargeted genus scan, Figure 3, spike-in check and the Amendment 48 filtered null | 600 ANCOM-BC2 fits on permuted phenotypes | Zenodo, `scan_null.tar.gz` (`results/scan_null/`) |
| ANCOM-BC2 Freedman-Lane p (Amendment 46) | Target values per permutation in `results/late_a46/ancombc2_fl/*_targets.tsv` (repository); full ANCOM-BC2 outputs per permutation | Repository and Zenodo, `late_a46_ancombc2_fl.tar.gz` |
| CLR Freedman-Lane p (9,999 permutations; Amendment 46 and 48) | Not stored per permutation: regenerated exactly from the fixed seeds (461001 on) by `scripts/112_a46_effect_size.R` | Counts b and B in `results/late_a46/effect_size/cohort_estimates.tsv` |
| Null-scan check (Amendment 47) and other late permutation results | Summaries per replicate | Repository, `results/late_a47/` and `results/late_*` |

## 5. Source data for tables and figures

| Item | Where |
|---|---|
| Every number in the manuscript with the result file it comes from | Repository, `results/manifest/` (`headline_numbers.tsv`, `table_numbers.tsv`, `draft17_numbers.tsv` to `draft20_numbers.tsv`), `results/late_a47/manifest_a47.tsv`, `results/late_a48/manifest_a48.tsv` |
| Figures 1, 2, 4, 5, Additional file 4 and the Table 4 data, with script and source for each | Repository, `results/figures_draft18/` and its `README.md` |
| Figure 3 | Repository, `results/figures_manuscript/fig3_calibration.*` (script 103; reads `results/scan_null/` from Zenodo) |
| Additional files 1 to 4 and their sources | Repository, `results/additional_files_draft18/` (scripts 141, 142, 144, 146) |

## 6. Not deposited

| Item | Size | Why, and where it is |
|---|---|---|
| Raw sequence reads | more than the 540 GB of processed reads | Public in ENA and SRA under the accessions above |
| Quality-filtered, host-depleted reads (`processed/`) | about 540 GB | Rebuilt from the public reads by the read-processing scripts; per-sample reports and logs are on Zenodo |
| Kraken 2 databases | 78 GB (Standard); 54.6 GB hash table (GTDB r207 full resolution) | Downloaded or rebuilt from public sources (`docs/DATABASES.md`); build logs kept on the analysis host (`/data/refs/logs/`) |
| Phenotype source files (`raw/phenotypes/`, supplements in `*.xlsx`) | small | Other authors' published files; obtain as in `docs/PHENOTYPES.md`. The Wallace et al. 2019 values in `harmonized/` keep that article's CC BY-NC 4.0 licence |
| Model-input dumps (`work/a46_inputs/`) | small | Hold per-animal phenotype values; regenerated by running scripts 21, 65 and 81 with `A46=1` (as in `docs/REPRODUCTION_CHECK.md`) |
| Outputs for PRJNA448333 | | Cohort dropped before analysis (its phenotypes could not be traced; Amendments 20, 24, 26, 39) |
