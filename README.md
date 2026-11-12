# Do published rumen microbiome links replicate in new herds?

A multi-cohort test, against a dated analysis plan, of whether rumen microbes reported to be associated with feed efficiency or methane emissions in cattle show the same association in independent animals.

Rushikesh R. Lagad, Goswami Lab, Department of Animal Science, University of Arkansas (principal investigator: Dr. Aranyak Goswami). Contact: rlagad@uark.edu

**Status:** prepared for submission to *Microbiome* (Lagad and Goswami), November 2026; not yet peer reviewed. The version for submission is tag `v1.0.1-submission`, archived on Zenodo with the files that are too large for this repository: https://doi.org/10.5281/zenodo.24170173 (all versions). The earlier tag `v1.0-submission` (Zenodo 10.5281/zenodo.24170174) lacks the Figure 1 generator and uses earlier figure file names; v1.0.1-submission supersedes it with no change to analyses or results

## The question

Many studies report rumen microbes linked to feed efficiency or to methane, usually from a single herd, and few of these links have been retested in other animals. This project takes the published claims as they were reported, fixes in advance how each one will be tested, and tests them in herds that played no part in the original discovery.

## Cohorts

| Cohort | Public accession | Animals analysed | Data | Trait |
|---|---|---|---|---|
| Nelore beef cattle | PRJNA987743 | 52 | shotgun metagenomes | feed efficiency (residual feed intake) |
| Charolais beef cattle | PRJEB70895 | 30 | metatranscriptomes and 16S | feed efficiency |
| Alberta beef cattle (Li et al. 2019, Microbiome 7:92) | PRJNA393057 | 567 | 16S | feed efficiency |
| Dairy cows (Difford et al. 2018, PLoS Genetics; Zhang et al. 2020) | PRJEB28065 | 749 | 16S, bacterial and archaeal | methane, g per day |
| Dairy cows, seven farms in Finland, Italy, Sweden and the UK (Wallace et al. 2019, Science Advances) | PRJNA517480 | 1,009 | 16S | methane, g per day; feed efficiency as a secondary analysis |

All sequence data and phenotype tables come from public archives and published supplements. Sequence files are not stored in this repository; the accessions above are where to download them.

## How the tests work

1. **Targets fixed in advance.** Each published association (a taxon, its direction and its trait) is written down before any of the new cohorts are analysed.
2. **Same processing for every cohort.** Reads are processed the same way within each data type, and each cohort's design (diet, farm, herd, sequencing batch) is included in its model.
3. **A check that each cohort can detect anything at all.** Before a cohort's results are read, the same pipeline must find a strong, expected difference in those animals (for example between diets or farms). All five cohorts pass.
4. **Calibration.** The phenotype is shuffled within each design group 200 times and the analysis repeated, so every result is judged against what the method produces on data with no real signal.
5. **Pooling across cohorts.** Effects are combined across the independent cohorts, corrected for the number of tests, and an association counts as replicated only if it holds in the pooled result and survives dropping each cohort in turn.
6. **Every change recorded.** The analysis plan is a dated document; any change to methods is added as a dated entry stating whether it came before or after the affected results were seen, including corrections of earlier mistakes.

## Findings (submitted version)

- **Feed efficiency:** none of the 27 published associations that could be tested replicated. For 17 of them the data rule out an effect as large as typically reported; for the other 10 the data were too imprecise to say.
- **Methane:** none of the 11 testable associations replicated, and for all 11 the data rule out even modest effects in the published direction.
- **Within one herd, not across herds:** clear associations do appear inside a single cohort (in the Wallace dairy cows, consistently across its seven farms), but they are absent from the other large dairy cohort. Effects of lactation stage and parity did not carry between the two dairy herds either, although both were processed with the same 16S pipeline and SILVA taxonomy.
- **Method:** one widely used differential abundance method (ANCOM-BC2) reported many false positives when the phenotype was shuffled, which is why every result here is calibrated. In a spike-in, calibration brought its false-positive rate in the large cohorts from about 25% to 5% while still recovering about nine in ten planted effects of a modest size. Three further methods (a CLR linear model, LinDA and ALDEx2) lead to the same conclusions.

Analyses run after the primary results had been seen are labelled as late analyses in the plan (Amendments 44 and 45) and change no replication label.

## Repository layout

| Folder | Contents |
|---|---|
| `plan/` | The analysis plan with all dated changes (`ANALYSIS_PLAN.md`) and the table of published targets (`target_keys.tsv`) |
| `scripts/` | Every step of the analysis, numbered in the order they were written; see `scripts/README.md` |
| `harmonized/` | Per-animal design and phenotype tables assembled from the public sources |
| `metadata/`, `PRJNA987743/`, `PRJEB70895/` | Run lists and download records |
| `profiles_16s/` | Genus count tables for the 16S cohorts |
| `results/` | All outputs: per-cohort estimates, shuffled-data runs, pooled results and summary tables |
| `logs/` | Run logs |
| `env/` | Software environments (conda exports and R session information) |
| `raw/` | Small reference files only (published supplementary tables, reference sequences) |

## Reproducing the analysis

The software environments are in `env/`: analysis in R 4.3.3 with ANCOM-BC2 2.4.0 and metafor 4.8.0, and 16S processing in R 4.5.3 with DADA2 1.38.0. Sequence data are downloaded from the accessions above with the download scripts in `scripts/`. Two inputs are not in this repository because of their size and are on Zenodo with the submitted version: the 600 permutation null fits behind the untargeted genus scan and the spike-in check (`results/scan_null/`, about 470 MB) and the GTDB genus and species count tables for the Nelore and Charolais metagenomes and metatranscriptomes (`profiles_gtdb_full/`); unpack both at the repository root. Read-level intermediates (`processed/`, about 540 GB) are not deposited and are rebuilt from the public reads. The order of the steps is described in `scripts/README.md`, and the analysis plan states for every step what was decided and when.

## Data sources and credit

This work re-uses data generated and shared by the authors of the studies listed in the cohort table and of the papers that reported the original associations, cited in full in the analysis plan. Please cite the original studies alongside any use of this work.

## Licence

Code (`scripts/`, `env/`) is released under the MIT licence (`LICENSE`). Results, tables, figures and documentation written for this project are released under CC BY 4.0 (`LICENSE-data`). Data and files from other sources (sequence archives, published supplements, phenotype tables, and the article figures and supplementary files kept as download records in `PRJNA987743/metadata/`) remain under the terms of their original sources; see "Data sources and credit".
