# Taxonomic profiling software and databases

Versions as reported on the build and analysis host (linux-mini-gini) on 11 November 2026, together with the build details recorded in the analysis plan and the build logs in `/data/refs/logs/` (not in this repository).

## Software

| Tool | Version | Source of the version |
|---|---|---|
| Kraken 2 | 2.17.2 | `kraken2 --version`; the same version is recorded at database build time in Amendment 8 |
| Bracken | 3.0.1 | `bracken -v` |

## Databases

| Database | Role | Content | Build | Bracken k-mer distributions |
|---|---|---|---|---|
| Kraken 2 Standard, `k2_standard_20240605` | Sensitivity profile set (`profiles/`) | Prebuilt index from `genome-idx.s3.amazonaws.com/kraken/` (49,138 taxa, 78 GB; Amendment 1) | Downloaded by `scripts/06a_fetch_refs.sh` | Shipped with the index: 50, 75, 100, 150, 200, 250 and 300 bp |
| GTDB release 207, capped (`k2_gtdb207`) | Sensitivity profile set (`profiles_gtdb/`) | Struo2 GTDB r207 library, 62,443 genomes, 7,745,746 contigs after the sequence-to-taxid map was regenerated and taxids renumbered (Amendment 8) | `scripts/10_build_gtdb_db.sh`; estimated hash table 308 GB, built capped at 85 GB (`--max-db-size 85000000000`) | k = 35, 150 bp (`scripts/10c_bracken_build.sh`, 30 October 2026; 71,195 rows) |
| GTDB release 207, full resolution (`k2_gtdb207_full`) | **Primary profile set** (`profiles_gtdb_full/`) | The 12,111 GTDB r207 genomes with at least 100 clade reads in any sample on the capped database, 1,522,143 contigs (Amendments 9 and 11) | `scripts/13_twopass_db.sh`; estimated hash table 54.6 GB, built uncapped under a 96 GB guard | k = 35, at 100, 136 and 150 bp (`scripts/13_twopass_db.sh`), plus any further read length built per cohort by `scripts/08d_profile_gtdb_full.sh` |

GTDB taxonomy: the Struo2 release 207 taxdump (`http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release207/`), with taxids renumbered to sequential integers because Kraken 2 dropped sequences whose taxid exceeded 2^31 (Amendment 8; `scripts/10b_renumber_taxids.py`).

The 16S cohorts are not profiled with Kraken 2. They use DADA2 with SILVA 138.1 (see `env/` and Amendments 17, 21 and 35).

## Suggested Methods wording

Shotgun and metatranscriptome reads were classified with Kraken 2 (v2.17.2) and abundances re-estimated with Bracken (v3.0.1; k = 35) against a full-resolution Kraken 2 database built from the 12,111 GTDB release 207 genomes detected in a first pass over all samples. The Kraken 2 Standard database (build 2024-06-05) and a size-capped database of all 62,443 GTDB r207 genomes were used as sensitivity profile sets.
