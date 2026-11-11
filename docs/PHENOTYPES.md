# Phenotype sources

The per-animal phenotype tables are not stored in this repository. `raw/phenotypes/` is excluded in `.gitignore`, and so are spreadsheets (`*.xlsx`, `*.XLSX`), because they are other authors' published supplements. This file says where each cohort's table comes from, which part of it is read, and which script reads it, so the inputs can be fetched again and checked against the md5 sums below.

The per-animal tables derived from these sources are committed in `harmonized/` and are what the analysis scripts read.

| Cohort | Accession | Trait used | Source file | Read by | Derived table (committed) |
|---|---|---|---|---|---|
| Nelore | PRJNA987743 | RFI and residual methane, continuous | Conteville et al. 2023, Supplementary Table 2 | `scripts/03b_harmonize_metadata.sh` | `harmonized/animal_metadata.tsv` |
| Charolais | PRJEB70895 | RFI class, high or low | INRAE dataset `mapping.tab` | `scripts/03b_harmonize_metadata.sh` | `harmonized/animal_metadata.tsv` |
| Alberta | PRJNA393057 | RFI, continuous | Li et al. 2019 (Microbiome 7:92), Additional file 1 | `scripts/64_harmonize_alberta.py` | `harmonized/alberta_animals.tsv` |
| Difford | PRJEB28065 | CH4, g per day | Difford et al. 2018, S4 Table | `scripts/82_harmonize_methane.py` | `harmonized/difford_animals.tsv` |
| Wallace | PRJNA517480 | CH4, g per day (primary) and g per kg DMI; RFI NRC (sensitivity) | Wallace et al. 2019, Data S1 | `scripts/82_harmonize_methane.py` | `harmonized/wallace_animals.tsv` |

## Nelore, PRJNA987743

- **Paper:** Conteville LC, da Silva JV, Andrade BGN, et al. Rumen and fecal microbiomes are related to diet and production traits in *Bos indicus* beef cattle. Front Microbiol 2023;14:1282851. doi:10.3389/fmicb.2023.1282851 (PMC10754987).
- **File:** Supplementary Table 2, saved as `PRJNA987743/metadata/Table_2.XLSX` (md5 `d14b103bcb39eee7cd389578f43b8343`). The replication targets come from Supplementary Table 7 of the same paper, `PRJNA987743/metadata/Table_7.XLSX` (md5 `23ffa41a483e03ee4e2262a8008ffe75`).
- **How to get it:** `scripts/02b_fetch_phenotypes.sh` downloads the supplement zip from Europe PMC (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10754987/supplementaryFiles`) and extracts it. PMC's direct links return a JavaScript check that curl cannot pass, so if Europe PMC fails, download Table 2 by hand from the article page.
- **What is read:** sheet `Sheet1`, header on row 2; rows with `Site` = rumen (52 animals). Columns `SampleID`, `Diet`, `Slaughter_day`, `Weight`, `RFI_GC`, `RFI_by_dieta`, `RME_GC`, `RME_by_dieta` and `#Raw_reads`.
- **Join:** `SampleID` equals the ENA `sample_alias` for all 52 animals. The script checks the join by read count: S2 raw reads / 2 / ENA read count must lie in [0.95, 1.10] for every animal (observed 1.0158 to 1.0436).

## Charolais, PRJEB70895

- **File:** `mapping.tab` from the INRAE dataset doi:10.57745/ACNRVQ (Recherche Data Gouv), committed as `PRJEB70895/metadata/mapping.tab` (md5 `3e87555a73920feb3648f19a25f706c7`). It is a sample sheet, not a spreadsheet, so it is in the repository.
- **How to get it:** `scripts/02b_fetch_phenotypes.sh`, from `https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/ACNRVQ`.
- **What is read:** columns `id`, `N Animal`, `Diet`, `Year`, `Sampling` and `RFI_1`. Only the extreme classes are used (`RFI_1` POSITIF = high RFI, NEGATIF = low RFI). No continuous RFI exists anywhere in the deposit.
- **Exclusions (plan section 2.2):** sample `20BAR-01326` (animal 8513, contradictory fields) and animal 8508 (no 16S). The 33 mid-RFI animals have no reads deposited. 30 animals remain.

## Alberta, PRJNA393057

- **Paper:** Li F, Li C, Chen Y, et al. Host genetics influence the rumen microbiota and heritable rumen microbial features associate with feed efficiency in cattle. Microbiome 2019;7:92. doi:10.1186/s40168-019-0699-1. The cohort is 709 cattle born in 2013 at the Kinsella ranch.
- **File:** Additional file 1 (Table S1) of the article, saved as `raw/phenotypes/Li2019_709_MOESM1.xlsx` (md5 `1f04d2e6a73b45d534bcfcec01d3b5aa`).
- **What is read:** sheet `Additional file 1`, header on row 3. Columns `Animal ID`, `Breed`, `Sex`, `Diet`, `Age (day)`, `Dry matter intake (DMI; kg/day)`, `Average daily gain (ADG; kg/day)` and `Residual feed intake (RFI; kg/day)`.
- **Join:** `Animal ID` equals the BioSample `source_material_id` of each PRJNA393057 run (`metadata/candidates/PRJNA393057.samples.tsv`). The primary analysis uses all 567 animals with a phenotype and a usable genus profile (Amendment 39).
- **Not used:** the PRJNA448333 steers (Li et al. 2019, Microbiome 7:6, the source of the Li target set). No per-animal phenotype for them is traceable in any public record (plan section 2.3, Amendments 20, 24 and 39).

## Difford, PRJEB28065

- **Paper:** Difford GF, Plichta DR, Løvendahl P, et al. Host genetics and the rumen microbiome jointly associate with methane emissions in dairy cows. PLoS Genet 2018;14(10):e1007580. doi:10.1371/journal.pgen.1007580 (PMC6200390).
- **File:** S4 Table of the article, saved as `raw/phenotypes/Difford2018_S4.bin` (md5 `a80247309c95e2e2013a740232e8b40f`). Despite the extension it is tab-separated text with CRLF line endings.
- **What is read:** all columns: `sample_id`, `Sequencing_Instrument`, `Sequencing_Batch`, `herd`, `Lactatation` (misspelled in the source), `days_in_milk` and `methane_grams_per_day`. 750 cows, 6 herds. Parity is coded 1, 2 or 3+.
- **Join:** ENA sample alias, with the submission prefix (`IC-#`, `NG-#`) and library tag (`_lib#`, `_Bact`, `_Arch`, `_V#`) removed, equals `sample_id`. The suffixes `_c` and `_J` are part of the cow ID (Amendment 28; `scripts/71_difford_map.py`).
- **Also kept:** S1 Table, `raw/phenotypes/Difford2018_S1.bin` (an .xlsx file, md5 `1e4eeda17bb1c449f550508d01824971`). It holds the authors' taxon estimates, was read only to check that the methane targets are independent of this cohort (Amendment 32), and supplies no phenotype. Amendment 19 also names a Dryad deposit for this cohort, but every Difford phenotype in the analysis comes from S4 Table through `scripts/82_harmonize_methane.py`.

## Wallace, PRJNA517480

- **Paper:** Wallace RJ, Sasson G, Garnsworthy PC, et al. A heritable subset of the core rumen microbiome dictates dairy cow productivity and emissions. Sci Adv 2019;5(7):eaav8391. doi:10.1126/sciadv.aav8391 (PMC6609165).
- **File:** Data S1, `aav8391_Data_S1.xlsx`, saved as `raw/phenotypes/Wallace2019_Data_S1.xlsx` (md5 `1d94d40e85e6b02b5343aac241caa5b6`). 1,016 cows on 7 farms in Finland, Italy, Sweden and the UK.
- **What is read:** sheet `Phenotype data`, data from row 4. Zero-based columns: 1 Cow ID, 4 farm code, 6 breed, 7 days in milk, 8 lactation number, 52 RFI NRC, 67 measurement method, 68 `Ruminal CH4 gas production, g/d`, 69 CH4 per kg DMI. The `Data key` sheet defines the columns.
- **Join:** cow key `<breed>_<farm>_S<Cow ID>`, where breed is NR (Nordic Red) or HOL. It matches all 1,012 sequenced cows (Amendment 34).
- **Caveats recorded in the plan:** CH4 is z-scored within farm because the measurement method differs by country. Intake in the UK and Italy was partly estimated with alkane markers, so CH4 per kg DMI and RFI NRC (DMI minus NRC 2001 predicted intake, not a regression RFI) are sensitivity traits only (Amendments 34 and 43).

## Re-creating `raw/phenotypes/`

1. Download each file above to the path given and check its md5 with `md5sum`.
2. Re-run the harmonising scripts: `bash scripts/03b_harmonize_metadata.sh`, `python3 scripts/64_harmonize_alberta.py` and `python3 scripts/82_harmonize_methane.py`. `82` needs `harmonized/wallace_runs.tsv` (from `scripts/78_wallace_map.py`), and `64` needs `profiles_16s/alberta/alberta_genus_counts.tsv`; both are committed.
3. The rebuilt tables should match the committed `harmonized/` files. At the time of writing their md5 sums were: `animal_metadata.tsv` `161c8bbc5cd31a10f6183b82fe105935`, `alberta_animals.tsv` `ca679b0cddc1374ce9de6fd40787dec6`, `difford_animals.tsv` `e1fb0ccef691a5d459a5a43552005850` and `wallace_animals.tsv` `ec518a40e99fc4dcbc747e2dd501e4a7`.
