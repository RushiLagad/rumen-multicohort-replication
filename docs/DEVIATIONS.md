# Deviations from the analysis plan

Source for Part 1 of Additional file 3. Part 2 of that file is plan/ANALYSIS_PLAN.md, unedited.


# How to read this file

Part 1 lists every change made to the analysis plan after version 1.0 (28 October 2026), in the order the amendments were written. For each change it gives the date, what changed, the timing stated in the plan, and when the change was made relative to the results. Part 2 is the full plan, reproduced unedited from `plan/ANALYSIS_PLAN.md` in the project repository.

Timing classes:

- **Before any association result.** Written before any target association had been estimated in any cohort.
- **Before the affected analysis.** Written before the analysis it governs was run, but after results from other cohorts or analyses had been seen.
- **After the affected results.** Written after results of the analysis it changes, or that motivated it, had been seen.
- **Record only.** A finding or correction to the record that changes no analysis rule.

Results first became available as follows. Nelore and Charolais RFI estimates, 2 November. Alberta RFI estimates, 3 November. Methane estimates, 6 November. Calibrated RFI estimates, 8 November. The de novo scan, 9 November. Farm, depth, second-estimator and Wallace RFI sensitivities, 10 November. Late analyses (Amendments 44 and 45), from 11 November.

# Part 1. Deviation table

| Amendment | Date | Change | Timing stated in the plan | Relative to results | Effect on labels or primary analysis |
|---|---|---|---|---|---|
| 1 | 2026-10-29 | Kraken2 Standard (k2_standard_20240605) fixed as the profiling database; 60% median classified-read adequacy test set | Before any profiling | Before any association result | None |
| 2 | 2026-10-29 | minimap2 replaces bowtie2 for host depletion; host-reference comparison on 8 Nelore samples | Pipeline deviation | Before any association result | None |
| 3 | 2026-10-29 | ARS-UCD2.0 used as host reference for every cohort; PRJNA987743 found deposited already host-depleted | After the reference comparison | Before any association result | None |
| 4 | Not stated | Text not present in the plan. Amendment 5 states that Amendment 4 reported the classified-read fractions of the first three Nelore samples (17.0, 19.4 and 16.0%), which Amendment 5 supersedes | Not stated | Before any association result (it precedes Amendment 5, 30 October) | None; no analysis rule depends on it |
| 5 | 2026-10-30 | Adequacy result corrected after a reporting bug (median 16.9%, fails 60%); Standard profile kept as sensitivity; GTDB r207 route chosen | After the adequacy result | Before any association result | Changes the database for the primary profiles |
| 6 | 2026-10-30 | Contamination screen for the Nelore Table 7 features specified | Pre-specified | Before any association result | None |
| 7 | 2026-10-30 | Contamination screen result: negative | Result | Record only | None |
| 8 | 2026-10-30 | GTDB r207 build defects found and corrected (taxid map, taxids above 2^31, capped table) | During the build | Before any association result | None |
| 9 | 2026-11-01 | Two-pass full-resolution database if adequacy fails; GTDB genus match keys | Written with 3 of 52 Nelore samples profiled | Before any association result | Defines the primary profile set |
| 10 | 2026-11-01 | Target key table committed (76 targets, 71 mapped) | Before any association test | Before any association result | None |
| 11 | 2026-11-02 | Two-pass genome selection threshold changed from at least 10 to at least 100 clade reads, because the first selection did not fit the table | Before the second-pass build | Before any association result | Changes a pre-specified rule |
| 12 | 2026-11-02 | 29 Nelore Table 7 genera with no selected genome declared untestable | Before the second-pass profile and any association test | Before any association result | 29 targets untestable |
| 13 | 2026-11-02 | Data handling and positive-control pass rules | Before any positive control or association test | Before any association result | None |
| 14, 14b, 14c, 14d | 2026-11-02 | Effect scale, pooling, classification and reporting rules where Section 10 was silent | Before any replication test or target estimate | Before any association result | Defines the labels |
| 15 | 2026-11-02 | Level 1 reproducibility analysis added | Added after the replication results were seen; committed before it was run | After the affected results | None (additional analysis) |
| 16 | 2026-11-02 | Amendment 7 misassignment comparison operationalised | Committed before it was run | Before the affected analysis | None |
| 17 | 2026-11-02 | Charolais 16S layer specified; descriptive only, never pooled | Before any Charolais 16S read was processed | Before the affected analysis | None |
| 18 | 2026-11-02 | Leave-one-animal-out influence check on Charolais mRNA estimates | Exploratory; before any leave-one-out fit was run | After the affected results | None (exploratory) |
| 19 | 2026-11-03 | Study expanded: Alberta added to the RFI arm; methane added as a co-primary outcome (Difford, Wallace) with new methane targets | Before any read or phenotype of a new cohort was downloaded | Before the affected analysis; after the Nelore and Charolais RFI results | Changes the scope of the study |
| 20 | 2026-11-03 | Provenance of the PRJNA448333 phenotypes found to be untraceable | No change to any analysis rule | Record only | None |
| 21 | 2026-11-03 | Paired-end 16S processing rules | Before any Alberta read was processed | Before the affected analysis | None |
| 22 | 2026-11-03 | Sample exclusion rules for 16S cohorts | Before any Alberta phenotype was joined | Before the affected analysis | None |
| 23 | 2026-11-03 | Alberta covariates restricted because diet, sex and collection date are confounded | After the design cross-tabs; before any Alberta association | Before the affected analysis | Defines the Alberta model |
| 24 | 2026-11-03 | Alberta primary analysis excludes 35 animals whose IDs match PRJNA448333 sample codes | Before any Alberta association was viewed (all-animal fits had completed, outputs unopened) | Before the affected analysis | Defines the Alberta primary (later reversed by 39) |
| 25 | 2026-11-03 | RFI pooling over Nelore, Charolais and Alberta | Before any Alberta estimate was viewed | Before the affected analysis | Defines the RFI pooling |
| 26 | 2026-11-03 | PRJNA448333 steers confirmed to belong to the Alberta cohort | No change to any analysis rule | Record only | None |
| 27 | 2026-11-04 | Wording corrections to Amendments 24 and 26 | No change to any rule or result | Record only | None |
| 28 | 2026-11-04 | Processing of the Difford and Wallace amplicons | Before any Difford or Wallace association | Before the affected analysis | None |
| 29 | 2026-11-04 | Difford reads subsampled to 50,000 pairs; explicit archaeal filter | Before any Difford association | Before the affected analysis | None |
| 30 | 2026-11-04 | Forward-read sensitivity table for Difford bacteria | Before any Difford association | Before the affected analysis | None |
| 31 | 2026-11-05 | Computational changes after a memory failure | Before any Difford association | Before the affected analysis | None |
| 32 | 2026-11-05 | Independence rule (Wallace M2 is re-analysis), Wallace without SE1 as sensitivity, M1 measurement | Before any methane association | Before the affected analysis | M2 in Wallace leaves the primary pooling |
| 33 | 2026-11-05 | Wallace read qualities, Wallace archaeal filter, Difford M1 gate | Before any Wallace association | Before the affected analysis | None |
| 34 | 2026-11-05 | Wallace phenotypes and covariates; disclosure of values displayed during a structure check | Before any Wallace phenotype was joined to microbiome data | Before the affected analysis | Defines the Wallace traits |
| 35 | 2026-11-05 | Error corrected: 16S family-level targets reassigned with SILVA families instead of the GTDB table | Correction (Alberta C2 had been estimated; r 0.054 became -0.037) | After the affected results | C2 estimates changed; no label changed |
| 36 | 2026-11-05 | Methane target rows and estimation details | Before any methane association | Before the affected analysis | None |
| 37 | 2026-11-06 | ANCOM-BC2 estimates calibrated by 200 within-design permutations | Methane: before any real methane fit. RFI: after the uncalibrated RFI results | Methane: before the affected analysis. RFI: after the affected results | No label changed |
| 38 | 2026-11-08 | Positive controls for the 16S cohorts; exploratory de novo scan | Committed before any of it was run | After the affected results (added after the primary results) | None (exploratory) |
| 39 | 2026-11-09 | Amendment 24 exclusion found to rest on a spurious ID match; Alberta all-animal analysis (567) becomes primary | After both Alberta versions had been estimated | After the affected results | Primary changed; labels identical under both versions |
| 40 | 2026-11-10 | Wallace leave-one-farm-out check of the scan signal | Committed before it was run; motivated by the de novo scan result | After the affected results | None (exploratory) |
| 41 | 2026-11-10 | Depth-adjusted Wallace sensitivity | Committed before it was run; motivated by an observed depth-methane association | After the affected results | None |
| 42 | 2026-11-10 | Second estimator (CLR linear model) for all targets. Pre-specified in Amendment 14 item 1 ("Sensitivity: t from lm on CLR"), not run at the time, and run late here | Committed before it was run | After the affected results (pre-specified analysis run late, after the primary results) | None in RFI; one methane label differs (N38, below) |
| 42 addendum | 2026-11-10 | CLR methane estimates calibrated by 200 permutations | Before running; after an uncalibrated CLR label change (N38) was seen | After the affected results | N38 Eubacterium "opposite sign" under calibrated CLR only |
| 43 | 2026-11-10 | Wallace RFI NRC sensitivity implemented as planned in Amendments 19 and 34 | Before any Wallace RFI estimate | Before the affected analysis | None |
| M1 source traits | 2026-11-10 | Descriptive note: the M1 source papers used methane yield and daily output | Recorded after the sensitivity results were seen | Record only | None |
| 44 | 2026-11-11 | Late analyses: LinDA and ALDEx2 (Section 6), Hedges' g for Charolais (Section 9, Amendment 14 item 1) and the Li et al. diversity claims (Section 3.1), all pre-specified and not run at the time; written to results/late_* only; record corrections (items in this file) | Committed before any of it was run | After the affected results (late analysis, run after all primary results were seen) | None: no label differs under LinDA or ALDEx2 (0 replicated in either arm); Hedges' g descriptive; none of the four diversity claims supported |
| 45 | 2026-11-11 | Late checks: positive controls for transfer (Difford vs Wallace covariate effects; Charolais diet effect across mRNA/GTDB and 16S/SILVA in the same animals), spike-in of known effects with and without calibration, and target-mapping concordance | Committed before any of it was run | After the affected results (late analysis, run after all primary results were seen; not pre-specified) | None: no label changes. Transfer: the Charolais diet effect carries across mRNA/GTDB and 16S/SILVA (weakly), dairy covariate effects do not carry between herds; spike-in: calibration lowers the false-positive rate from about 25% to 5% in the large cohorts and recovers about 90% of planted effects at calibrated r near 0.17; target abundances agree across platforms where a strong effect exists (rho 0.57 to 0.91), median 0.24 |
| 46 | 2026-11-12 | Response to external statistical review: CLR partial correlation as the effect-size estimand with HKSJ intervals and a coverage check that gates exclusion claims; Freedman-Lane covariate-preserving permutation; genus scan at 99,999 (shotgun) and 9,999 (16S) permutations; permutation p as (b + 1)/(B + 1); target strata and reporting fixes | Committed before any of it was run | After the affected results (late analysis, run after all primary and late results were seen; not pre-specified) | None: no label changes (replication criterion and Amendment 37 labels unchanged). Coverage check passed (per-cohort 0.944 to 0.966; pooled HKSJ 0.944, 0.92 at tau = 0.10) after a sign bug in the first coverage run was fixed and rerun; CLR exclusion of r = 0.30 in the published direction 9 of 27 RFI and 7 of 11 methane targets (ANCOM-BC2-derived: 17 and 11); scan with 99,999 or 9,999 Freedman-Lane permutations: 0 RFI genus discoveries, Wallace 60 methane genera that do not carry to Difford; 3 primary permutation p of 0 become 0.005 |

## Pre-specified items not run as planned

| Plan section | Pre-specified | What was done | Amendment |
|---|---|---|---|
| 6, Models | ALDEx2 and LinDA as sensitivity methods beside ANCOM-BC2 | Not run at the time; run late, after all primary results were seen; no label changed | 44 |
| 14 item 1 | CLR linear model as a sensitivity estimator ("Sensitivity: t from lm on CLR") | Not run at the time; run late as the second estimator | 42 |
| 9, Pooling; 14 item 1 | Hedges' g for the Charolais binary contrast | Not run at the time; run late, descriptive only | 44 |
| 3.1, Primary target set | Li et al. diversity claims (bacterial Shannon evenness and Shannon, archaeal Simpson and inverse Simpson), "tested separately" | Not run at the time; run late, two-sided; none supported | 44 |
| 5, Read processing, item 4 | HUMAnN3 functional profiling of the shotgun layers | Never run; no functional results are reported, so the Section 1 secondary question (does moving from DNA-based taxonomy to expressed function change the answer) is not addressed | None |
| 2.3, Cohorts | PRJNA448333 (Li et al.) as a conditional baseline cohort | Dropped: its phenotypes could not be traced to a citable source (Amendments 20, 24, 26, 39) | 20 |

## Departures from Section 12 ("What will not be done")

| Amendment | Section 12 rule concerned | Departure |
|---|---|---|
| 19 | "No cohort added or removed after results are seen"; "No target added to any set after results are seen" | Alberta, Difford and Wallace, and the methane targets, were added after the Nelore and Charolais RFI results had been seen |
| 20 | "No cohort added or removed after results are seen" | PRJNA448333, conditional in Section 2.3, was dropped after the Nelore and Charolais RFI results because its phenotypes could not be traced |
| 37 | "No change of correction method after seeing results" | Permutation calibration was adopted for the RFI arm after its uncalibrated results (the methane arm was calibrated before any methane fit) |
| 38 | Fixed target list (Section 3, "will not be extended after results are seen") | An untargeted genus scan across whole tables, outside the fixed targets; exploratory and never merged with the replication labels |
| 40 | "No subgroup or interaction analysis not named in this document" | Leave-one-farm-out and per-farm fits in Wallace |

Plan version 1.0 is titled for feed efficiency in beef cattle. Methane entered through Amendment 19. Its Section 13 deviation log was left empty, and every deviation is recorded as a dated amendment instead.



