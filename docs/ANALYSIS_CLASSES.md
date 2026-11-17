# Analysis classification (Amendment 48.10)

Every analysis reported in the paper, classed by when its rule was written relative to the results it affects. Source for the classification table in Additional file 3. The timing of each amendment is as recorded in `docs/DEVIATIONS.md` (Part 1 of Additional file 3); this table adds the class and nothing else.

Classes (Amendment 48.10):

- **Confirmatory.** Plan version 1.0 (28 October 2026), or an amendment committed before the affected results were seen.
- **Planned sensitivity.** Named as a sensitivity analysis before the affected results were seen. "Run late" marks one that was named in advance but run only after all primary results had been seen (Amendments 42 and 44).
- **Post hoc.** Written after the affected results had been seen, whether or not it was committed before it was run.

Results first became available as follows (from `docs/DEVIATIONS.md`): Nelore and Charolais RFI estimates, 2 November 2026; Alberta RFI, 3 November; methane, 6 November; calibrated RFI, 8 November; de novo scan, 9 November; farm, depth, second-estimator and Wallace RFI sensitivities, 10 November; late analyses (Amendments 44 to 46) from 11 November; Amendment 47 on 13 November; Amendment 48 on 17 November.

| # | Analysis | Where reported | Amendment(s) | Committed | Relative to the affected results | Class |
|-|------------|------|------|----|-------|------|
| 1 | Read processing, profiling databases, host depletion, genus tables | Methods; Table 1 | Plan 1.0; 1, 2, 3, 5, 8, 9, 11, 21, 28 to 31, 33 | 29 Oct to 5 Nov | Before any association result in the cohort concerned | Confirmatory |
| 2 | Target sets, key table and mapping rules | Methods; Fig. 1b; Additional file 1 | Plan 1.0 (Section 3); 10, 12, 19 (methane targets), 36 | 1 to 5 Nov | Before any association result for the targets concerned | Confirmatory |
| 3 | Detection checks (positive controls in each primary table) | Results; Table 1 | 13 | 2 Nov | Before any positive control or association test | Confirmatory |
| 4 | Feed efficiency replication test, Nelore and Charolais: ANCOM-BC2, REML pooling, Benjamini-Hochberg, replication label (uncalibrated) | Results; Table 3 | Plan 1.0 (Sections 6 to 10); 14, 14b to 14d, 25 | 2 Nov | Before any association result | Confirmatory |
| 5 | Alberta added to the feed efficiency arm (all-animal analysis as finally used) | Results; Fig. 2a | 19, 21 to 25 (rules); 39 (all-animal analysis made primary) | 3 Nov (rules); 9 Nov (switch) | Rules before any Alberta read; switch after both Alberta versions had been estimated | Rules confirmatory; the all-animal primary is post hoc (39) |
| 6 | Alberta without the 35 ID-matched animals | Table 3 versions; Additional file 2 | 24 | 3 Nov | Before any Alberta association was viewed | Planned sensitivity (primary until Amendment 39) |
| 7 | Methane arm: Difford and Wallace, ANCOM-BC2 with 200-permutation calibration, pooling and labels | Results; Fig. 2b; Table 3 | 19, 28 to 34, 36, 37 (methane) | 3 to 6 Nov | Before any methane association | Confirmatory |
| 8 | Permutation calibration of the feed efficiency arm (calibrated W, r and labels) | Results; Table 2; Fig. 3a | 37 (RFI) | 6 Nov | After the uncalibrated RFI results | Post hoc |
| 9 | Wallace without farm SE1 | Fig. 4; Table 3 versions | 32 | 5 Nov | Before any methane association | Planned sensitivity |
| 10 | Methane per kg dry matter intake (Wallace) | Fig. 4; Table 3 versions | 34 item 2 | 5 Nov | Before any Wallace phenotype was joined | Planned sensitivity |
| 11 | Wallace as a fourth feed efficiency cohort (RFI NRC) | Table 3 versions; Additional file 2 | 19, 34, 43 | 3 to 10 Nov | Before any Wallace RFI estimate | Planned sensitivity |
| 12 | Charolais 16S layer (descriptive, never pooled) | Additional file 2 | 17 | 2 Nov | Before any Charolais 16S read | Planned sensitivity |
| 13 | Family-level 16S target mapping from SILVA lineages | Methods; Additional file 1 | 35 | 5 Nov | After Alberta C2 had been estimated (correction) | Post hoc |
| 14 | Level 1 reproducibility re-test in the Nelore discovery animals | Results | 15 | 2 Nov | After the replication results | Post hoc |
| 15 | Charolais leave-one-animal-out influence check (exploratory) | Results | 18 | 2 Nov | After the Charolais results | Post hoc |
| 16 | 16S positive controls; untargeted genus scan (200 permutations) | Results; Fig. 3 | 38 | 8 Nov | After the primary results | Post hoc |
| 17 | Wallace leave-one-farm-out and per-farm fits | Results; Fig. 5 (farm-robust genera marked); Additional file 2 | 40 | 10 Nov | After the de novo scan result | Post hoc |
| 18 | Depth-adjusted Wallace analysis | Fig. 4; Table 3 versions | 41 | 10 Nov | After an observed depth-methane association | Post hoc |
| 19 | Second estimator: CLR linear model | Table 3 versions; Fig. 3b | Named in 14 item 1; run under 42 | 2 Nov (named); 10 Nov (run) | Named before any result; run after the primary results | Planned sensitivity, run late |
| 20 | Calibration of the CLR methane estimates | Table 3 versions | 42 addendum | 10 Nov | After an uncalibrated CLR label change (N38) | Post hoc |
| 21 | LinDA and ALDEx2 | Results; Table 3 versions | Plan 1.0 Section 6; run under 44 | 28 Oct (named); 11 Nov (run) | Named before any result; run after all primary results | Planned sensitivity, run late |
| 22 | Hedges' g for the Charolais contrast; Li et al. diversity claims | Results; Additional file 2 | Plan 1.0 Sections 3.1 and 9; 14 item 1; run under 44 | 28 Oct and 2 Nov (named); 11 Nov (run) | Named before any result; run after all primary results | Planned sensitivity, run late |
| 23 | Transfer positive controls, planted effects (spike-in), target-mapping concordance | Results | 45 | 11 Nov | After all primary results | Post hoc |
| 24 | CLR partial r estimand, Freedman-Lane permutation, HKSJ intervals and coverage check, high-resolution scan, evidence ledger and overlap audit | Methods; Results; Additional files 1 and 2 | 46 | 12 Nov | After all primary and late results | Post hoc |
| 25 | Interval rule (modified Knapp-Hartung) chosen by coverage, CLR labels, null-scan calibration, four-cohort RFI under CLR | Methods; Fig. 2; Table 3 | 47 and addenda | 13 Nov | After all primary, late and Amendment 46 results | Post hoc |
| 26 | Seven-claim evidence table, patterns and secondary criteria | Table 4; Additional file 4 | 48.1 | 17 Nov | After all primary, late, Amendment 46 and 47 results | Post hoc |
| 27 | ANCOM-BC2 configuration audit and filtered null calibration; Charolais structural zeros | Results; Table 2 | 48.2 | 17 Nov | As 26 | Post hoc |
| 28 | Harmonized Wallace vs Difford comparison with a days-in-milk positive control | Results; Fig. 5 | 48.3 | 17 Nov | As 26 | Post hoc |
| 29 | Wallace farm random effects and prediction intervals | Results; Additional file 2 | 48.4 | 17 Nov | As 26 | Post hoc |
| 30 | CLR transformation grid (12 cells) | Results; Additional file 2 | 48.5 | 17 Nov | As 26 | Post hoc |
| 31 | Charolais leave-one-animal-out under CLR | Results; Additional file 2 | 48.6 | 17 Nov | As 26 | Post hoc |
| 32 | Minimum detectable r, published effects, attenuation (arithmetic) | Results; Discussion; Additional file 2 | 48.7 | 17 Nov | As 26 | Post hoc |
| 33 | Mapping-confidence classes and MC1-only re-tabulation | Results; Additional file 1 | 48.8 | 17 Nov | As 26 | Post hoc |
| 34 | Clade-ratio robustness grid and influence refits | Results; Additional file 2 | 48.9 | 17 Nov | As 26 | Post hoc |

Notes.

- The study expansion (Amendment 19) was written after the Nelore and Charolais RFI results and before any read or phenotype of a new cohort was downloaded. It is a departure from plan Section 12, listed in Additional file 3, and its analyses are classed as confirmatory for the new cohorts only (rows 5, 7).
- Two post hoc amendments define parts of the primary analysis: the permutation calibration of the feed efficiency arm (row 8, Amendment 37) and the all-animal Alberta analysis (row 5, Amendment 39). Both are reported beside the versions they replaced (Table 3). No target is labelled replicated in any version, confirmatory or post hoc.
- Rows 24 to 34 are labelled in the manuscript and in their output files as late analyses, run after all primary results were seen.
