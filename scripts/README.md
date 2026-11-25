# Scripts

Scripts are numbered in the order they were written, which follows the order of the analysis. Each one starts with a comment saying what it does and how to run it.

| Numbers | Stage |
|---|---|
| 00 to 14 | Cohort checks, downloads, metadata, shotgun and metatranscriptome profiling, target mapping |
| 20 to 25 | Detection checks, target estimates and pooling for the first cohorts, re-testing in the discovery animals |
| 30 to 35 | Charolais 16S processing, leave-one-animal-out checks |
| 40, 41, overnight, status | Batch runners and status helpers |
| 50 to 66 | Alberta 16S download, processing and estimates; three-cohort feed efficiency pooling |
| 70 to 87 | Methane cohorts: read processing, methanogen clade assignment, phenotypes, estimates, pooling |
| 88 to 100 | Permutation nulls (shuffled phenotypes) and rescaling of W by their SD, searches across all genera, farm-by-farm checks, sequencing depth check, second method, sensitivity analyses, summary tables |
| 101 to 104 | Manuscript figures (102 Figure 1, 103 Figures 2 to 4) and Additional files 1 and 2 (104, Python) |
| 105 to 107 | Late analyses, run after all primary results were seen (Amendment 44): LinDA and ALDEx2 fits and permutations (105), their pooling (106), Charolais Hedges' g and the Li et al. diversity indices (107) |
| 108 to 110 | Late checks (Amendment 45): positive controls for transfer (108), target-mapping concordance (109), spike-in summary (110); fits run through 105 |
| 111 | Figure 1 schematic, draft 17 (111_figure1_schematic.py; 111_render_svg.js renders the SVG) |
| 112 to 119 | Amendment 46: CLR partial correlation, Freedman-Lane permutation and pooling (112), interval coverage simulation (113), high-resolution genus scan (114), runner (115), Additional file 3 (116), farm report (117), evidence ledger (118), number manifest (119) |
| 120 to 125 | Amendment 47: pooled-interval rule chosen by coverage (120), CLR labels (121), null-scan check (122), runner (123), number manifest (124), image metadata removal (125) |
| 126, 127 | Manuscript table numbers (126) and draft 17 numbers (127) |
| 128 to 138 | Amendment 48: seven claims (128), ANCOM-BC2 configuration and null (129), dairy comparison (130), Wallace farms (131), CLR transformation grid (132), Charolais leave-one-animal-out (133), power and attenuation (134), mapping classes (135), clade ratio (136), runner (137), number manifest (138) |
| 139 to 145 | Draft 18 and 19 reporting: Figures 2, 4, 5 and Additional file 4 (139), Figure 1 (140), evidence ledger (141), Additional files 1 to 4 (142), draft 18 numbers (143), methane source audit (144), draft 19 numbers (145) |
| 146 | Draft 20 reporting: normal-theory p of the rescaled W beside the empirical permutation p, and draft 20 numbers |
| 147 | Clean-environment reproduction check (docs/REPRODUCTION_CHECK.md) |
| 148 | Sample map: every deposit run to its BioSample, animal and count-table column, with the reason when a run was not analysed (results/sample_map/) |
| lib_replication.R | Shared pooling and labelling functions |
| lib_clr.R | Shared functions for the second method (centred log-ratio with linear models) |
| lib_late.R | Shared functions for the late analyses (Amendment 44), called from 21, 65 and 81 |
| lib_a46.R, lib_a47.R, lib_a48.R | Shared functions for Amendments 46, 47 and 48 |

Long runs were started as background jobs with memory limits, and every job can be resumed: it skips any output that already exists.
