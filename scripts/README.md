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
| 88 to 100 | Calibration with shuffled phenotypes, searches across all genera, farm-by-farm checks, sequencing depth check, second method, sensitivity analyses, summary tables |
| lib_replication.R | Shared pooling and labelling functions |
| lib_clr.R | Shared functions for the second method (centred log-ratio with linear models) |

Long runs were started as background jobs with memory limits, and every job can be resumed: it skips any output that already exists.
