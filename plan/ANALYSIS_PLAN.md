# Pre-specification: multi-cohort replication of rumen microbiome
# associations with feed efficiency in beef cattle

**Version 1.0 — 2026-10-28**
Rushikesh R. Lagad, Aranyak Goswami
Department of Animal Science, University of Arkansas

This document is written and committed **before any read processing or
association testing** on the cohorts described. Every analytical choice below
is fixed as of this commit. Deviations are recorded in Section 13 with dates
and reasons, never by silent edit.

---

## 1. Question

Do published rumen microbiome associations with feed efficiency in beef cattle
replicate in independent cohorts, once each cohort's sensitivity is calibrated
against a design-matched positive control?

Secondary: does moving from DNA-based taxonomy to expressed function change the
answer, where both layers exist?

The study is framed as replication, not discovery. A null result is an
anticipated and interpretable outcome, which is why positive controls and a
power statement are specified in advance rather than added in response to
findings.

---

## 2. Cohorts, frozen

Metadata is frozen at `harmonized/animal_metadata.tsv` and
`harmonized/runs.tsv`, checksums in `harmonized/CHECKSUMS.md5`.

### 2.1 Nelore, PRJNA987743 — primary
- 52 *Bos indicus* steers, Brazil, rumen shotgun metagenomes
- Continuous RFI (`RFI_GC`, `RFI_by_dieta`) and continuous residual methane
  emission (`RME_GC`, `RME_by_dieta`), from Supplementary Table S2
- Diet: Conventional 26, By-products 26
- Slaughter batch: 2 dates, 13 animals per diet per batch (fully crossed)
- `Weight` (Heavy/Light) partially aliased with batch; never modelled jointly
  with batch
- **Join provenance:** S2 `SampleID` equals ENA `sample_alias` exactly, all 52.
  Independently verified by read count: S2 `#Raw_reads` / 2 / ENA `read_count`
  falls in [1.0158, 1.0436] for every animal, a 2.8% band consistent with light
  filtering before deposit. A swapped label would fall outside it.
- `RFI_GC` is centred to zero within slaughter batch, so their contemporary
  group was the batch. It is **not** centred within diet, and diet variances
  differ about threefold (sd 1.86 by-products vs 1.10 conventional). Pooled
  analyses therefore use `RFI_by_dieta` or within-cohort z-scores.

### 2.2 Charolais, PRJEB70895 — cross-breed replication
- 30 Charolais bulls, France, abattoir-sampled rumen content
- Binary extreme RFI only. No continuous RFI exists anywhere in the deposit;
  the `RFI` column holds "Positive"/"Negative" strings on 31 of 80 rows and no
  numeric values.
- Diet (corn silage / grass silage) x year (2019 / 2020), approx. 4 per cell
- Layers: 16S amplicon and mRNA, both present for all 30 animals
- **Exclusions, decided before analysis:**
  - `8513-N-M-2-A20` (`20BAR-01326`): the name encodes N and M, the fields say
    POSITIF and GS, animal 8513 is the only ID appearing in both years, and
    every other 2020 animal is 9xxx or 190827. Its BioSample
    (SAMEA114804964) carries no attributes, so the archive cannot adjudicate.
    Excluded as unresolvable.
  - `8508`: mRNA deposited but no 16S; its mapping row is a Sonde sample.
    Excluded to keep the paired-layer set identical.
  - The 33 mid-RFI (`RFI_1 = "0"`) animals and all Sonde samples have **no
    reads deposited**. A mid-RFI monotonicity control and a within-animal
    sampling-method comparison are therefore not possible in this cohort.
- Animal `8574` has two mRNA runs sharing one sample accession
  (SAMEA115345548). Runs are concatenated, not treated as replicates.

### 2.3 Li et al., PRJNA448333 — baseline, conditional
Named by the handling editor as the baseline. Included with explicit caveats.

- 48 animals, 144 runs (Metagenome, total RNA, mRNA), Lacombe Alberta, 2014
- **Breed** is recoverable: a pure function of the `source_material_id` prefix
  (1-4 = HYB n=28, 5-6 = CHAR n=9, 7 = ANG n=11). Join to archive is 48/48.
- **Phenotypes are not traceable.** BioSample records carry only animal number,
  collection date, location and host species. No RFI, DMI, ADG, breed, sex or
  diet field exists in the archive. The 35 RFI values in
  `rumen448333/input/sample_metadata.csv` come from a source we have not
  identified, and no script in that project builds the file.
- **The deposit does not match the published design.** Li et al. (2019,
  *Microbiome* 7:6) describe 48 steers, 16 per breed, 8 high-RFI and 8 low-RFI
  in each. The deposit partitions 28/9/11, and the local metadata records
  heifers and bulls as well as steers.
- **Breed is perfectly confounded with collection date.** ANG all on
  2014-10-03; CHAR on 09-19 and 09-26; HYB across 07-18 to 08-22. No date
  contains more than one breed. Breed effects in this cohort are not
  identifiable from batch effects. This is a property of the original design,
  not of our reanalysis.
- **Status:** contributes (a) the breed-or-batch positive control, and (b) the
  RFI test **only if** per-animal phenotypes are obtained from the authors.
  A request was sent 2026-10-28. If no traceable phenotype is obtained, this
  cohort contributes the positive control and the reprocessed reads only, and
  the RFI replication rests on the two verified cohorts.

---

## 3. Pre-specified replication targets

The target list is fixed here and will not be extended after results are seen.

### 3.1 Primary target set — Li et al. 2019 (*Microbiome* 7:6)
Their complete reported RFI associations, with the rank at which each was
reported and the direction stated in their paper.

| # | feature | rank | breed reported in | direction | their p |
|---|---|---|---|---|---|
| T1 | unnamed RF19 (text: RF16) | 16S clade | ANG | higher in L-RFI | 0.045 |
| T2 | Firmicutes | phylum | CHAR | higher in L-RFI | 0.042 |
| T3 | Chloroflexi | phylum | CHAR | — | 0.046 |
| T4 | unclassified Methanomassiliicoccaceae | family | CHAR | higher in H-RFI | 0.014 |
| T5 | *Blautia* | genus | HYB | higher in L-RFI | 0.035 |
| T6 | *Methanobrevibacter smithii* | species | HYB | higher in L-RFI | 0.040 |

Diversity claims, tested separately: bacterial Shannon evenness (CHAR,
p=0.03), archaeal Simpson and inverse Simpson (CHAR, p=0.03), bacterial
Shannon (HYB, p=0.04).

**Provenance of these claims, recorded for the discussion.** All six derive
from t tests on arcsine-square-root relative abundances, run within each breed
on 8 high-RFI versus 8 low-RFI animals, across dozens of taxa at phylum and
genus level, with **no multiple-testing correction**. Breed was not a covariate.
Reported p-values span 0.014 to 0.046.

T1 is expected to be untestable: RF19 is a Greengenes 16S clade label with no
counterpart in a genome-based reference database. If it cannot be mapped, that
is reported as a finding, not treated as a missing value.

### 3.2 Secondary target set — Nelore published associations
The 64 features in `PRJNA987743/metadata/Table_7.XLSX`, with their reported
coefficients, standard errors and q-values for `RFI_GC` and `RME_GC`, in rumen
and faeces. Tested for replication in the other cohorts.

### 3.3 Tertiary target set — canonical literature genera
*Methanobrevibacter*, Succinivibrionaceae, *Fibrobacter*, *Ruminococcus*,
*Prevotella*, *Butyrivibrio*. Named in our 2026 review as the recurrent
candidates. Treated as a separate family of six tests.

Each target set is corrected within itself. Sets are not pooled for correction.

---

## 4. Feature space, fixed

- **Primary: genus level.** Short-read species assignment at these depths is
  not reliable, and a 5,467-feature species-level space at n=35 cannot support
  any correction (BH requires p < 1.8e-4 at the top rank; the observed minimum
  was 4.7e-4).
- **Sensitivity: species level**, reported as such.
- Prevalence filter: retained in a cohort if present in at least 10% of that
  cohort's samples.
- Zero replacement: `zCompositions::cmultRepl`, CZM.
- Transform: centred log-ratio, computed **within cohort**, after filtering.
- Targets specified at phylum or family rank (T2, T3, T4) are tested at the
  rank at which they were reported, by aggregating counts to that rank before
  CLR.

---

## 5. Read processing, fixed

Identical for all cohorts, except where the assay forbids it.

1. `fastp`: adapter and quality trimming, per-sample JSON retained as the QC
   record. No separate FastQC/MultiQC step.
2. Host removal: ARS-UCD2.0 for *Bos taurus* cohorts. For Nelore
   (*Bos indicus*), host removal is run against ARS-UCD2.0 **and** against a
   *Bos indicus* assembly, and the retained-read difference is reported. If it
   exceeds 2% of reads, the *Bos indicus* mapping is used for that cohort and
   the choice is reported as a cohort-specific deviation.
3. Taxonomic profiling: Kraken2 plus Bracken.
4. Functional profiling: HUMAnN3 for shotgun layers.
5. Charolais 16S is profiled separately and never merged with shotgun counts;
   it enters only as its own cohort-level effect estimate.

**Database choice is a declared open decision.** The Li et al. reanalysis
returned dozens of *Faecalibacterium* strain-level features among its
top-ranked taxa in rumen content, which indicates a reference database weighted
toward human gut genomes. The database must be rumen-inclusive and is fixed by
amendment to this document before any profiling run, with the chosen release
and date recorded. It is not selected after seeing association results.

---

## 6. Models, per cohort

RFI is standardised to mean 0 and sd 1 **within cohort** before modelling, so
coefficients are per-SD and comparable across cohorts whose RFI derivations and
units differ (Li et al. sd 0.465; Nelore `RFI_GC` sd 1.522).

- **Nelore:** `CLR(taxon) ~ diet + batch + RFI_z`. Marginal terms. Methane
  tested in a separate model of identical form with `RME_z`.
- **Charolais:** `CLR(taxon) ~ diet + year + RFI_class`. 16S and mRNA modelled
  separately, never pooled.
- **Li et al., if phenotypes are obtained:** within-breed tests reproducing
  their design, and a pooled model with a single combined breed-batch term,
  since the two are not separable. Both reported.

Composition-aware differential abundance: ANCOM-BC2 as primary, with ALDEx2 and
LinDA as pre-specified sensitivity methods. ANCOM-BC2 is run with
`p_adj_method = "BH"`. Where `group` is set, structural-zero handling is
reported, because it interacts with group size.

Dispersion is tested with `betadisper` before any location test (PERMANOVA), and
a significant dispersion difference is reported alongside the location result
rather than ignored.

---

## 7. Positive controls, declared in advance

Each cohort must demonstrate that the pipeline can detect a strong effect in
those same animals before its RFI result is interpretable.

| cohort | positive control | status |
|---|---|---|
| Nelore | diet, Conventional vs By-products | to be run |
| Charolais | diet (CS/GS) and year | to be run |
| Li et al. | breed-or-batch (confounded, reported as such) | 42 taxa at q<0.05, min p 9.8e-12, species level |

**Interpretation rule, fixed now:** a cohort's RFI null is reported as
informative only if that cohort's positive control yields at least one taxon at
BH q < 0.05 in the same feature space and model family. If a cohort fails its
own positive control, its RFI result is reported as uninformative rather than
as evidence of absence.

The Li et al. control is a breed-or-batch contrast and will be labelled that way
throughout. It calibrates sensitivity; it does not establish a breed effect.

---

## 8. Multiple-testing correction

- **Replication tests:** Benjamini-Hochberg within each target set separately
  (6 tests for set 1, 64 for set 2, 6 for set 3).
- **Discovery scans:** Benjamini-Hochberg across the full genus-level feature
  space, reported separately from the replication tests and never merged with
  them.
- Nominal uncorrected p-values are reported alongside, explicitly labelled, so
  that comparison with the original publications is possible.

---

## 9. Pooling across cohorts

Random-effects meta-analysis, REML, on per-SD standardised effects.

- Continuous-RFI cohorts contribute standardised regression coefficients.
- Charolais contributes Hedges' *g* from its binary contrast, converted to a
  comparable scale, with the conversion stated and the extreme-group selection
  noted as inflating *g* relative to an underlying correlation.
- Reported for every target: pooled estimate, 95% CI, tau-squared, I-squared,
  and a forest plot.
- Leave-one-cohort-out is run for every pooled estimate.

---

## 10. Replication criteria, fixed before results

A target counts as **replicated** only if all three hold:

1. The effect has the same sign as originally reported.
2. The pooled 95% CI excludes zero after BH within its target set.
3. It survives leave-one-cohort-out, i.e. it is not carried by a single cohort.

A target counts as **not replicated** if the pooled CI includes zero in a
cohort set whose positive controls passed.

A target counts as **untestable** if it cannot be mapped into the feature space
(expected for T1) or if every cohort holding it failed its positive control.

---

## 11. Power, computed

Two-sided alpha 0.05, correlation test.

| n | min detectable r at 80% power | power at r=0.2 | r=0.3 | r=0.4 |
|---|---|---|---|---|
| 30 (Charolais) | 0.492 | 0.18 | 0.36 | 0.60 |
| 35 (Li et al.) | 0.458 | 0.21 | 0.42 | 0.67 |
| 52 (Nelore) | 0.380 | 0.29 | 0.58 | 0.84 |
| 87 (pooled, if Li et al. joins) | 0.296 | 0.46 | 0.81 | 0.97 |

Charolais, simulated as 15 vs 15 drawn from the top and bottom 16% of 100
animals, has power 0.36 at underlying r=0.20, 0.69 at 0.30 and 0.92 at 0.40.
Extreme-group selection concentrates the contrast, so this 30-animal arm has
more power for a true effect than the 35-animal baseline.

**Stated consequence.** No cohort here can detect r below about 0.3 reliably.
Any null reported in this study is a null with respect to effects of roughly
r >= 0.3, and will be stated in those terms. Claims of absence beyond that
range are not supported by these data and will not be made.

---

## 12. What will not be done

- No target added to any set after results are seen.
- No cohort added or removed after results are seen. PRJEB31266 is excluded for
  reasons of scale and heterogeneity and will not be added later.
- No change of primary feature space after seeing results. Genus is primary;
  species is sensitivity.
- No change of correction method after seeing results.
- No subgroup or interaction analysis not named in this document.
- No reporting of uncorrected p-values as evidence of association.
- Li et al. phenotypes will not be reconstructed by inference, rank-matching, or
  any procedure other than receipt from the authors or a citable published
  source.

---

## 13. Deviations log

| date | section | change | reason |
|---|---|---|---|
| | | | |

---

## 14. Data provenance record

| file | md5 | source |
|---|---|---|
| harmonized/animal_metadata.tsv | see CHECKSUMS.md5 | this pipeline |
| harmonized/runs.tsv | see CHECKSUMS.md5 | this pipeline |
| PRJNA987743/metadata/Table_2.XLSX | | Europe PMC, PMC10754987, supplement S2 |
| PRJEB70895/metadata/mapping.tab | 3e87555a73920feb3648f19a25f706c7 | Dataverse doi:10.57745/ACNRVQ |
| rumen448333/metadata/biosample_attrs.tsv | | ENA BioSample XML, PRJNA448333 |
| rumen448333/input/sample_metadata.csv | | **PROVENANCE UNKNOWN** — not used for any phenotype claim until resolved |

---

## Amendment 1, 2026-10-29: database decision (Section 5)

Kraken2 Standard, build k2_standard_20240605, 49,138 taxa, 78 GB, fixed as the
profiling database for all cohorts. Composition recorded before any profiling:

| rumen genus | refs | gut genus | refs |
|---|---|---|---|
| Prevotella | 34 | Bacteroides | 57 |
| Fibrobacter | 8 | Faecalibacterium | 22 |
| Ruminococcus | 11 | Akkermansia | 11 |
| Methanobrevibacter | 11 | | |
| Succinivibrio | 3 | | |
| Methanomassiliicoccus | 1 | | |
| Ruminobacter | 0 (absent) | | |
| Anaerovibrio | 0 (absent) | | |

The build is weighted toward human gut genomes. It is adopted because it is
applied identically to every cohort, which preserves the validity of a
replication comparison at reduced sensitivity, and because a rumen-specific
build is a multi-day rebuild that may not change the conclusion.

**Declared consequences, before seeing any result:**
1. Targets T4 (unclassified Methanomassiliicoccaceae, 1 reference genome) and
   the canonical Succinivibrionaceae (3) are flagged a priori as
   low-confidence. A null for either is reported as uninformative about the
   taxon, not as failure to replicate.
2. Ruminobacter and Anaerovibrio are absent and will not be reported at all.

**Pre-specified adequacy test.** After the first cohort is profiled, the median
per-sample classified-read fraction is computed. If it falls below 60%, the
database is judged inadequate for rumen content, this amendment is superseded,
and a rumen-inclusive build is constructed before any association testing. The
threshold is fixed here, before any profiling output exists.

## Amendment 2, 2026-10-29: pipeline deviations (Section 5)

- minimap2 -x sr replaces bowtie2 for host depletion. Index builds in 37 s
  against several hours, accuracy for depletion of short reads is equivalent.
- The ARS-UCD2.0 versus UOA_Brahman_1 comparison runs on 8 Nelore samples
  rather than all 52. The quantity estimated is a per-genome mapping rate,
  which 8 samples resolve well inside the 2% decision boundary.


## Amendment 3, 2026-10-29: host reference and the as-deposited state of PRJNA987743

**Reference decision.** ARS-UCD2.0 (Bos taurus) is used for every cohort.
Comparison on 8 Nelore samples, 400,000 read pairs each: mean retained-read
difference between ARS-UCD2.0 and UOA_Brahman_1 was 0.005%, maximum 0.018%,
far below the 2% threshold set in Section 5. A single reference across cohorts
also removes a between-cohort technical factor.

**PRJNA987743 was deposited already host-depleted.** The ENA read count equals
Supplementary Table S2's `#Non-host reads` divided by two for all 52 rumen
samples, mean ratio 1.0000 with standard deviation 0.0000, while the ratio to
`#Raw_reads` is 0.9771. The authors report 2.08% host removal, and our pipeline
finds a further 0.4% median residual. The submitters deposited filtered reads.

Consequences, recorded before analysis:
1. Host depletion is close to a no-op for Nelore. The step is still run, so the
   procedure stays uniform across cohorts, but the cohorts differ in their
   as-deposited state. Li et al. and Charolais deposits have not been checked
   and may contain raw reads.
2. Non-host read yield is therefore NOT a comparable QC metric across cohorts
   and will be reported per cohort, never pooled or compared between them.
3. The reference comparison above is uninformative about whether ARS-UCD2.0 and
   UOA_Brahman_1 differ on raw Bos indicus rumen data, because only ~0.3% host
   DNA remained for them to differentiate. The decision stands, but it rests on
   there being nothing left to separate rather than on demonstrated equivalence.
   If a future cohort supplies raw Bos indicus reads, the comparison is rerun.


## Amendment 5, 2026-10-30: adequacy test, corrected and complete

Amendment 4 reported 17.0, 19.4 and 16.0 percent from the first three samples.
The full 52-sample Nelore profile with Kraken2 Standard gives a classified-read
fraction of **median 16.9%, range 15.3 to 19.4%**. The Amendment 1 threshold of
60% FAILS.

A reporting bug in 08_profile.sh initially printed this as "median 100.0%, PASS".
The adequacy calculation globbed kraken/*.report, which also matched the Bracken
reports; those carry no `unclassified` line, so the code recorded them as 100%
classified and the median fell on the Bracken files rather than the Kraken2 ones.
The script has been corrected to exclude Bracken outputs. The erroneous PASS was
never acted on.

**The failure is depth, not misassignment.** Genus composition of the classified
fraction is characteristically ruminal:

| taxon | % of classified | taxon | % of classified |
|---|---|---|---|
| Prevotella | 11.04 | Bacteroides | 2.18 |
| Fibrobacter | 2.24 | Faecalibacterium | 0.98 |
| Methanobrevibacter | 2.49 | Ruminococcus | 0.96 |
| Butyrivibrio | 1.17 | | |

Prevotella exceeds Bacteroides fivefold and Fibrobacter exceeds Faecalibacterium
by 2.3x, the opposite of what database-driven misassignment to human gut taxa
would produce. The classified minority is real rumen signal; the database simply
lacks the uncultured majority of the community.

Consequences:
1. A rumen-inclusive database is still required before association testing, per
   Amendment 1. The threshold is not revised.
2. The Standard-database profile (2,725 genera, 11,605 species across 52 samples)
   is retained as a pre-specified sensitivity analysis rather than discarded. The
   classified fraction varies by only 4 percentage points across samples, so the
   loss is approximately non-differential and between-sample CLR comparisons on
   that subset remain interpretable.
3. Database route selected: GTDB release 207 (62,443 genomes) from the Struo2
   distribution, built locally with a size cap, since the pre-built hash.k2d is
   287 GB against 122 GB of available RAM.

