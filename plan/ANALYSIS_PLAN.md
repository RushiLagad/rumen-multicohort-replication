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


## Amendment 6, 2026-10-30: contamination screen, pre-specified

Supplementary Table 7 of the Nelore source paper (Front Microbiol 2023,
PMC10754987) reports 64 features associated with RFI or residual methane. Their
habitat distribution is asymmetric:

- Rumen / RFI_GC, 22 genera: 17 are soil or environmental organisms
  (Actinomadura, Actinoplanes, Anaeromyxobacter, Bosea, Caulobacter,
  Cellulomonas, Comamonas, Gordonia, Lysobacter, Nonomuraea, Pseudonocardia,
  Rubrobacter, Sinorhizobium, Sorangium, Sphingobium, Brevibacterium,
  Burkholderia); 4 are genuine rumen taxa.
- Feces / RFI_GC, 7 genera: 5 are marine (Aliivibrio, Idiomarina, Marinomonas,
  Photobacterium, Psychromonas).
- Feces / RME_GC, 22 genera: include deep-sea, hot-spring, freshwater bloom and
  insect-endosymbiont organisms.
- Rumen / RME_GC, 13 genera: biologically coherent, 12 of 13 Lachnospiraceae or
  related anaerobes.

Several rumen RFI genera (Bosea, Burkholderia, Caulobacter, Comamonas,
Sphingobium, Brevibacterium, Gordonia, Pseudonocardia) are documented reagent
and extraction-kit contaminants (Salter et al. 2014, BMC Biology 12:87).

**Screen, criteria fixed before examination.** For each Table 7 feature, in our
harmonized Nelore profile: (a) median relative abundance across 52 samples;
(b) Spearman correlation of CLR abundance with total non-host depth;
(c) association with slaughter batch. A feature is flagged
contamination-consistent if median relative abundance < 0.01% AND it shows a
negative depth correlation or a batch association at BH q < 0.05 across these
64 tests.

Applied to the published feature list only. It does not alter the replication
targets fixed in Section 3, and is reported as a property of the source study.


## Amendment 7, 2026-10-30: contamination screen result, NEGATIVE

The Amendment 6 screen was run against the Nelore Standard-database genus
profile (2,725 genera x 52 samples, 26 animals per slaughter batch).

**The pre-specified criterion is not met. No Table 7 feature is flagged
contamination-consistent.**

- Relative abundance: rumen RFI genera median 0.065%, rumen RME genera median
  0.150%. **Zero of 21 RFI genera fall below the 0.01% threshold.** Burkholderia
  0.189%, Selenomonas 0.241%, Comamonas 0.116%. These are not the trace levels
  characteristic of reagent contamination in low-biomass samples.
- Depth correlation: the implementation correlated log raw counts with total
  depth rather than CLR-transformed abundance, so all 34 genera returned rho
  0.54 to 0.82 at p < 1e-4 as an arithmetic consequence. **This test was invalid
  and is reported as uninformative, not as evidence either way.**
- Batch association: 12 of 21 rumen RFI genera versus 2 of 13 rumen RME genera
  reach nominal p < 0.05 (Fisher exact OR 7.3, p = 0.030). **No individual
  genus survives BH correction across the 34 tests (minimum q = 0.099).**

**Conclusion recorded as specified.** The hypothesis that the source study's RFI
associations reflect reagent contamination is not supported. The habitat-based
observation in Amendment 6 stands as an observation about the published feature
list, and the group-level batch asymmetry is reported with its effect size and
its failure to survive correction. Neither is claimed as a contamination finding.

**Alternative hypothesis, to be tested with the GTDB database.** The soil and
marine genera may instead reflect misassignment: with Kraken2 Standard
classifying only 16.9% of rumen reads, unclassifiable sequences from uncultured
rumen organisms can be assigned to the nearest represented relative. If so,
their apparent abundance should fall substantially when profiled against a
database with genuine rumen coverage. This is a pre-specified comparison:
for the 64 Table 7 features, relative abundance under Kraken2 Standard versus
GTDB r207, same samples, same pipeline otherwise. A median reduction exceeding
50% in the environmental genera, with the Lachnospiraceae set stable, supports
misassignment.


### Amendment 8 (2026-10-30): GTDB r207 database build — defects found and corrected
1. The Struo2-distributed `seqid2taxid.map` did not match the Struo2 `library.tar.gz` headers (0 of 7.7M contigs assigned; the 15:13 build produced an 80 GB table containing 0 sequences). Map regenerated from the 62,444 `prelim_map*.txt` files shipped inside the library; deduplicated to 7,745,746 contigs.
2. kraken2 (conda build reported as 2.17.2) silently dropped every sequence whose leaf taxid exceeded 2^31 (3-genome test: 362/519 counted; 31,330 of 62,443 genomes affected). All taxids renumbered to sequential integers (root = 1, max 401,815) in nodes.dmp, names.dmp and the map; `scripts/10b_renumber_taxids.py`; old-to-new table at `taxonomy/taxid_remap.tsv`; originals in `taxonomy/orig/`. Verified 519/519 on the 3-genome test and 62,443/62,443 on a one-stub-per-genome build.
3. Library used as shipped (`*.fna`, unmasked); the `.fna.masked` twins carry different headers and no masking characters and were not usable. The Standard DB is masked; noted as a difference in the Standard-vs-GTDB comparison.
4. Capacity: estimate_capacity 7.70e10 distinct minimizers; cap 85 GB = 2.125e10 cells. The DB therefore retains a hash-subsampled minority of minimizers; adequacy is judged by the pre-registered 60% test (Amendment 1), not assumed.
5. systemd-oomd killed the 11:48 build (pressure-based kill of the user session); disabled on the build host.
6. Acceptance gate added (`scripts/11_verify_db.sh`): processed sequences must equal map lines; hash.k2d > 80 GB; kraken2-inspect and a smoke classification before the 5-sample adequacy test. Final build launched 15:56, build_db from 16:05:42.

### Amendment 9 (2026-11-01, written with 3 of 52 Nelore samples profiled: 54.1, 47.1, 51.3% classified): adequacy branch and genus match keys
**Context.** The capped GTDB r207 database (Amendment 8) retains roughly a quarter of library minimizers. A 136 bp read carries ~40 minimizers, so reads from genomes present in the library are still classified; the unclassified fraction is dominated by reads from lineages absent from GTDB r207 representatives or too divergent for exact 35-mer matches. A full-minimizer table recovers part of the latter (Smith et al. 2022). The pre-registered test (Amendment 1: median classified fraction across Nelore samples ≥ 60%) is evaluated on all 52 Nelore samples from the run started 2026-11-01 16:39.

**Branch, fixed before the result is known.**
- If median ≥ 60%: the capped database is the primary database. No further build.
- If median < 60%: two-pass, full-resolution database.
  1. Screen = the completed capped run over all 131 samples (3 cohorts).
  2. Genome set = every library genome whose species has ≥ 10 clade-level reads in the Kraken2 report of at least one sample, any cohort. Threshold fixed here; not tuned.
  3. Rebuild from that genome set with no `--max-db-size` cap (full minimizers), same k/l/taxonomy (renumbered, Amendment 8), Bracken distributions at 136, 150 and 100 bp.
  4. Acceptance: all selected genomes counted (processed sequences = contigs of the set); the Amendment 1 test re-run on the new database and reported whether it passes or not. If it also fails, the full-resolution database is still primary (it dominates the capped one on sensitivity) and the shortfall is reported as a limitation; no third build.
  5. The capped-database profiles are retained as a sensitivity analysis (database-choice comparison alongside Standard vs GTDB).
  Rationale for the subset: reads from excluded genomes cannot contribute meaningfully, since the screen saw all 62,443 genomes with enough minimizers to detect presence at the 10-read level; exclusion removes capacity, not signal. Option rejected: memory-mapping the 287 GB Struo2 table (download time, random-access throughput, and an unverified build).

**Genus match keys (GTDB suffixes).** GTDB splits NCBI genera into suffixed genera (e.g. Ruminococcus_B/_E/_G, Methanobrevibacter_A/_C, Treponema_D/_F). Rules: (a) the primary feature space keeps GTDB genera as named; (b) a published target genus named in NCBI nomenclature matches the set of GTDB genera whose name equals the target or the target followed by `_` and a suffix, with counts summed across that set for the replication test; (c) phylum-level targets use current names (Bacillota for Firmicutes, Chloroflexota for Chloroflexi) and their GTDB equivalents where GTDB names differ (Bacillota = Firmicutes + Firmicutes_A/_B/_C… in GTDB r207, summed). The target-to-key table is written to `plan/target_keys.tsv` before association testing and committed.

### Amendment 10 (2026-11-01, before any association test): target key table committed
`plan/target_keys.tsv` built by `scripts/12_target_keys.py` from the Amendment 9 rule against the names present in k2_gtdb207. Result: 76 targets; 71 mapped; 4 Nelore genera with no GTDB r207 counterpart (Sorangium, Novisyntrophococcus, Saccharothrix, Desulfosediminicola) are untestable per Section 3 and reported as such; T1 untestable by design.
1. **Site strata.** 29 of the 64 Nelore Table 7 associations were reported in faeces. All 64 remain in the set and are corrected together, as pre-registered. Results are reported in two strata, rumen-site (35) and faecal-site (29); a faecal association tested in rumen data is a cross-site claim and the headline replication count for the Nelore set uses the rumen stratum.
2. **Composite keys.** Where the suffix rule maps one target to several GTDB taxa, counts are summed before CLR and the result is reported as the composite (e.g. Eubacterium = 14 GTDB genera; Firmicutes = Firmicutes + Firmicutes_A..H; Methanobrevibacter smithii = Methanobrevibacter_A smithii + smithii_A). Composite membership is fixed in the table.
3. Targets appearing in more than one set (Blautia in Li T5 and Nelore N40; Eubacterium for RFI and RME) are separate tests, corrected within their own set, as Section 3 states.

**Amendment 9 result (2026-11-01 22:02):** capped GTDB r207, Nelore n=52, median classified 53.0% (range 45.1–64.5%). FAIL against 60%. Two-pass branch executed as specified (`scripts/13_twopass_db.sh`); outputs in `profiles_gtdb_full/`, capped profiles retained in `profiles_gtdb/`.

**Amendment 11 (2026-11-02, before the second-pass build):** The Amendment 9 screen (species with >=10 clade reads in any sample) selected 58,360 of 62,443 genomes (93%) and so cannot shrink the library to fit the 96 GB table. This was found from the selection count alone, before any second-pass build or adequacy test, so the correction does not depend on the outcome. GTDB r207 holds one representative genome per species, so species and genome counts coincide. The threshold is replaced by the smallest power of ten whose selection fits the table by genome count (at most 15,000 genomes, 24% of the library, against a table that holds about 31% of the full minimizer estimate). Maximum clade reads per species across the 131 samples: >=10 gives 58,360; >=100 gives 12,111; >=1,000 gives 4,802; >=10,000 gives 2,016. The rule is therefore >=100 clade reads in any sample (12,111 genomes). All other terms of Amendment 9 stand, including the 60% median adequacy test, a single rebuild only, and reporting the shortfall if the test fails. The first extraction under the old threshold was aborted before pruning.

**Amendment 12 (2026-11-02, before the second-pass profile and before any association test):** Under the Amendment 9/11 database, 29 Nelore Table_7 genera (N05, N07, N09, N12, N14, N16, N18, N24, N25, N26, N27, N29, N43, N44, N46, N48, N50, N51, N52, N53, N54, N56, N57, N58, N59, N60, N61, N63, N64) have genomes in GTDB r207 but none among the 12,111 selected, because no species in them reached 100 clade reads in any of the 131 samples on the capped database. Being unmappable into the primary feature space, they are untestable under the existing Section 3 definition. They are reported as their own class, "untestable: below detection on GTDB", distinct from the 4 genera with no GTDB counterpart and from T1. The replication rate is reported twice: among testable targets (42 of the 71 mapped), and as a lower bound with these 29 counted as not replicated. Their detection on the capped database, which holds all 62,443 genomes, is reported descriptively (maximum reads and prevalence per cohort) and serves as the Amendment 7 misassignment comparison; no association test is run on them. Coverage of every target is in `plan/target_coverage.tsv` (`scripts/14_target_coverage.py`).

**Amendment 13 (2026-11-02, before any positive control or association test was run):** Data handling and interpretation fixed in advance. (1) Charolais animal 8574 has two mRNA runs; their counts are summed as technical replicates, giving 30 animals and 30 mRNA samples. (2) The Charolais positive control (Section 7, "diet (CS/GS) and year") is passed if at least one genus has BH q < 0.05 for the diet term or the year term of `diet + year + RFI_class`; both terms are reported. The Nelore control is the diet term of `diet + batch + RFI_z`. (3) Controls use ANCOM-BC2 (`p_adj_method = "BH"`, `group` set to the contrast, structural zeros on, `prv_cut = 0.10` matching the plan's prevalence filter, other settings at package defaults) and, for the same terms, PERMANOVA (marginal terms, 9,999 permutations, Aitchison distance on CZM-replaced CLR) and betadisper; seed 20261102. (4) The Nelore RFI column is `rfi_value` (sd 1.529 in the 52 animals, against 1.522 quoted in Section 6); RFI is z-scored within cohort, so the difference does not enter the model. (5) Primary controls use the full-resolution profiles (`profiles_gtdb_full`); the capped profiles (`profiles_gtdb`) are the sensitivity set and are used for development runs. Script: `scripts/20_positive_controls.R`.

**Amendment 9 result, second pass (2026-11-02 14:25):** Full-resolution GTDB r207 database built from the 12,111 genomes selected under Amendment 11 (1,522,143 contigs, 54 GB table, uncapped, verification passed). Nelore, n=52: median classified fraction 61.8% (range 53.9 to 72.6%). PASS against the 60% threshold of Amendment 1. This database gives the primary profile set (`profiles_gtdb_full`). The capped GTDB profiles (median 53.0%, FAIL) and the Standard RefSeq profiles (16.9%, FAIL) are retained as sensitivity sets, and all three results are reported. No third build is made.

**Positive control results (2026-11-02), primary profiles (`profiles_gtdb_full`):** Nelore, diet: 637 of 2,887 genera at BH q < 0.05 (min q 1.7e-4), PERMANOVA R2 0.040, p 0.001, betadisper p 0.35: PASS. Charolais mRNA, diet: 456 of 2,761 genera (min q 1.8e-8), R2 0.106, p 1e-4, betadisper p 0.82; year: 113 genera (min q 3.7e-6), R2 0.067, p 1e-4, betadisper p 0.50: PASS. The capped profiles (sensitivity set) also pass for both cohorts; in the capped Charolais year contrast the dispersion test is significant (betadisper p 0.006) and is reported alongside. The Li et al. control awaits breed labels. Under Section 7 the Nelore and Charolais RFI results are therefore interpretable.

**Amendment 14 (2026-11-02, before any replication test was run):** Specification where Section 10 is silent; no new criterion for "replicated" is added. (1) Effect scale: per target, cohort and RFI term, t = lfc/se from ANCOM-BC2 in the Section 6 model (RFI_z in Nelore, high versus low RFI_class in Charolais); effect = partial r = t/sqrt(t^2 + df), df = n minus fitted coefficients; pooled on Fisher z with variance 1/(n - 3 - q), q = other covariates. Hedges' g is also reported for Charolais; extreme-group inflation is stated, not corrected. Sensitivity: t from lm on CLR. (2) Sign: positive means more abundant with higher RFI, so "higher in L-RFI" and "lower with higher RFI" are negative, and "higher in H-RFI" and "higher with higher RFI" positive. T3 has no stated direction, so criterion 1 is not applied to it and its label says so. (3) Cohort roles: Li-set targets (T2 to T6) and canonical targets are tested in Nelore and Charolais mRNA (and Li if phenotypes arrive). Nelore-set targets are tested in Charolais mRNA only, Nelore being their discovery cohort; the Nelore re-analysis is reported as reproducibility and is not pooled. Of the 64, the 21 mapped rumen RFI targets form the replication stratum; the 7 mapped faeces RFI targets are tested in rumen data and reported as a separate cross-site stratum; the 32 mapped RME targets cannot be tested outside Nelore (no methane data) and are classed "untestable: phenotype not available", the 12 rumen ones being reported for Nelore reproducibility only. (4) Pooling: REML with Wald z as in Section 9; with k = 2 the fixed-effect and Hartung-Knapp estimates are shown as sensitivity. With k = 1 there is no pooling; a target whose single estimate has the stated sign and BH-adjusted p < 0.05 is labelled "single-cohort support", not "replicated", since criterion 3 cannot be evaluated. Criterion 2 is met when the BH-adjusted pooled p < 0.05. (5) BH uses the Section 8 family sizes (6, 64, 6), untestable targets entering with p = 1; the same correction over testable targets only is a sensitivity. (6) Charolais 16S is reported beside the mRNA result and not pooled, as the 30 animals are the same. (7) Any "not replicated" is stated as a null for effects of roughly r >= 0.3 (Section 11).

**Amendment 14b (2026-11-02): classification rules, committed before any target estimate.**
1. Cohort roles. Sets L and C: Nelore and Charolais mRNA (expected k=2). Set N: Charolais mRNA only (k=1). Faeces-site targets are tested on rumen profiles and labelled cross-site. RME targets and the 29 targets of Amendment 12 are untestable.
2. Model. ANCOM-BC2 on the primary profile (profiles_gtdb_full), prv_cut=0.10, group=NULL, struc_zero=FALSE, other arguments at package defaults. Term of interest: rfi_z (Nelore), high_RFI vs low_RFI (Charolais). Genus targets are collapsed over their gtdb_keys before fitting; family, phylum and species targets are fitted on the table aggregated to that rank. Targets with identical key sets share one feature. Partial overlaps stop the script and are resolved by amendment.
3. Missing estimates. A target removed by the 10% prevalence rule, or given no estimate, in a cohort contributes nothing there. k is the number of cohorts with an estimate. BH family sizes (6, 64, 6) do not change; untestable targets enter with p=1.
4. Labels. Concordant: pooled sign equals published sign (positive = more abundant with higher RFI). Significant: BH q<0.05 within the set.
- replicated: k>=2, concordant, significant, criterion 3 holds.
- single-cohort support: k=1, concordant, significant; or k>=2, concordant, significant, criterion 3 fails.
- opposite sign: significant but not concordant. Counts as not replicated, reported separately.
- not replicated: every other testable target.
- untestable: k=0, unmapped, or no phenotype.
5. Criterion 3. For k>=2 each cohort is left out in turn and the target re-estimated from the rest. Concordance and BH q<0.05 (same family size) must hold in every leave-one-out estimate. At k=2 each cohort alone must be significant. A nominal variant (p<0.05, uncorrected) is a sensitivity column and sets no label. k=1 targets have no criterion 3 and cannot exceed single-cohort support.
6. Pooling as Amendment 14. Hartung-Knapp and fixed-effect p-values at k=2 are sensitivities and set no label.

**Amendment 14c (2026-11-02): reporting additions, committed before any target estimate. Labels of Amendment 14b are unchanged.**
1. Scale. Nelore uses a continuous contrast (rfi_z), Charolais a high vs low RFI class contrast. Partial r from the two designs is not on one scale. Cohort-specific estimates are primary; the pooled r is labelled approximate. Labels rest on direction and BH significance only.
2. BH sensitivity. q-values are also computed within testable targets only, and reported beside the pre-registered q. Labels use the pre-registered q.
3. Cross-site. Faeces-site targets carry a cross_site flag and are not described as replications of the faecal association.
4. p-value check. ANCOM-BC2 p-values are reported beside the t-based p. A tenfold disagreement is flagged.
5. Large-effect flag. A not-replicated target whose 95% CI includes both 0 and |r|=0.3 in the published direction is flagged, so absence of replication is not read as absence of effect.
6. Assay. Nelore is DNA shotgun and Charolais is mRNA. Heterogeneity is not attributed to breed or diet alone.
7. Few-feature ranks. If the phylum or family ANCOM-BC2 logs show convergence problems or implausible estimates, those targets are reported as exploratory and a CLR-based check is added by amendment before interpretation.

**Amendment 14d (2026-11-02), committed before any target was classified.**
1. Targets without a stated direction (T3 and canonical C1 to C6) follow the T3 rule of Amendment 14: criterion 1 is not applied. Concordance is replaced by agreement of sign: the pooled estimate and every leave-one-out estimate must share one sign. Labels carry "(direction not stated)". "Opposite sign" cannot occur for them.
2. Disclosure. The Nelore phylum fit log, printed to check that ANCOM-BC2 ran, showed T2 (W 0.24) and T3 (W -0.31), both null. No other target estimate had been viewed before this amendment.
3. Report columns. coverage is the minimum over cohorts of keys present divided by keys listed; partial_coverage flags values below 1. Hedges' g for Charolais uses the adjusted t with group sizes 14 and 16.

**Amendment 15 (2026-11-02): Level 1 reproducibility, added after the replication results were seen and committed before it was run.**
1. Purpose. The 64 Nelore-set targets (Section 3.2) were tested in Charolais only. To read that result, the same targets are re-tested in Nelore, where they were discovered. This is a reproducibility check of the published associations under this pipeline, not an independent test, and is reported as such. It is not specified in Sections 1 to 14 and is declared post hoc with respect to the replication results.
2. Targets and family. The same 64 targets, no additions. BH over the same family of 64, untestable targets entering with p = 1. A correction over testable targets only is a sensitivity. Faecal-site targets are tested on rumen profiles and labelled cross-site, as in Amendment 14b.
3. Models. ANCOM-BC2 on profiles_gtdb_full with the settings of Amendment 14b item 2. RFI_GC targets: diet + batch + RFI_z. RME_GC targets: diet + batch + RME_z, with RME_z the standardised methane_value. Each target is tested once, for its own phenotype.
4. Labels. reproduced: BH q < 0.05 and the published sign. opposite sign: q < 0.05 and the other sign. not reproduced: every other testable target. untestable: below detection, filtered, or no genome-based counterpart.
5. Reporting. Proportion of testable targets with the published sign, with an exact binomial 95% CI. Spearman correlation between the published coefficient and this estimate, descriptive only because the scales differ. No uncorrected p-values are reported as evidence (Section 12).
6. Caveat stated in advance. The targets were selected from these same animals, so sign agreement is expected to exceed chance even under no true effect. A high reproduction rate is not independent confirmation.
7. Level 3 (transfer) is not specified here. Any such analysis gets its own amendment before it is run.

**Amendment 16 (2026-11-02): the Amendment 7 misassignment comparison, operationalised and committed before it was run.**
1. Data. Nelore genus profiles for the same 52 samples: Kraken2 Standard (profiles/) and capped GTDB r207 (profiles_gtdb/). The capped database holds every genome, so a drop is a real test. The two-pass database is not used: it holds only genomes the screen found, so absent genera would fall to zero by construction (Amendment 12).
2. Metric. Per feature, the median across samples of the genus share of classified reads in each profile (the metric of Amendment 6). Standard: the row whose name equals the Table 7 name. GTDB: summed over the gtdb_keys of plan/target_keys.tsv. Reduction = 1 - GTDB median / Standard median; a feature absent in GTDB has reduction 1. Features with no Standard counts and the unmapped genera are listed and excluded.
3. Groups, from Amendment 6. Environmental: the 17 rumen RFI soil or environmental genera and the 5 marine faecal genera it names. Stable comparison set (Amendment 7 "Lachnospiraceae set"): the rumen RME targets. Other groups are reported, not tested.
4. Criterion. Amendment 7 left "stable" undefined; it is defined here. Misassignment is supported if the median reduction in the environmental group exceeds 50% and the median reduction in the stable set is below 25%. Medians carry bootstrap 95% CIs over features (2,000 resamples, seed 20261102). Otherwise the result is reported as not supporting, with the full distribution.
5. Secondary, descriptive. The two databases classify different fractions of reads (16.9% vs 53.0%), so shares of classified reads can move for reasons unrelated to misassignment. The difference between the two group medians is reported beside the primary result.
6. Limit. A drop supports misassignment as an explanation for the source study's environmental genera. It does not show that the source study's associations are false.

**Amendment 16 result (2026-11-02): the pre-specified criterion is NOT met.** Environmental group (22 named, 21 included): median reduction 94% (bootstrap 95% CI 90 to 96%); all 21 are reduced by at least 63%. Stable set (13 rumen RME targets, 12 included): median reduction 85% (CI 22 to 95%), above the 25% limit, so the criterion fails. The difference between the group medians is 9 points. Stable-set values range from an increase of 185% (Eubacterium) to a reduction of 97%, and three of 12 increased. The comparison therefore does not distinguish misassignment from database and taxonomy change: shares of classified reads fall when classification spreads over 12,423 genera instead of 2,725, and GTDB reorganises Lachnospiraceae genera, which plausibly explains the stable set's spread. The stable set also contains one marine genus (Winogradskyella, reduction 76%). Reported as not supported by the pre-specified criterion, and not refuted. Per-feature table: results/misassignment/standard_vs_gtdb.tsv. No further comparison is planned; any later analysis of this question is exploratory and needs its own amendment before it is run.

**Amendment 17 (2026-11-02, before any Charolais 16S read was processed):** Specification of the Charolais 16S layer, which the plan defines only as separate and unpooled. (1) Role: the same 30 animals as the Charolais mRNA layer, so it is an assay sensitivity check on that estimate. It is never pooled with the mRNA layer, with Nelore or with shotgun counts, receives no replication label and does not enter the Section 10 classification. (2) Reads: the 341F primer CCTAYGGGRBGCASCAG and the 806R primer GGACTACHVGGGTWTCTAAT, each preceded by a variable spacer, are removed with cutadapt (-e 0.12, -O 12, pairs lacking either primer discarded, minimum length 100). Only forward reads are analysed: after trimming, reads of about 230 and 225 bp must span an insert of about 428 bp, an expected overlap of about 25 bp, so merging would discard pairs by length. Merged pairs are not analysed. (3) DADA2 with the default error model, per-sample inference and consensus chimera removal; maxN 0, maxEE 2, truncQ 2, PhiX removed. Forward truncation length is the last position (capped at 220) at which the median Q pooled over the 30 samples is at least 30. It uses read quality only, and no phenotype file is loaded before the genus table exists. (4) Taxonomy: assignTaxonomy, minBoot 50, against the local SILVA 138.1 nr reference in /data/references/silva (mothur format, converted to a DADA2 training file without edits). (5) Genus names: brackets removed, a trailing " group" removed, a trailing space plus digits removed. Names ending in _ge, and "uncultured", "unidentified" and unassigned, are dropped, and only classified reads enter the table. A target key matches only by exact equality with these names. GTDB suffixed keys (for example Eubacterium_A) do not match, so coverage is partial for such targets and is reported as such. (6) Ranks: genus and family only. Family comes from the same genus lineage used for the other layers. Species (not resolved by about 230 bp of V3-V4) and phylum (SILVA and NCBI phylum names differ) are not tested. (7) Estimation is unchanged from Amendment 14b (ANCOM-BC2, diet + year + RFI_class, prv_cut 0.10, partial r, df 26), written to results/charolais16s. Reported per target: coverage, r, p, sign agreement with the mRNA estimate, and the rank correlation of r between the two layers. No threshold is added. Any change to these rules after the genus table exists needs a new amendment.

**Amendment 18 (2026-11-02, exploratory, before any leave-one-animal-out fit was run):** Influence check on the Charolais mRNA estimates. It is exploratory, is not part of the pre-registered classification, and creates or changes no label. (1) Design: each of the 30 animals is removed in turn and the Amendment 14b fit is repeated unchanged (same model, filters and target definitions; n 29, df 25) for the genus, family and phylum targets. Species targets are omitted because each species fit is slow. Nelore is not repeated because the question concerns the n 30 cohort. (2) Reported per target: the full-data r, the range of r over the 30 fits, the largest absolute change and the animal causing it, the number of fits whose sign differs from the full-data sign, and, for targets with p below 0.05 on all animals (exact t on W, as for k = 1), the share of fits keeping p below 0.05. (3) No threshold is added. Targets whose sign flips under any single deletion are listed in the supplement as influenced by single animals, and the paper states that this check was exploratory. (4) Any further use of these fits needs a new amendment.

**Amendment 17 result (2026-11-03):** The forward-read DADA2 run (truncation 220) kept 98.5% of the 2,434,300 input read pairs through primer removal and gave 8,628 ASVs; after the naming rule 187 genera remain, and a median of 66% of non-chimeric reads are classified to genus. Of 34 genus targets, 11 matched a 16S genus by exact name and 9 were estimable in 16S; 22 were estimable in the Charolais mRNA layer and 9 in both. Among those 9, the signs agreed for 4 (4.5 expected by chance), the rank correlation of r between layers was 0.25, and none had p below 0.05 in 16S. At family rank 1 of 2 targets matched (C2: 16S r about 0.17, mRNA r 0.35, same sign). The 16S layer neither supports nor contradicts the mRNA estimates: coverage is low and both layers have little power at n 30. It is reported descriptively and assigns no label. Per-target table: results/charolais16s/comparison_with_mrna.tsv.

**Amendment 18 result (2026-11-03, exploratory):** All 30 leave-one-animal-out fits completed for genus, family and phylum (26 estimable targets). No target had p below 0.05 in the full data. The median of the largest single-deletion change in r was 0.089 (maximum 0.207: N01, genus, animal 9597). Seven targets changed sign under at least one deletion; all seven have a full-data |r| below 0.06. Five targets have |r| of 0.2 or more (N01, N10, N28, T4, C2): none changes sign, four stay within about 0.1 of their full-data value, and N01 ranges from 0.05 to 0.38. No label was changed. Per-target table: results/loo/influence_charolais_mrna.tsv.

**Amendment 19 (2026-11-03, before any read or per-animal phenotype of a new cohort was downloaded):** Expansion to a multi-cohort study with methane as a co-primary outcome (decision with Dr. Goswami: option B). Only run-level metadata from ENA has been inspected (metadata/candidates/).
(1) Inclusion rule for a cohort: rumen sample; public raw reads; a public per-animal phenotype file located before analysis; design covariates documented in the paper; at least 25 animals; no animal shared with another included cohort. Overlap is checked by birth year, site and IDs and recorded in metadata/candidates/.
(2) Cohorts and roles. RFI: Nelore PRJNA987743 (52, shotgun, existing); Charolais PRJEB70895 (30, mRNA and 16S, existing); Alberta PRJNA393057 (Li et al. 2019 Microbiome 7:92; 706 runs, 16S, RFI on 572 per Table S1; born 2013 at Kinsella, so independent of the 2014-born PRJNA448333 steers). Methane: Wallace et al. 2019 PRJNA517480 (1,016 cows, 7 farms, phenotypes in Data S1); Difford et al. 2018 / Zhang et al. 2020 PRJEB28065 (750 cows, 6 herds, S4 Table and Dryad); Nelore RME (existing). Wallace RFI or feed conversion is a sensitivity analysis only, because intake was partly estimated with alkane markers. SRUC beef (PRJEB10338, PRJEB21624, PRJEB31266) enters only if a per-animal methane file and at least 250 metagenomes are located before analysis; ENA currently lists 291 metagenomes against 359 reported. Nelore 16S PRJNA638250 is the same 52 animals: a 16S against shotgun calibration, never an independent cohort. PRJNA448333 enters only if its phenotypes arrive.
(3) Processing, blind to phenotypes. Shotgun: unchanged (Kraken2 and Bracken, two-pass GTDB). 16S: the Amendment 17 pipeline per cohort and marker. Primers are identified from read starts and recorded before processing. Pairs are merged when the expected overlap after primer removal is at least 40 bp; otherwise forward reads only. SILVA 138.1, minBoot 50, and the Amendment 17 genus-name rule. Phenotype files are downloaded with the reads but joined only after every genus table of that cohort is written.
(4) Targets. RFI: the existing pre-registered sets, unchanged. Methane primary set, chosen because the literature reports each in more than one country: M1, log ratio of the Methanobrevibacter gottschalkii clade (M. gottschalkii, M. millerae, M. thaueri) to the M. ruminantium clade (M. ruminantium, M. olleyae), higher with higher methane; M2, Succinivibrionaceae, lower with higher methane; M3, total archaeal relative abundance, two-sided with no direction stated, because the literature is split. Methane secondary set, reported in a single region: Sharpea, Kandleria and Megasphaera, each lower with higher methane. A target the cohort cannot resolve (for example M1 without species-level archaea) is untestable there. The existing Nelore RME targets are tested in the new methane cohorts as a further stratum.
(5) Per-cohort model: ANCOM-BC2 as in Amendment 14b, with covariates taken from each paper's design before phenotypes are joined. Alberta: diet + breed + sex + collection date. Wallace: farm (which carries country, breed, diet and methane method) + parity + days in milk where available. Difford/Zhang: herd + parity + lactation stage where available. A covariate missing from the public file is dropped and the drop is recorded. Phenotype: RFI as a z-score; methane as each cohort's own primary published trait, with units recorded. Effects are converted to partial r as before.
(6) Pooling: REML random effects per target, separately for RFI and for methane. Assay (16S, shotgun, mRNA) enters a meta-regression as a sensitivity analysis only. BH correction within each set. Labels unchanged. Criterion 3 at k of 3 or more: every leave-one-cohort-out pooled estimate is concordant with q below 0.05. At k = 2 the existing rule applies. Wallace is also run leave-one-farm-out, because farm and methane method coincide.
(7) Pre-specified contrast: the share of methane primary targets labelled replicated is compared with the share of RFI targets labelled replicated, with exact 95% CIs. This is descriptive; no threshold is added.
(8) Level 3 transfer: a signature score per animal (signed mean CLR of the set's targets; nothing fitted), tested only in cohorts not yet analysed (Alberta for RFI; Wallace and Difford/Zhang for methane) with 10,000 permutations within diet, farm or herd cells. Nelore and Charolais receive a dry run with shuffled phenotypes only, to test the code.
(9) Exploratory, labelled as such: rerun the original published methods (relative abundance with LEfSe or Spearman, no covariates) in Nelore and Alberta to show how covariate omission produces associations. Release harmonised tables, the effect-size atlas, a power calculator and code with a Zenodo DOI. Power is recomputed once the final cohort set is known, before any fit.

**Amendment 20 (2026-11-03): provenance finding for the PRJNA448333 phenotypes; no change to any analysis rule.** A March 2026 file, /data/rumen448333/input/sample_metadata.csv, gives RFI, DMI, ADG and FCR for 35 of the 48 PRJNA448333 steers. All 35 RFI values equal, to four decimals, the RFI of the animal with the same ID number in Additional file 1 of Li et al. 2019 Microbiome 7:92 (the Alberta cohort, 709 cattle born in 2013 at Kinsella; file raw/phenotypes/Li2019_709_MOESM1.xlsx, md5 1f04d2e6a73b45d534bcfcec01d3b5aa). The PRJNA448333 steers were born in 2014 and sampled at Lacombe, so these phenotypes belong to different animals that share ID numbers. The breed counts (11, 9, 28 against 16 per breed) and the bulls and heifers (against steers only) in that file follow from the same mismatch. Consequences: (1) every earlier analysis of PRJNA448333 against RFI, DMI, ADG or FCR, including the March and August 2026 reanalyses and the CausalRFI pilot, is invalid and is not used; collection-date results that used no phenotype are unaffected; (2) PRJNA448333 phenotypes remain unavailable, and Amendment 19 item 2 stands; (3) the Alberta phenotype file was opened only to count exact RFI equalities by ID number; no Alberta phenotype value was displayed and nothing was joined to microbiome data, so Amendment 19 item 3 holds for the Alberta cohort.

**Amendment 21 (2026-11-03, before any Alberta read was processed): paired-end 16S rules.** Read starts in PRJNA393057 (three runs inspected, sequences only) show a 27F-type forward primer (GAGTTTGATCMTGGCTCAG, without its leading A) and a 519R-type reverse primer (CCGCGGCKGCTGGCAC), 2 x 300 bp, about 90 bp of expected overlap after primer removal, so pairs are merged under Amendment 19 item 3. Where the Amendment 17 rules are silent for paired reads: (1) cutadapt removes both primers (-e 0.12, -O 12, pairs lacking either primer discarded). (2) Truncation: R1 at the last position, capped at 280, where the pooled median Q over all samples is at least 30; R2 at the last position, capped at 260, where it is at least 25. If the expected overlap after truncation falls below 40 bp, R2 is lengthened one position at a time until it reaches 40. (3) maxEE 2 for each read, truncQ 2, maxN 0, PhiX removed; mergePairs with minOverlap 20 and maxMismatch 0; consensus chimera removal; SILVA 138.1 with minBoot 50 and the Amendment 17 genus-name rule. (4) Fallback: if the median sample keeps fewer than 50% of its filtered pairs after merging, the cohort is processed with forward reads only under Amendment 17 instead; the switch and its trigger value are reported. (5) The bacterial amplicon cannot measure archaea reliably: archaeal targets are untestable in this cohort unless an archaeal amplicon deposit is found before analysis. These rules apply to every later paired-end 16S cohort, with the primers recorded per cohort in the same way.

**Amendment 22 (2026-11-03, while Alberta DADA2 was running and before any Alberta phenotype was joined): sample exclusion rules for 16S cohorts.** (1) Wrong marker: a run whose primer-matched fraction after trimming is below 10% and whose most frequent read starts are a non-16S marker (for example fungal ITS or eukaryotic 18S primers) is excluded as a mis-deposited amplicon. In Alberta (PRJNA393057) this applies to SRR9056754 and SRR9056234, whose read starts carry ITS2 and eukaryotic 18S primer sequences and which kept 0.3% and under 10% of pairs. (2) Minimum depth: a sample with fewer than 1,000 non-chimeric reads after DADA2 is excluded from that cohort's genus table. (3) Both rules apply to every 16S cohort, use read data only, and the excluded runs are listed with the reason in results/16s/<project>/excluded_runs.tsv. Animals whose only run is excluded are untestable in that cohort. Primer-matched fractions from 10% to 50% (16 Alberta samples) are retained, subject to rule 2.

**Amendment 23 (2026-11-03, after the Alberta design cross-tabs and before any Alberta association was estimated): covariates that cannot be separated.** The Alberta cross-tabs (counts only; scripts/64_harmonize_alberta.py) show that the Amendment 19 item 5 covariates are confounded: all 46 bulls are on Diet 3, Diet 4 is heifers only and Diet 2 steers only, and the Diet 4 animals are exactly those sampled on 5 and 13 February 2014, so the model diet + breed + sex + collection date is not of full rank. Rule, for Alberta and any later cohort with the same problem: (1) covariates that are collinear with one another are replaced by a single factor, management group, formed from their observed combinations (Alberta: diet x sex x collection date); (2) a management group with fewer than 5 animals with the phenotype is excluded from that cohort's fit and listed; (3) the remaining covariates stay as separate terms (Alberta: breed), and the fitted design must be of full rank, checked before the fit; (4) the RFI effect is then estimated within management groups, so diet, sex and sampling-date differences cannot enter it. Alberta model: management group + breed + RFI_z.

**Amendment 24 (2026-11-03, before any Alberta association was viewed; corrects Amendment 20):** Li et al. 2019 Microbiome 7:92 states that its 709 cattle were born in 2013 at the Roy Berg Kinsella Research Ranch and were sampled by oro-gastric tube at 292.9 ± 0.6 days of age. Li et al. 2019 Microbiome 7:6 states that its 48 steers were selected from a Kinsella herd of 738 born in 2014 and sampled at slaughter at Lacombe, but ENA dates those samples to July to October 2014, which fits a 2013 birth better. The corresponding author has written that the lab has no 2013 cohort. Whether the 48 steers are a subset of the 709 is therefore unresolved. (1) Correction to Amendment 20: its conclusion that the two sets are different animals rested on the stated birth years and is withdrawn. Whether the March 2026 phenotype file joined the right animals is unresolved: all 41 PRJNA448333 ID numbers present in Table S1 fall in the matching breed block, but 22 of them are recorded there as heifers or bulls, against a steers-only description. The consequences of Amendment 20 stand: earlier PRJNA448333 phenotype analyses are not used, and PRJNA448333 phenotypes remain unavailable until the lab confirms the ID mapping. (2) Independence: because the Li-set targets were derived from PRJNA448333, the primary Alberta analysis excludes every animal whose ID equals a PRJNA448333 source_material_id followed by A (list in harmonized/alberta_overlap_448333.txt, built from ENA metadata). The full Alberta cohort is a sensitivity analysis, labelled alberta_all. (3) This rule replaces the independence statement for Alberta in Amendment 19 item 2. (4) Disclosure: Alberta fits on all 567 animals under Amendment 23 completed at 16:45 on 2026-11-03 before this amendment. Their output files were not opened. They are moved unread to results/alberta/sealed_20261103/ and are not used. Both modes are rerun from scratch under this amendment.

**Amendment 25 (2026-11-03, before any Alberta estimate was viewed): RFI pooling with Alberta.** (1) The RFI arm is pooled over Nelore, Charolais mRNA and Alberta, with the Amendment 14 and 19 rules unchanged; Alberta enters as its primary analysis (PRJNA448333 ID numbers excluded), and the same pooling with alberta_all is the sensitivity analysis. (2) The Fisher-z variance 1/(n − 3 − q) uses each cohort's own number of covariate columns q, which equals 1/(df − 1). For Nelore and Charolais (q = 2) this is identical to the value used so far; for Alberta (q = 7 in the primary analysis, 8 in the sensitivity analysis) it is the correct value, where a fixed q = 2 would understate the variance slightly. (3) Charolais 16S stays unpooled. Wallace RFI (sensitivity only) and PRJNA448333 (only if its phenotypes arrive) are added later under Amendment 19. (4) Results go to results/replication_rfi3/. The two-cohort classification in results/replication/ is kept unchanged as the pre-expansion result.

**Amendment 26 (2026-11-03): the PRJNA448333 steers are part of the Alberta cohort; no change to any analysis rule.** The corresponding author of both Li et al. 2019 papers confirmed by email (3 November 2026, 17:53) that the 48 PRJNA448333 animals are the same animals as in the PRJNA393057 cohort, selected from it for slaughter; PRJNA393057 is rumen fluid from live animals (amplicon sequencing) and PRJNA448333 is rumen digesta collected at slaughter. Consequences: (1) the Amendment 24 primary Alberta analysis, which excludes these animals, is the correct independent test of the Li-set targets, and its result (0 of 27 replicated; identical labels with and without them) stands; (2) PRJNA448333 cannot serve as an independent replication cohort; (3) whether the March 2026 phenotype file joined the right records depends on whether the SRA source_material_id values equal the Table S1 Animal IDs without the suffix A, which has been asked of the authors; until they answer, PRJNA448333 phenotypes stay unused and Amendment 20 item 1 is neither confirmed nor withdrawn; (4) any within-animal comparison of live rumen fluid and slaughter digesta would be a new, exploratory analysis needing its own amendment.

**Amendment 27 (2026-11-04): corrections to the wording of Amendments 24 and 26; no change to any analysis rule or result.** (1) Amendment 24 item (1) states that all 41 PRJNA448333 ID numbers present in Table S1 fall in the matching breed block. The breed assignment used for that check came from the March 2026 joined file, which was itself built by matching ID numbers to Table S1, so the agreement follows from the join and is not evidence that the join is correct. That sentence is withdrawn. The other finding in the same item stands: Table S1 records 22 of these 41 IDs as heifers or bulls, against a steers-only description of PRJNA448333. (2) Amendment 26 item (1) calls the primary Alberta analysis "the correct independent test of the Li-set targets". The corresponding author confirmed that the 48 animals belong to the PRJNA393057 cohort but has not confirmed which Table S1 IDs they are. The primary analysis excludes animals by ID-number equality (Amendment 24 item 2), which is unverified, and the 22 heifers or bulls among the matched IDs suggest that it may not identify the 48 steers. Some of the 48 may therefore remain in the primary fit, and some excluded animals may be unrelated. Item (1) is replaced by: the primary Alberta analysis applies the pre-registered exclusion rule, and its independence from the Li et al. discovery animals is not established until the ID mapping is confirmed. (3) The conclusion does not depend on this: the primary and alberta_all classifications are identical (Amendment 25 result), and at most 48 of the 529 animals in the primary fit could be discovery animals. (4) If the authors supply the mapping, the primary analysis is rerun with the confirmed exclusion under a new amendment; if they do not, the current rule stays and the uncertainty is reported.

**Amendment 28 (2026-11-04, before any Difford or Wallace association was estimated): processing of the Difford (PRJEB28065) and Wallace (PRJNA517480) amplicons.** Only run metadata, read starts, ID fields and herd by sequencing-batch counts have been viewed; no phenotype value has been displayed or joined to microbiome data. Difford. (1) Markers: each run is bacterial V1-V3 (27F AGAGTTTGATCCTGGCTCAG / 534R ATTACCGCGGCTGCTGG) or archaeal V4-V6 (S-D-Arch-0519-a-S-15 CAGCMGCCGCGGTAA / S-D-Arch-1041-a-A-18 GGCCATGCACCWCCTCTC), as published in Difford et al. 2018. The read-start census (scripts/70_marker_census.py) assigns 1,302 bacterial and 1,015 archaeal runs and agrees with the Bact or Arch tag of every run that carries one. (2) Animals: runs are linked to S4 Table cows through the ENA sample alias with the submission prefix (IC-# or NG-#) and the library tag (_lib#, _Bact, _Arch, _V#) removed. The suffixes _c and _J are part of the cow ID; 81 numbers occur with and without a suffix and are different cows. This links 1,676 runs, giving bacterial data for 750 of 750 cows and archaeal data for 735. The remaining 641 runs (cows not in S4, the a and b libraries, and the Standard_A and Standard_B samples) are not used. Map: harmonized/difford_runs.tsv (scripts/71_difford_map.py). (3) Processing: cutadapt with the marker's forward primer required on R1 and reverse primer on R2 (--discard-untrimmed); DADA2 paired-end under the Amendment 21 rules, with error models learned separately for each flowcell (24 flowcells; MiSeq and HiSeq runs are never pooled for error learning). Sequence tables from several runs of the same cow and marker are summed before chimera removal. Amendment 22 exclusions apply. Taxonomy and genus naming as for Alberta. (4) Archaeal estimand: the archaeal amplicon measures abundance within the archaea, not within the whole community. Difford archaeal targets are estimated from that table and reported separately; they are not pooled with whole-community estimates. Bacterial targets use the bacterial table. (5) Model: herd + parity (S4 Lactatation, factor) + days in milk (linear) + sequencing batch (S4 Sequencing_Batch, factor) + methane_z, where methane_z is methane_grams_per_day as a z-score. Sequencing batch is added because the platform and flowcell differ between batches. Sequencing_Instrument is dropped if it is determined by batch. The design must be of full rank, checked before the fit; if it is not, herd and batch are combined into a management group under Amendment 23. Wallace. (6) Markers are taken from the library_name field (Prokaryotes, Bacteria, Archaea, Protozoa, Fungi; rounds SR1 and SR2). Each label has one dominant read start in the census, with median lengths of 252, 256, 88, 226 and 190 bp respectively. (7) The primary table for all bacterial and archaeal targets is Prokaryotes (ProkA, V4, MiSeq), the paper's main prokaryote marker, which measures abundance relative to the whole prokaryote community as the shotgun and universal 16S cohorts do. Sensitivity analyses: the Bacteria marker for bacterial targets, and the Archaea marker (HiSeq, about 88 bp) for archaeal targets, flagged as short-read. Fungi and Protozoa are not processed, because no pre-registered target uses them. (8) Processing: the deposited reads are already merged, have no primers and carry a flat quality score, so quality-based truncation and expected-error filtering do not apply and denoising rests on abundance alone; this is reported as a limitation. DADA2 single-end with maxN = 0, a length window of the marker's census median plus or minus 10 bp, error models learned per marker and sequencing round, and denoising per run. SR1 and SR2 tables are summed per cow before chimera removal. Amendment 22 exclusions apply. Taxonomy and genus naming as for Alberta. (9) Model per Amendment 19 item 5 (farm + parity + days in milk where Data S1 provides them). Methane is z-scored within farm, because the measurement method differs by country (breath sampling in the UK, GreenFeed in Italy and Sweden, respiration chambers in Finland).

**Amendment 29 (2026-11-04, before any Difford association was estimated): read subsampling for Difford and an explicit archaeal filter.** (1) A test on 8 runs per marker showed that full-depth DADA2 (about 600,000 read pairs per run) would take several days. Each trimmed Difford run is therefore subsampled at random to 50,000 read pairs (all pairs if fewer), with mates kept together and the seed fixed by run accession (scripts/75_subsample_pairs.py), before filtering and denoising. Cows with several runs keep the sum of their subsampled runs; differences in depth are handled by the model as before. (2) In the archaeal table, ASVs that SILVA does not assign to the kingdom Archaea are removed before genus aggregation, which implements Amendment 28 item 4; the test found an off-target share of 0.4%. (3) The test outputs were used only to check run time, truncation, merge rate and off-target share; no association was estimated and no phenotype was loaded.

**Amendment 30 (2026-11-04, before any Difford association was estimated): forward-read sensitivity table for Difford bacteria.** The 8-run test gave a bacterial median merge rate of 0.541, close to the Amendment 21 fallback threshold of 0.5. With R1 truncated at 250 bp the expected overlap is 39 bp, so pairs from taxa with longer V1-V3 amplicons are more likely to fail to merge, which can bias the merged table. (1) The Amendment 21 rule still decides the primary Difford bacterial table. (2) A forward-read table (R1 only, truncation by the Amendment 21 item 4 rule) is also always produced for Difford bacteria, and every Difford bacterial target is estimated from it as a sensitivity analysis; a target whose label differs between the two tables is reported as such. (3) The archaeal table (test merge rate 0.98) is unaffected. (4) The processing script's minimum-genus check is set to at least one genus for the archaeal table, which holds only a few methanogen genera. Qualities are not binned (22 to 23 distinct quality symbols on HiSeq and MiSeq), so error models are learned as usual.

**Amendment 31 (2026-11-05, before any Difford association was estimated): computational changes after a memory failure.** (1) The full Difford run started on 4 November exhausted the workstation's memory about six minutes after both DADA2 jobs began (system journal: memory pressure from 14:07, allocation failures from 14:12) and froze it; nothing beyond trimming and subsampling was produced. The trimmed and subsampled files are complete (1,788 bacterial and 1,564 archaeal files) and are reused. (2) The median quality profile that sets truncation (Amendment 21) is estimated from 1,000 randomly sampled reads per run instead of 5,000, still over all runs, so that the profile is not held in memory for about 4.5 million reads at once. (3) Within each flowcell, runs are denoised and merged one at a time with the flowcell's error models; DADA2 treats samples independently (no pooling), so results are unchanged and only memory use differs. (4) Each marker runs as a separate job with a 50 GB memory cap and lowered priority.

**Amendment 32 (2026-11-05, before any methane association was estimated in Difford or Wallace): independence of methane targets, M1 measurement, and record corrections.** Sources checked for the methane targets: Shi et al. 2014 and Kittelmann et al. 2014 (New Zealand sheep), Danielsson et al. 2017 (Swedish dairy cows, Lovsta research herd, SLU), Wallace et al. 2015 (UK beef), Kamke et al. 2016 (New Zealand sheep), Wallace et al. 2019, Difford et al. 2018 (text and S1 Table) and Zhang et al. 2020. (1) Independence rule: a target is not counted as replication in a cohort whose own published analysis reported an association between that target and methane, or contributed to choosing the target. Such a test is labelled re-analysis, reported, and enters only a sensitivity pooling. Applied: M2 (Succinivibrionaceae) in Wallace is re-analysis, because Wallace et al. 2019 reported Succinivibrionaceae as linked to methane in this cohort and the evidence summary behind Amendment 19 cited it. M3 in Difford is untestable (Amendment 19 item 4, Amendment 28 item 4) and was in any case examined by Zhang et al. 2020. Difford et al. 2018 name Sporobacter, Sphaerochaeta and Bacteroidales as their methane-associated taxa; their S1 Table lists estimates for Methanobrevibacter, Succinivibrionaceae genera and Sharpea, which did not contribute to target choice, so Difford remains independent for its testable targets. The S1 Table was viewed before any Difford fit. The same rule is applied to the Nelore RME publications before any Nelore methane test of the new targets. (2) Possible overlap: Danielsson et al. 2017 studied cows at Lovsta; Wallace et al. 2019 does not name its Swedish farm (SE1, 100 Nordic Red). Unless records show the animals differ, M1 and M2 in Wallace are also estimated without SE1 as a sensitivity analysis. (3) Pooling: a target's replication label (Amendment 19 item 6) uses only its independent cohorts; a target with fewer than three is reported as independent evidence in the cohorts it has and cannot be labelled replicated. (4) M1 measurement: RIM-DB (Seedorf et al. 2014) was not downloadable, so Methanobrevibacter ASVs are assigned by identity to NCBI RefSeq 16S type-strain sequences of the M. gottschalkii clade (M. gottschalkii, M. millerae, M. thaueri) and the M. ruminantium clade (M. ruminantium, M. olleyae). An ASV is assigned to a clade only if its best match to that clade has at least one fewer mismatch than its best match to the other clade; otherwise it is unassigned. Before any fit, an in silico check on the reference sequences confirms that each cohort's amplified region separates the two clades. M1 is untestable in a cohort if the region does not separate them or if fewer than 70% of the cohort's Methanobrevibacter reads are assigned. Per animal, M1 = log(gottschalkii-clade reads + 0.5) - log(ruminantium-clade reads + 0.5), analysed by linear regression with the cohort's covariates and converted to partial r; animals with fewer than 100 assigned Methanobrevibacter reads are excluded from M1. (5) Corrections: Amendment 31 item 3 should read that results are unchanged apart from DADA2's usual run-to-run variation from random read sampling in error learning. The truncation profile is now computed from per-position quality histograms, which gives medians identical to the previous method (checked on 100 runs, 0 positions differ) without holding all sampled reads in memory; filtering uses 8 processes.

**Amendment 33 (2026-11-05, before any Wallace association was estimated): Wallace read qualities, Wallace archaeal filter, and the Difford M1 gate.** (1) Correction to Amendment 28 item 8: the Wallace reads do not carry a flat quality score. Thirty runs per marker (2,000 reads each) show Phred scores from 3 to 80, with 67% (Bacteria), 71% (Prokaryotes) and 95% (Archaea) of bases above Q41 and under 0.6% below Q20, consistent with quality values recomputed during read merging. The earlier statement rested on too few reads and is withdrawn. (2) Wallace processing therefore uses DADA2's standard error model, which estimates the empirical error rate for each reported quality value from the data, with error models learned per marker and sequencing round as before, and the expected-error filter used for all other cohorts (maxEE = 2, truncQ = 2) together with the length window of Amendment 28 item 8. Reads removed by the filter are reported. (3) In the Wallace Archaea table, ASVs that SILVA does not assign to the kingdom Archaea are removed before genus aggregation, as for Difford (Amendment 29 item 2). (4) Difford M1 gate (Amendment 32 item 4): 240 Methanobrevibacter ASVs (64.6% of archaeal reads); 100% of Methanobrevibacter reads are assigned to a clade, and 99.7% are at least 10 bases closer to one clade's type strains than to the other's. M1 is therefore testable in Difford. The gate for Wallace is run on the Prokaryotes table after processing.

**Amendment 34 (2026-11-05, before any Wallace phenotype was joined to microbiome data): Wallace phenotypes, and a disclosure.** (1) Source: Wallace et al. 2019 Data S1 (aav8391_Data_S1.xlsx, from PMC6609165; raw/phenotypes/Wallace2019_Data_S1.xlsx, md5 1d94d40e85e6b02b5343aac241caa5b6), sheet "Phenotype data", 1,016 cows. Sequenced cows are linked by farm code and Cow ID (sample alias NR_FI1_S981 = farm FI1, Cow ID 981); all 1,012 sequenced cows match on both, and Cow IDs are unique. (2) Methane, primary: CH4 in g/d (column "Ruminal CH4 gas production, g/d"), the trait recorded on every farm and the unit used for Difford, z-scored within farm because the measurement method differs by country (column "Measurement method"). Sensitivity: CH4 in g/kg DMI, labelled as depending on intake estimated with alkane markers in the UK and Italy. (3) Covariates: farm (factor; it carries country, breed, diet and methane method), parity from "Lactation number" coded 1, 2, or 3 or more, and days in milk (linear). Cows missing the trait or a covariate are excluded from that fit and counted. (4) RFI sensitivity (Amendment 19 item 2): "RFI NRC", defined in the Data key as dry matter intake minus the NRC (2001) predicted intake, z-scored within farm. It is not a regression-based RFI and rests on alkane-estimated intake in the UK and Italy; it is reported as sensitivity only. (5) Disclosure: a structure check on 5 November matched numeric cells against numeric cow IDs and therefore printed three raw values from each of several Data S1 columns (ear tag, breed, days in milk, lactation number, live weight, diet composition, digestibility, milk urea, CO2 intensity, plasma creatinine). No methane value (columns 68 to 71) and no RFI value (column 52) was displayed, nothing was joined to microbiome data, and the displayed values cannot show any association. Later checks print names and counts only.

**Amendment 35 (2026-11-05): correction of family assignment in 16S cohorts.** (1) Error: family-level targets in the 16S cohorts (Alberta, Charolais 16S) were aggregated with the GTDB genus-to-family table built for the shotgun cohorts (results/replication/genus_lineage.tsv). 102 of Alberta's 361 SILVA genera, holding 43% of reads, are absent from that table and were pooled as unassigned_family. Of the two family-level targets, C2 (Succinivibrionaceae) therefore covered 4 genera and 0.86% of Alberta reads instead of 6 genera and 5.46%; Succinivibrionaceae_UCG-001 and UCG-002 were missing. T4 (Methanomassiliicoccaceae) has no Alberta reads under either mapping. Genus-level targets, and the shotgun cohorts, are unaffected. (2) Correction: in 16S cohorts a genus's family is taken from that cohort's own SILVA assignment (scripts/84_silva_lineage.R); a genus whose ASVs fall in more than one SILVA family takes the family holding most of its reads. (3) The Alberta family-level fits (primary and alberta_all) and the RFI pooling (Amendment 25) are rerun. The earlier outputs are kept unchanged in results/alberta/superseded_gtdb_family/ and results/replication_rfi3/superseded_gtdb_family/ and are reported as superseded. The Charolais 16S family fit is rerun the same way. (4) The methane cohorts use the corrected rule from the start. (5) This correction is made after the Amendment 25 results were viewed. It arose from a mapping error found while preparing the methane analysis, not from those results, and both versions are reported.

**Amendment 35 result (2026-11-05):** with SILVA families, Alberta C2 (Succinivibrionaceae) covers 6 genera and 5.46% of reads; Alberta r = -0.037 (primary; alberta_all -0.027), previously 0.054 (0.027). Pooled C2 r = 0.184 (95% CI -0.104 to 0.444, q = 0.54), previously 0.205; label unchanged (not replicated). No other target changed in either mode. The RFI arm still has 0 of 27 testable targets replicated. Superseded outputs remain in results/alberta/superseded_gtdb_family/ and results/replication_rfi3/superseded_gtdb_family/.

**Amendments 28 to 33, processing results (2026-11-05):** Difford bacteria: 894 runs, 750 cows; truncation R1 251, R2 260 (expected overlap 40 bp); median merge rate 0.424, below the Amendment 21 threshold, so the forward-read table (truncation 220 bp) is the primary and only bacterial table, and the Amendment 30 merged-versus-forward comparison does not arise; 55,420 ASVs, 357 genera, median 77% of reads classified to genus; 1 cow excluded under Amendment 22. Difford archaea: 735 cows, truncation 277/260, merge rate 0.989, 565 ASVs, 7 genera, 7 cows excluded; non-archaeal reads 0.2%. Wallace Prokaryotes: 2,022 runs, 1,012 cows, 12,392 ASVs, 451 genera, 79% classified, 0 excluded. Wallace Archaea: 1,012 cows, 726 ASVs, 0.6% non-archaeal reads. Wallace Bacteria: 1,877 of 1,920 runs passed filtering; processing in progress. M1 gate passed in both Difford (Amendment 33) and Wallace Prokaryotes (45 Methanobrevibacter ASVs, 100% assigned, 96.6% with a clade gap of at least 10 bases). Families for all 16S cohorts come from each cohort's SILVA assignment (Amendment 35).

**Amendment 36 (2026-11-05, before any methane association was estimated): methane target rows and estimation details.** (1) The Amendment 19 methane targets are added to plan/target_keys.tsv as set Methane (M1 clade log-ratio; M2 Succinivibrionaceae, family; M3 total archaea) and set MethaneSecondary (Sharpea, Kandleria, Megasphaera, genus), each with its published direction; BH correction is applied within each set, as for the RFI sets. (2) In 16S cohorts, targets are matched by exact name to SILVA genera, and families are formed from each cohort's SILVA lineage (Amendment 35). (3) M3 is estimated in the Wallace Prokaryotes table as one feature formed by summing all genera that SILVA assigns to the kingdom Archaea, fitted by ANCOM-BC2 alongside the other features. (4) M1 is estimated by linear regression of log((gottschalkii-clade reads + 0.5)/(ruminantium-clade reads + 0.5)) on the cohort covariates and the methane z-score (Amendment 32 item 4), from the Difford archaeal table and the Wallace Prokaryotes table; cows with fewer than 100 clade-assigned reads are excluded. (5) Covariates: Difford, herd + parity (1, 2, 3 or more) + days in milk + sequencing batch, with Sequencing_Instrument dropped if it is determined by batch and the design checked for full rank (Amendment 28 item 5); Wallace, farm + parity + days in milk (Amendment 34). Methane: Difford, CH4 g/d as a z-score; Wallace, CH4 g/d z-scored within farm. (6) Dry run: before any real fit, every cohort script is run with methane permuted at random within herd or farm (fixed seed), writing only to results/methane_dryrun/; this tests the code and is reported as such. (7) The Nelore RME targets (set Nelore, phenotype RME_GC) are tested in the new cohorts as a further stratum, matched by name as in (2).

**Amendment 36 dry-run result (2026-11-05 to 06):** scripts/81_estimate_methane.R ran end to end in shuffle mode on difford_arch, difford_bact and wallace_prokaryotes (designs of full rank; Difford 728 and 749 cows, Wallace 1,009; M1 from 728 and 853 cows with at least 100 clade-assigned reads). Thirty within-farm or within-herd permutations per table showed that the ANCOM-BC2 test statistic W is calibrated for prevalent targets (null SD of W 0.74 to 1.07 for M2, Pseudobutyrivibrio, Lachnospira; 0.67 to 0.87 for M1 and M3) but inflated for sparse targets (null SD 1.4 to 2.2 for Coprococcus, Eubacterium and Blautia at 17 to 59% prevalence; 4.66 for Sharpea in Wallace at 12% prevalence, where 67% of null permutations gave |W| > 1.96). The permutation itself is sound: extreme values fall on different targets under different seeds.

**Amendment 37 (2026-11-06, before any real methane fit): empirical calibration of ANCOM-BC2 statistics.** (1) For every cohort table and target, the phenotype is permuted 200 times within farm, herd or management group (seeds 1001 to 1200), the model is refitted, and the null standard deviation of W is recorded. (2) The calibrated statistic is W_cal = W / max(1, SD_null): a target whose null SD exceeds 1 is deflated, and a target whose null SD is below 1 is left unchanged, so calibration can only reduce significance. (3) W_cal replaces W in the conversion to partial r, in pooling (Amendment 14 and 25 rules), and in all labels; the Fisher-z variance is unchanged. A permutation p-value, the share of null |W| at or above the observed |W|, is reported beside the parametric one. (4) The rule applies to every cohort, including the RFI arm; the Nelore, Charolais and Alberta fits receive shuffle modes (permutation within diet or management group), the RFI arm is re-pooled with calibrated statistics, and the earlier classification is kept and reported as the uncalibrated version. (5) Real methane fits run only after the 200-permutation nulls exist for that table. (6) Wallace Bacteria (1,012 cows, 40,349 ASVs, 539 genera, 76% classified) is the sensitivity table for bacterial targets; the Wallace Archaea marker (88 bp) is not used for any target because M1 needs clade-level assignment that the marker's unknown primers do not permit, and M3 needs a whole-community table.

**Amendment 37 result, methane arm (2026-11-06):** 200-permutation nulls for all four tables; real fits and pooling run once, after the nulls. Primary pooling over the independent cohorts (Difford and Wallace; Wallace M2 as re-analysis): M1 clade log-ratio r = 0.030 (95% CI -0.020 to 0.079; Difford 0.032, Wallace 0.028), q = 0.47; M2 Succinivibrionaceae (Difford only) r = 0.006, q = 0.88; M3 total archaea (Wallace only) r = 0.032, q = 0.47; Sharpea r = -0.034, q = 0.55; Kandleria absent and Megasphaera below prevalence in both cohorts; all Nelore RME targets q = 1. No methane target is replicated or reaches single-cohort support. Calibration factors ranged from 1.0 (M1, M2, M3) to 5.3 (Sharpea in Wallace, raw W -9.38, calibrated -1.78, permutation p 0.07). Wallace M2 re-analysis: calibrated W -3.94 (lower with higher methane, as Wallace et al. 2019 reported); the independent Difford estimate is null. Lachnospira (N42) in Wallace: calibrated W -4.50, permutation p < 0.005, direction opposite to the Nelore report; Difford null. Outputs in results/methane/ and results/methane_pooled/.

**Methane sensitivity results (2026-11-06):** Wallace M2 (re-analysis) is stable across all Wallace subsets and traits (raw W -4.12 full, -4.17 without SE1, -4.21 for CH4 per kg DMI, -3.55 on the Bacteria marker) and null in Difford (W 0.16). Lachnospira likewise (-7.30, -5.90, -6.89, -5.98; Difford 0.04 after calibration). M1 is null for CH4 g/d in both cohorts (Wallace 0.83, no SE1 0.58, Difford -0.65) but W = 3.91 for CH4 per kg DMI in Wallace, a pre-registered sensitivity trait available only in Wallace; recorded as exploratory. M3 on the Bacteria marker is not interpretable (the bacterial primer under-amplifies archaea) and is not used. Pooling with Wallace M2 included: r = -0.061 (95% CI -0.186 to 0.066), q = 0.35, not replicated.

**Amendment 37 result, RFI arm (2026-11-08).** Item 4 done. 200 within-design permutations per table (Nelore within diet; Charolais mRNA within diet x year; Alberta within management group), seeds 1001 to 1200, 12 tables, 2,400 fits, 0 failures (scripts/88_rfi_null.sh, 2026-11-06 6:00 pm to 2026-11-07 8:13 pm; results/rfi_null/). Calibrated pooling in scripts/89_pool_rfi_calibrated.R, outputs in results/replication_rfi3_calibrated/; the uncalibrated outputs in results/replication_rfi3/ are kept and are now the superseded version, the calibrated classification is the primary one for the paper.
Null SD of W by table (median over targets): nelore_species 0.86, nelore_genus 1.05, nelore_family 1.16, nelore_phylum 1.01, charolais_mrna_species 0.90, charolais_mrna_genus 0.81, charolais_mrna_family 0.99, charolais_mrna_phylum 0.75, alberta_genus 0.92, alberta_family 1.12; maximum 1.66 (Fibrobacter, Alberta). ANCOM-BC2 is close to nominal on these denser tables, unlike the sparse methane targets (Amendment 37 result, methane arm). 17 of 50 cohort-level W values were shrunk, by at most 40 percent.
Result: 0 of 76 labels change between uncalibrated and calibrated pooling; 0 of 27 testable targets replicate; primary and sensitivity (Alberta all 567) agree on every label; 10 of 27 not-replicated targets still have a CI that includes the published r = 0.3, unchanged. Largest movements: Methanomassiliicoccaceae (T4) pooled r -0.279 to -0.256, CI -0.457 to -0.031, q 0.09 to 0.16, driven by Charolais mRNA (p_perm 0.015), Nelore p_perm 0.19; Succinivibrionaceae (C2) r 0.184 to 0.157, CI -0.105 to 0.399, q 0.65; Anaeromyxobacter (N01) r 0.257 to 0.195. Largest single-cohort permutation signals: Ruminococcus in Nelore (p_perm 0.015) and Succinivibrionaceae in Nelore (p_perm 0.030), both absent in the other cohorts. T4 is the closest thing to a signal in the RFI arm and is reported as such, not as a replication.

**Amendment 38 (2026-11-08, committed before any of it was run): positive controls for the three 16S cohorts, and a de novo discovery scan.**
1. Gap. Section 7 requires a positive control per cohort before its result is read as informative. Alberta, Difford and Wallace have none on record; until they pass one, their RFI and methane nulls are formally uninformative. This amendment supplies them.
2. Controls, one per cohort table, on the same animals, genus feature space (10 percent prevalence) and model family as the real fit (scripts 65 and 81), with the blocking factor replaced by the control contrast and the phenotype term retained: Alberta, diet (4 levels), diet + breed + RFI_z; Difford bacterial and archaeal tables, herd (6 levels), herd + parity + dim + CH4_z, labelled herd-or-batch because herd is confounded with sequencing batch and instrument; Wallace prokaryotic table, farm (7 levels), farm + parity + dim + CH4_z, labelled farm-or-breed because each farm has one breed and one methane method. ANCOM-BC2 with the real fits' settings (group NULL, struc_zero FALSE, prv_cut 0.10). Pass: at least one genus with BH q < 0.05 for any level of the control term. PERMANOVA (Aitchison distance, marginal, 1,999 permutations) and betadisper on the control term are reported alongside, as for Nelore and Charolais. A cohort that fails keeps its results but they are reported as uninformative (Section 7). Script 90_positive_controls_16s.R; outputs results/positive_controls/16s/. No phenotype value is printed.
3. De novo scan, exploratory, never merged with the replication labels (Section 8). For each real fit that already wrote a full ANCOM-BC2 table (Nelore genus, Charolais mRNA genus, Alberta genus primary, Difford bacterial, Wallace prokaryotic; the methane tables for CH4_z, the others for RFI), the phenotype coefficient is read for every genus (target aggregates TGT_* excluded), BH is applied across all genera of that table, and the number of genera at q < 0.05 and at nominal p < 0.05 is reported. Cross-cohort overlap is reported only between tables sharing a taxonomy namespace (GTDB: Nelore and Charolais; SILVA: Alberta, Difford, Wallace), as the count of genera nominally significant with the same sign in both tables against the hypergeometric expectation, and sign concordance over all shared genera (Spearman of W). Empirical calibration: 20 within-design permutations per table with the full ANCOM-BC2 table kept (results/scan_null/), giving the false-positive rate at nominal p < 0.05 and the ratio of median |W| to 0.674 per table. Twenty is enough for a table-level rate, not for per-feature SDs, and is labelled so. Script 91_denovo_scan.R; the permutations reuse 21, 65 and 81 with KEEP_CSV=1.
4. Reading rule, fixed now. If a cohort has within-cohort signal (genera at q < 0.05 after calibration) and little or no cross-cohort overlap beyond expectation, the paper says so: associations exist within herds and do not transfer. If no cohort has within-cohort signal, the paper says the phenotypes are not resolvable at taxon level at these sample sizes, and the replication nulls are read in that light. Neither reading changes any replication label.

**Amendment 38 item 2 result (2026-11-08).** All four 16S positive controls pass (results/positive_controls/16s/, logs/positive_controls/). Alberta, diet: 46 of 97 genera at BH q < 0.05 (min q 1.2e-14), PERMANOVA R2 0.15, p 5e-4. Difford bacterial, herd-or-batch: 77 of 104 (min q 6.2e-20), R2 0.15, p 5e-4. Difford archaeal, herd-or-batch: 3 of 3 (min q 8.4e-29), R2 0.14, p 5e-4. Wallace prokaryotic, farm-or-breed: 48 of 104 (min q 9.1e-61), R2 0.215, p 5e-4. Group dispersion differs in all four (betadisper p 5e-4), reported alongside. Under Section 7 the Alberta RFI and the Difford and Wallace methane results are informative. The 16S feature space after the 10 percent prevalence filter is about 100 genera per table (3 for the archaeal table), the same space the real fits used.

**Amendment 38 item 3 result, first pass (2026-11-09).** Nulls: 20 full-table permutations per RFI genus table (seeds 2001 to 2020, results/scan_null/), 200 per methane table (results/methane_dryrun/). Statistic: W = lfc/se throughout, as in the target work; ANCOM-BC2 2.4.0 p-values are not a normal transform of W (checked: W -14.2 with p 2e-6, and p not monotone in |W|), so they are used only for the nominal counts and the null false-positive rate. Deviation from item 3 recorded: the table-level inflation factor is replaced by the pooled empirical null of W (all genera, all permutations), and by a per-genus null where 200 permutations exist, because the null is heavy-tailed, not merely wide.
ANCOM-BC2 under the null: false-positive rate at nominal p < 0.05 of 0.089 (Nelore), 0.058 (Charolais mRNA), 0.225 (Alberta), 0.105 (Difford bacterial), 0.258 (Wallace prokaryotic); null W SD 1.31, 1.40, 1.96, 1.31, 2.37; null |W| 99th percentile 3.6, 4.6, 6.8, 4.6, 8.5 (2.6 if normal); largest number of genera at BH q < 0.05 in a single permutation 678, 1, 24, 21, 45.
RFI, real data against the pooled empirical null: 0 genera at BH q < 0.05 in Nelore (2,869 genera, top VMEW01 W -14.2, q_emp 0.75), Charolais mRNA (2,658, top Alkaliflexus W -27.7, q_emp 0.35) and Alberta (91, top Ruminobacter W 8.7, q_emp 0.41). Nelore x Charolais overlap of empirical-null hits: 10 observed vs 4.4 expected (p 0.01), 5 same-sign, Spearman of W over 2,654 shared genera 0.026; the same-sign genera (Sulfurifustis, Rudanella, Fen-1039, JAGLBA01, Erythrobacter) are environmental taxa, read as misassignment rather than biology.
Methane: Wallace 28 genera at per-genus permutation q < 0.05 (6 against the pooled null; top Dialister lfc -0.61 per SD CH4, W -23.3, q_perm 0.022; FD2005, Sphaerochaeta, Elusimicrobium, NED5E9 lower and Lactobacillus higher with methane). Difford 0 by either null. Across the 80 shared genera the two cohorts' W are negatively correlated (Spearman -0.30, p 0.007). Cross-trait pairs (Alberta RFI vs the methane tables) show no overlap.
Reading (item 4): RFI takes the second reading, no within-cohort genus-level signal in any of three cohorts under a calibrated test. Methane takes the first, within-cohort signal in Wallace that does not transfer to Difford. Both exploratory; no replication label changes. Pending before the RFI statement is final: 200 permutations per RFI genus table for a per-genus null (seeds 2021 to 2200), because the pooled null is conservative (Wallace: 6 pooled vs 28 per-genus) and "no signal" must not rest on a conservative null.

**Amendment 39 (2026-11-09): the Amendment 24 overlap exclusion was based on a spurious ID match; Alberta primary becomes the all-animal analysis.**
1. Evidence. ENA sample metadata for PRJNA448333 (raw/metadata/PRJNA448333_ena_runs.tsv): the 48 steers were sampled at Lacombe on seven dates, 2014-07-18 to 2014-10-03, and the first digit of each source_material_id (101 to 712) corresponds one to one with the slaughter date (1 = 18 July, 7 = 3 October), with within-batch numbers 01 to 12. The IDs are slaughter-batch sample codes, not animal tags. Table S1 of Microbiome 7:92 uses tag numbers with an A suffix (for example 906A). The number match used in Amendment 24 (101 to 101A) therefore linked sample codes to unrelated animals, which is why it produced heifers and bulls.
2. Birth year. Slaughter at finishing weight in July to October 2014 means the 48 were born in 2013, not 2014 as the 7:6 Methods state; this agrees with Dr. Guan's email (2026-11-06) that the 48 are a subset of the 709-animal 2013 cohort. The 48 are therefore present in PRJNA393057 but unidentifiable from public records. Difford of sample type also noted: the 48 are rumen digesta at slaughter, the 709 are rumen fluid by oro-gastric tube at about 293 days.
3. Change. The all-animal Alberta analysis (tag _all, 567 animals, previously the sensitivity) becomes the primary for the Li-set targets in Alberta; the 35-exclusion analysis (529 animals) is kept as superseded. The overlap with the Li discovery animals is reported as a bounded limitation: at most 48 of 567 Alberta animals (8.5 percent) are discovery animals for the Li set; for the Canonical and Nelore sets there is no overlap issue. Under both 66 and 89 the two versions agreed on every label, so no classification changes. A mapping list from the Guan lab, if it comes, replaces the bound with an exact exclusion under a further amendment.
4. Disclosure. While checking the Table S1 ID format, breed, sex and diet of two Alberta animals were displayed (906A and 904A: ANG, Heifer, Diet_4). No RFI value was displayed.
