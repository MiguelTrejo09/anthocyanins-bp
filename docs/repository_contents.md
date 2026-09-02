# What is in this repository

## Scripts (`R/`)

Ten scripts, numbered in execution order. `run_all.R` detects which inputs are
present and runs what it can. `00_setup.R` holds the shared configuration:
paths, the graphics device, the visual identity, and the effect-size functions
used by every analysis.

Scripts 01 to 03 cover the bibliometric arm and require the database exports.
Scripts 04 to 08 run from `data/extraction_BP.csv` alone, so the meta-analytic
results are reproducible without access to commercial databases.

## Data (`data/`)

| File | Rows | Contents |
|---|---|---|
| `extraction_BP.csv` | 77 | One row per comparison per outcome for the 39 trials contributing data. Column definitions in `docs/extraction_dictionary.md` |
| `excluded_at_extraction.csv` | 40 | Eligible trials that could not contribute data, with the reason for each |
| `rob2_assessments.csv` | 80 | Consensus RoB 2 judgements by domain. Column D6 holds the period and carryover domain and reads `NA` for parallel-group trials |
| `screening_decisions.csv` | 1 523 | Independent decisions of three reviewers for the bibliometric corpus |
| `screening_decisions_meta.csv` | 204 | The same for the meta-analytic corpus |
| `interrater_agreement.csv` | 14 | Fleiss' and pairwise Cohen's kappa for both arms and for risk of bias |
| `registry/clinicaltrials_search.csv` | 11 | ClinicalTrials.gov export used to look for trials without a locatable publication |

## Documentation (`docs/`)

`search_strings.md` gives the verbatim queries for all five database searches and
the search date. `extraction_dictionary.md` defines every column of the
extraction table and the conventions applied. `prisma_counts.csv` holds the study
flow in machine-readable form. `data_provenance.md` states which stage of the
review each data file documents.

## Figures (`figures_final/`)

The figures as published, at 600 dpi. Regenerating them requires running the
pipeline; they are included here so that the published versions remain available
independently of software versions.

## Reproducing the meta-analysis

```r
source("R/04_meta_analysis.R")
source("R/05_rob_figures.R")
source("R/06_forest_plots.R")
source("R/08_literature_comparison.R")
```

Expected output: systolic mean difference −2.43 mmHg (95 % CI −3.56 to −1.30,
k = 40, I² = 58.8 %) and diastolic −1.55 mmHg (−2.54 to −0.56, k = 37,
I² = 64.3 %).
