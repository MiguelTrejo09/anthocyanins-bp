# Anthocyanins and blood pressure in adults

Data and code accompanying a systematic review and meta-analysis of randomized
controlled trials, reported together with a bibliometric analysis of the wider
literature.

**Manuscript:** *[title, journal, DOI once published]*
**Protocol:** https://doi.org/10.5281/zenodo.21829758
**PROSPERO:** CRD420261472107

---

## What this repository contains

| Path | Contents |
|---|---|
| `R/` | Analysis scripts, numbered in the order they run |
| `data/` | Extraction table, screening decisions, risk-of-bias assessments, agreement statistics |
| `docs/` | Verbatim search strings, extraction dictionary, PRISMA counts, data provenance |
| `figures_final/` | Figures as published, 600 dpi TIFF |
| `outputs/` | Written by the scripts when the pipeline is run |
| `figures/` | Written by the scripts when the pipeline is run |

## Reproducing the analysis

```r
# working directory: repository root
source("R/run_all.R")
```

Missing packages are installed by `R/00_setup.R` on first run. The versions used
for the published analysis are written to `outputs/sessionInfo.txt` at the end of
the pipeline.

Scripts 01 to 03 cover the bibliometric arm and require the raw database exports
in `data_raw/` (see below). Scripts 04 to 08 run from `data/extraction_BP.csv`
alone, so the meta-analytic results can be reproduced without them. `run_all.R`
detects which inputs are present and skips what it cannot run.

| Script | Purpose |
|---|---|
| `00_setup.R` | Paths, palette, graphics device, shared effect-size functions |
| `01_merge_dedup.R` | Import three databases, merge, deduplicate |
| `02_screening_kappa.R` | Inter-rater agreement for both arms, consensus corpora |
| `03_bibliometric.R` | Indicators, Bradford zones, networks, thematic map |
| `04_meta_analysis.R` | Pooled effects, heterogeneity, subgroups, sensitivity analyses |
| `05_rob_figures.R` | Risk-of-bias traffic light and summary plots |
| `06_forest_plots.R` | Forest and funnel plots |
| `07_decoupling_figure.R` | Attention against cumulative evidence |
| `08_literature_comparison.R` | Post hoc analyses |

## Analytical choices

**Crossover trials** are analysed as paired comparisons. The variance of the mean
difference is computed as (SD²int + SD²ctrl − 2·r·SDint·SDctrl) / n, with the
within-subject correlation taken from the trial where reported and set at 0.5
otherwise; sensitivity analyses run at 0.3 and 0.7. Treating these comparisons as
independent groups would overstate their variance.

**Multi-arm trials** in which two or more active arms share one control group have
that control's sample size divided across the arms, so control participants are
not counted more than once. Trials providing a separate control for each active
arm are treated as independent comparisons.

**Missing values** are left empty. Where a mean or dispersion measure could not be
obtained or validly derived from the source article, the comparison is dropped by
the pipeline and written to a log. No value has been imputed.

**Derived effect estimates.** Two comparisons whose effect was recovered from a
reported between-group difference or an exact p value, rather than from arm-level
means, carry the tag `[EFFECT RECOVERED]` in the notes column and are excluded in
a sensitivity analysis.

## What is not included, and why

Full texts of the included trials are not redistributed; they are the copyright of
their publishers. Every trial is identified by DOI in `data/extraction_BP.csv` and
`data/excluded_at_extraction.csv`.

The raw Web of Science and Scopus exports are not redistributed either, since
their licences do not permit it. `docs/search_strings.md` gives the verbatim
queries and the search date, so they can be regenerated with institutional access.
Place them in `data_raw/` as `wos*.txt`, `scopus*.csv` and `pubmed*.txt`.

No individual participant data were sought or held.

`data/registry/clinicaltrials_search.csv` holds the ClinicalTrials.gov export used to look for registered trials without a locatable publication. The WHO International Clinical Trials Registry Platform, also specified in the protocol, could not be accessed.

## Post hoc analyses

The analyses in `08_literature_comparison.R` were not specified in the registered
protocol and are reported as post hoc in the manuscript: exclusion of trials
shorter than two weeks, and reanalysis under the intervention-format
classification used by earlier syntheses.

## Citation

See `CITATION.cff`. Please cite the manuscript; the archived release carries its
own DOI for citing this version of the code.

## Licence

Code is released under the MIT licence (`LICENSE`). Data files are released under
CC BY 4.0 (`LICENSE-data`).
