# Extraction dictionary

`data/extraction_BP.csv` holds one row per comparison per outcome. A trial
reporting both systolic and diastolic pressure contributes two rows; a trial with
two active arms contributes four.

| Column | Definition |
|---|---|
| `study_id` | Identifier from the screening index. Suffixes a and b denote separate comparisons from one trial |
| `author`, `year`, `doi` | Citation identifiers |
| `design` | `parallel` or `crossover`. Determines how the variance is computed |
| `corr_r` | Within-subject correlation for crossover trials, where the trial reports it. Empty otherwise; the pipeline then assumes 0.5 |
| `shared_ctrl_id` | Shared text for rows of one trial that compare several active arms against the same control group. Empty for single comparisons |
| `source` | Botanical or product source of the anthocyanins |
| `form` | `purified`, `extract`, `juice` or `food` |
| `anthocyanin_dose_mg_day` | Daily anthocyanin dose in mg, where the trial quantified it. Empty where only a volume or food mass was reported |
| `dose_quantified` | `y` or `n`. Recorded for every trial; the proportion is an outcome of the review |
| `duration_weeks`, `duration_band` | Duration in weeks, and the band: acute (<2), short (2–4), medium (4–12), long (>12) |
| `bp_status` | `normotensive`, `prehypertensive`, `hypertensive` or `mixed` |
| `population` | Brief description of the participants |
| `outcome` | `SBP` or `DBP` |
| `n_int`, `mean_int`, `sd_int` | Sample size, mean and standard deviation of the anthocyanin arm |
| `n_ctrl`, `mean_ctrl`, `sd_ctrl` | The same for the control arm |
| `rob2_overall` | Consensus RoB 2 judgement |
| `notes` | How a standard deviation was derived, and any other departure from direct extraction |

## Conventions

Post-intervention values were used throughout, except where a trial reported only
change from baseline; the type taken is recorded in `notes` for every comparison.

Where a trial reported standard errors, standard deviations were obtained as
SE × √n. Where confidence intervals were reported, standard deviations were
derived from the interval width. Every derivation is recorded in `notes`.

**Cells are left empty where a value could not be obtained or validly derived.**
Those comparisons are dropped by the pipeline and written to
`outputs/log_dropped_rows.csv`. No value has been imputed.

For crossover comparisons, `n_int` and `n_ctrl` hold the number of participants,
not the sum of observations across both periods.

## Comparisons recovered by derivation

Two trials reported an effect estimate without arm-level means. Both were recovered rather than excluded. They carry the tag `[EFFECT RECOVERED]`
in `notes` and were excluded in a sensitivity analysis. Standard deviations
derived from standard errors are ordinary extraction and are not tagged.

One reported a between-group difference with its confidence interval. The
standard error was obtained from the interval width and converted to a group
standard deviation; the means were coded as 0 and the negative of the difference,
so that their difference reproduces the reported effect.

The other reported an effect with an exact p value and no measure of dispersion.
The standard error was obtained from the corresponding t statistic. Its design is
recorded as `parallel` although the trial used a crossover design, because the
derived variance already reflects the paired design.

---

# Screening and agreement files

`data/screening_decisions.csv` holds the independent decisions of the three
reviewers for the 1,523 deduplicated records of the bibliometric arm.
`data/screening_decisions_meta.csv` holds the same for the 204 records of the
meta-analysis arm. Both include the exclusion reason recorded for excluded
records.

`data/interrater_agreement.csv` reports Fleiss' kappa and pairwise Cohen's kappa
for study selection in each arm, and for the risk-of-bias judgements. The
selection values can be recomputed from the two screening files;
`R/02_screening_kappa.R` does so.
