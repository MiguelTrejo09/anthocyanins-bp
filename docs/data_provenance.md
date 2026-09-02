# Data provenance

Each file in `data/` and the stage of the review it documents.

| File | Rows | Stage documented |
|---|---|---|
| `screening_decisions.csv` | 1 523 | Title, abstract and full-text screening of the bibliometric corpus by three independent reviewers, with the exclusion reason recorded for excluded records |
| `screening_decisions_meta.csv` | 204 | The same for the meta-analytic corpus |
| `interrater_agreement.csv` | 14 | Fleiss' and pairwise Cohen's kappa for study selection in both arms and for the risk-of-bias judgements |
| `rob2_assessments.csv` | 80 | Consensus RoB 2 judgements by domain for every eligible trial. Column D6 holds the period and carryover domain and reads `NA` for parallel-group trials |
| `extraction_BP.csv` | 77 | One row per comparison per outcome for the trials contributing data. Column definitions in `docs/extraction_dictionary.md` |
| `excluded_at_extraction.csv` | 40 | Eligible trials that could not contribute data, with the reason |

## Chain of numbers

The counts below are reproducible from the files above and are reported in the
PRISMA flow diagram. `docs/prisma_counts.csv` holds them in machine-readable form.

**Bibliometric arm.** 2 581 records identified across three databases; 1 058
duplicates removed; 1 523 screened; 1 288 excluded; 235 included.

**Meta-analytic arm.** 288 records identified; 84 duplicates removed; 204
screened; 115 excluded; 89 assessed at full text; 9 excluded; 80 eligible. Of
those, 40 excluded at the extraction stage, 1 with an active comparator analysed
separately, and 39 contributing 41 independent comparisons across 77 rows.

## Agreement

Study selection reached almost perfect agreement in both arms (Fleiss' κ = 0.983
and 0.973). Agreement on risk-of-bias judgements was poor (Fleiss' κ = 0.130),
reflecting systematic differences in the interpretation of two domains rather
than random variation; discrepancies were resolved by discussion against criteria
agreed in advance. Both the individual assessors' judgements and the consensus are
included so that the resolution can be traced.
