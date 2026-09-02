# =============================================================================
# 08_literature_comparison.R — post hoc analyses
#
# Input:  data/extraction_BP.csv
# Output: outputs/literature_comparison_{SBP,DBP}.csv
#
# Post hoc analyses examining the handling of crossover trials and the
# classification of intervention formats. Neither was specified in the registered
# protocol and both are reported as post hoc in the manuscript.
# =============================================================================
source("R/00_setup.R")

dat <- prepare_extraction(log_dropped = FALSE)

row_for <- function(d, label, mode = "paired") {
  blank <- function(note, k = 0) tibble::tibble(
    analysis = label, k = k, MD = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
    p = NA_real_, I2 = NA_real_, method = NA_character_, note = note)
  if (is.null(d) || nrow(d) < 2) return(blank("k < 2, not estimable",
                                              if (is.null(d)) 0 else nrow(d)))
  es <- compute_es(d, mode = mode)
  ft <- fit_rma(es)
  if (is.null(ft$model)) return(blank("not estimable", nrow(es)))
  m <- ft$model
  tibble::tibble(analysis = label, k = m$k, MD = round(m$b[1], 3),
                 CI_lower = round(m$ci.lb, 3), CI_upper = round(m$ci.ub, 3),
                 p = signif(m$pval, 4), I2 = round(m$I2, 1), method = ft$method,
                 note = ifelse(m$ci.lb > 0 | m$ci.ub < 0, "significant", "not significant"))
}

for (oc in c("SBP", "DBP")) {
  d <- dat[dat$outcome == oc, ]
  res <- dplyr::bind_rows(
    row_for(d, "Primary analysis (crossover paired)"),
    row_for(d, "Crossover treated as independent groups", mode = "independent"),
    row_for(d[d$form == "purified", ], "Purified anthocyanins only"),
    row_for(d[d$form != "purified", ], "Anthocyanin-rich products only"),
    row_for(d[d$form %in% c("juice", "food"), ], "Whole foods and juices only"),
    row_for(d[d$form == "extract", ], "Standardised extracts only"),
    row_for(d[d$duration_weeks >= 2 | is.na(d$duration_weeks), ], "Excluding acute trials"),
    row_for(d[d$bp_status == "hypertensive", ], "Hypertensive participants only"),
    row_for(d[d$bp_status == "normotensive", ], "Normotensive participants only"))
  readr::write_csv(res, file.path(DIR_OUT, paste0("literature_comparison_", oc, ".csv")))
  message("\n", oc); print(as.data.frame(res))
}
