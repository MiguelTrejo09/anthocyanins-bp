# =============================================================================
# 04_meta_analysis.R — pooled effects, heterogeneity, subgroups, sensitivity
#
# Input:  data/extraction_BP.csv
# Output: outputs/meta_{SBP,DBP}_summary.csv, _bias.csv, _influence.csv,
#         outputs/meta_{SBP,DBP}_subgroup_*.csv,
#         outputs/meta_{SBP,DBP}_sens_*.csv,
#         outputs/pooled_estimates.csv
#
# Crossover comparisons are analysed as paired data and multi-arm trials share
# their control across arms; both corrections are applied in 00_setup.R.
# =============================================================================
source("R/00_setup.R")

dat <- prepare_extraction()
message("Comparisons available: ", nrow(dat))

run_meta <- function(d, outcome, r_assumed = R_ASSUMED, tag = "", figures = TRUE) {
  if (nrow(d) < 2) { message("  ", outcome, " ", tag, ": k < 2, skipped"); return(NULL) }
  es  <- compute_es(d, r_assumed)
  ft  <- fit_rma(es)
  if (is.null(ft$model)) { message("  ", outcome, " ", tag, ": not estimable"); return(NULL) }
  m   <- ft$model
  pr  <- predict(m)
  suffix <- if (tag == "") "_summary" else paste0("_sens_", tag)

  readr::write_csv(tibble::tibble(
    outcome = outcome, analysis = ifelse(tag == "", "primary", tag),
    method = ft$method, k = m$k,
    MD_mmHg = round(m$b[1], 3), CI_lower = round(m$ci.lb, 3), CI_upper = round(m$ci.ub, 3),
    p = signif(m$pval, 4), tau2 = round(m$tau2, 4), I2 = round(m$I2, 1),
    Q = round(m$QE, 2), Q_p = signif(m$QEp, 4),
    PI_lower = round(pr$pi.lb, 3), PI_upper = round(pr$pi.ub, 3),
    n_crossover = sum(es$design == "crossover"),
    n_parallel  = sum(es$design != "crossover")),
    file.path(DIR_OUT, paste0("meta_", outcome, suffix, ".csv")))

  if (figures && tag == "") {
    bias <- tryCatch({
      rt <- metafor::regtest(m); bt <- metafor::ranktest(m)
      tibble::tibble(outcome = outcome,
                     egger_z = round(rt$zval, 3), egger_p = signif(rt$pval, 4),
                     begg_tau = round(bt$tau, 3), begg_p = signif(bt$pval, 4),
                     note = ifelse(m$k < 10, "k < 10: tests underpowered", ""))
    }, error = function(e) tibble::tibble(outcome = outcome, note = conditionMessage(e)))
    readr::write_csv(bias, file.path(DIR_OUT, paste0("meta_", outcome, "_bias.csv")))

    inf <- tryCatch(influence(m), error = function(e) NULL)
    if (!is.null(inf)) {
      readr::write_csv(
        cbind(study = paste(es$author, es$year), as.data.frame(inf$inf)),
        file.path(DIR_OUT, paste0("meta_", outcome, "_influence.csv")))
    }
  }
  list(model = m, es = es)
}

message("\nPrimary analysis")
sbp <- run_meta(dat[dat$outcome == "SBP", ], "SBP")
dbp <- run_meta(dat[dat$outcome == "DBP", ], "DBP")
for (nm in c("SBP", "DBP")) {
  o <- get(tolower(nm))
  if (!is.null(o)) {
    m <- o$model
    message(sprintf("  %s: k = %d, MD = %+.2f mmHg (95%% CI %+.2f to %+.2f), p = %.4f, I2 = %.1f%%",
                    nm, m$k, m$b[1], m$ci.lb, m$ci.ub, m$pval, m$I2))
  }
}

# ---- Sensitivity: assumed within-subject correlation ------------------------
message("\nSensitivity: assumed correlation for crossover trials")
for (r in c(0.3, 0.5, 0.7)) {
  tg <- paste0("r", sub("\\.", "", format(r)))
  run_meta(dat[dat$outcome == "SBP", ], "SBP", r_assumed = r, tag = tg)
  run_meta(dat[dat$outcome == "DBP", ], "DBP", r_assumed = r, tag = tg)
}

# ---- Sensitivity: exclusions ------------------------------------------------
message("\nSensitivity: exclusions")
# Comparisons whose effect estimate was recovered from a reported difference or
# p value, rather than extracted from arm-level means, carry the tag
# [EFFECT RECOVERED] in the notes column and are excluded in a sensitivity
# analysis. Standard deviations derived from standard errors are ordinary
# extraction and are not flagged.
recovered <- unique(dat$study_id[grepl("[EFFECT RECOVERED]", dat$notes, fixed = TRUE)])
for (oc in c("SBP", "DBP")) {
  d <- dat[dat$outcome == oc, ]
  run_meta(d[!d$study_id %in% recovered, ], oc, tag = "excl_recovered")
  run_meta(d[d$rob2_overall != "high", ], oc, tag = "excl_high_rob")
  run_meta(d[d$form == "purified", ], oc, tag = "purified_only")
}

# ---- Pre-specified subgroup analyses ----------------------------------------
subgroup <- function(obj, moderator, outcome) {
  if (is.null(obj) || !moderator %in% names(obj$es)) return(invisible(NULL))
  es <- obj$es
  es[[moderator]] <- factor(es[[moderator]])
  if (nlevels(es[[moderator]]) < 2) return(invisible(NULL))
  m <- tryCatch(metafor::rma(yi, vi, mods = as.formula(paste("~", moderator)),
                             data = es, method = "REML"), error = function(e) NULL)
  if (is.null(m)) return(invisible(NULL))
  lv <- es |> dplyr::group_by(.data[[moderator]]) |>
    dplyr::summarise(k = dplyr::n(), .groups = "drop")
  readr::write_csv(tibble::tibble(
    outcome = outcome, moderator = moderator,
    QM = round(m$QM, 3), QM_p = signif(m$QMp, 4),
    R2 = ifelse(is.null(m$R2), NA, round(m$R2, 1)),
    levels = paste(lv[[1]], lv$k, sep = " = ", collapse = "; ")),
    file.path(DIR_OUT, paste0("meta_", outcome, "_subgroup_", moderator, ".csv")))
}
message("\nSubgroup analyses")
for (mod in c("form", "bp_status", "duration_band", "dose_quantified", "source")) {
  subgroup(sbp, mod, "SBP"); subgroup(dbp, mod, "DBP")
}

# ---- Estimates for the decoupling figure ------------------------------------
if (!is.null(sbp) && !is.null(dbp)) {
  readr::write_csv(tibble::tibble(
    outcome = c("SBP", "DBP"),
    MD = c(sbp$model$b[1], dbp$model$b[1]),
    CI_lower = c(sbp$model$ci.lb, dbp$model$ci.lb),
    CI_upper = c(sbp$model$ci.ub, dbp$model$ci.ub),
    PI_lower = c(predict(sbp$model)$pi.lb, predict(dbp$model)$pi.lb),
    PI_upper = c(predict(sbp$model)$pi.ub, predict(dbp$model)$pi.ub),
    k = c(sbp$model$k, dbp$model$k),
    I2 = c(sbp$model$I2, dbp$model$I2)),
    file.path(DIR_OUT, "pooled_estimates.csv"))
}
message("\nDone.")
