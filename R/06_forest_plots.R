# =============================================================================
# 06_forest_plots.R — forest and funnel plots
#
# Input:  data/extraction_BP.csv
# Output: figures/Fig_forest_{SBP,DBP}.tiff, Fig_funnel_{SBP,DBP}.tiff
#
# Drawn with ggplot2 so that they share the visual style of the bibliometric
# figures. The estimates are those produced by script 04.
# =============================================================================
source("R/00_setup.R")

dat <- prepare_extraction(log_dropped = FALSE)
MID <- 2   # minimally important difference, mmHg

forest_plot <- function(outcome) {
  d  <- dat[dat$outcome == outcome, ]
  es <- compute_es(d)
  m  <- fit_rma(es)$model
  if (is.null(m)) return(NULL)

  df <- tibble::tibble(
    study  = paste(es$author, es$year),
    design = ifelse(es$design == "crossover", "Crossover", "Parallel"),
    yi = es$yi, lo = es$yi - Z95 * sqrt(es$vi), hi = es$yi + Z95 * sqrt(es$vi),
    weight = 1 / (es$vi + m$tau2)) |>
    dplyr::mutate(weight_pct = 100 * weight / sum(weight)) |>
    dplyr::arrange(yi) |>
    dplyr::mutate(study = factor(make.unique(study), levels = make.unique(study)))

  rng <- range(c(df$lo, df$hi), na.rm = TRUE); pad <- diff(rng) * 0.34
  label <- ifelse(outcome == "SBP", "systolic blood pressure", "diastolic blood pressure")

  ggplot(df, aes(x = yi, y = study)) +
    annotate("rect", xmin = -MID, xmax = MID, ymin = -Inf, ymax = Inf,
             fill = PAL_SEQ[1], alpha = .55) +
    geom_vline(xintercept = 0, colour = INK, linewidth = .4) +
    geom_vline(xintercept = c(m$ci.lb, m$ci.ub), colour = PAL_SEQ[4],
               linetype = "22", linewidth = .3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = MUTED, linewidth = .45) +
    geom_point(aes(size = weight_pct, fill = design), shape = 22,
               colour = "white", stroke = .35) +
    scale_size(range = c(1.6, 5.4), name = "Weight (%)") +
    scale_fill_manual(values = c(Parallel = PAL_CAT[1], Crossover = PAL_CAT[2]),
                      name = "Design") +
    geom_text(aes(x = rng[2] + pad * .95,
                  label = sprintf("%+.2f [%+.2f, %+.2f]", yi, lo, hi)),
              hjust = 1, size = 2.5, colour = INK) +
    annotate("point", x = m$b[1], y = 0.25, shape = 23, size = 4.6,
             fill = PAL_SEQ[5], colour = "white", stroke = .4) +
    annotate("segment", x = m$ci.lb, xend = m$ci.ub, y = 0.25, yend = 0.25,
             colour = PAL_SEQ[5], linewidth = .8) +
    annotate("text", x = rng[2] + pad * .95, y = 0.25, hjust = 1, size = 2.8,
             fontface = "bold", colour = PAL_SEQ[5],
             label = sprintf("%+.2f [%+.2f, %+.2f]", m$b[1], m$ci.lb, m$ci.ub)) +
    annotate("text", x = rng[1] - pad * .1, y = 0.25, hjust = 0, size = 2.8,
             fontface = "bold", colour = PAL_SEQ[5],
             label = sprintf("Random-effects model (REML), k = %d", m$k)) +
    scale_x_continuous(limits = c(rng[1] - pad * .15, rng[2] + pad),
                       breaks = scales::pretty_breaks(7)) +
    scale_y_discrete(expand = expansion(add = c(1.4, .6))) +
    labs(title = sprintf("Effect of anthocyanins on %s", label),
         subtitle = sprintf("Pooled MD = %+.2f mmHg (95%% CI %+.2f to %+.2f) | I2 = %.0f%% | tau2 = %.2f",
                            m$b[1], m$ci.lb, m$ci.ub, m$I2, m$tau2),
         x = sprintf("Mean difference in %s (mmHg)   favours anthocyanins <-- | --> favours control", outcome),
         y = NULL,
         caption = sprintf("Crossover comparisons analysed as paired data (assumed within-subject r = %.1f). Shaded band marks the %d mmHg minimally important difference.\nSquare size is proportional to study weight.", R_ASSUMED, MID)) +
    theme_pub() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = rel(0.78)))
}

funnel_plot <- function(outcome) {
  d  <- dat[dat$outcome == outcome, ]
  es <- compute_es(d)
  m  <- fit_rma(es)$model
  if (is.null(m)) return(NULL)
  eg <- tryCatch(metafor::regtest(m), error = function(e) NULL)
  df <- tibble::tibble(yi = es$yi, se = sqrt(es$vi),
                       design = ifelse(es$design == "crossover", "Crossover", "Parallel"))
  ci <- tibble::tibble(se = seq(0, max(df$se, na.rm = TRUE) * 1.05, length.out = 200)) |>
    dplyr::mutate(lo = m$b[1] - Z95 * se, hi = m$b[1] + Z95 * se)

  ggplot() +
    geom_ribbon(data = ci, aes(y = se, xmin = lo, xmax = hi), fill = PAL_SEQ[1], alpha = .8) +
    geom_line(data = ci, aes(x = lo, y = se), colour = PAL_SEQ[3], linewidth = .35) +
    geom_line(data = ci, aes(x = hi, y = se), colour = PAL_SEQ[3], linewidth = .35) +
    geom_vline(xintercept = m$b[1], colour = PAL_SEQ[4], linetype = "22", linewidth = .5) +
    geom_point(data = df, aes(x = yi, y = se, fill = design), shape = 21,
               size = 2.6, colour = "white", stroke = .35) +
    scale_fill_manual(values = c(Parallel = PAL_CAT[1], Crossover = PAL_CAT[2]),
                      name = "Design") +
    scale_y_reverse(expand = expansion(mult = c(.04, .02))) +
    labs(title = sprintf("Funnel plot, %s", outcome),
         subtitle = if (!is.null(eg))
           sprintf("Egger test: z = %.2f, p = %.3f | k = %d", eg$zval, eg$pval, m$k)
         else sprintf("k = %d", m$k),
         x = sprintf("Mean difference in %s (mmHg)", outcome), y = "Standard error",
         caption = "Shaded region is the 95% pseudo-confidence area around the pooled estimate.") +
    theme_pub()
}

for (oc in c("SBP", "DBP")) {
  k <- sum(dat$outcome == oc)
  p <- forest_plot(oc)
  if (!is.null(p)) save_fig(p, paste0("Fig_forest_", oc, ".tiff"), 9, max(5.5, 0.22 * k + 3.2))
  f <- funnel_plot(oc)
  if (!is.null(f)) save_fig(f, paste0("Fig_funnel_", oc, ".tiff"), 6.8, 6.2)
  message("Figures for ", oc, " written.")
}
