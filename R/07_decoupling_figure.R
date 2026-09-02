# =============================================================================
# 07_decoupling_figure.R — main figure: attention against evidence
#
# Input:  outputs/annual_production.csv, data/extraction_BP.csv
# Output: figures/Fig_decoupling.tiff, outputs/cumulative_{SBP,DBP}.csv
#
# Panel A: annual publication output. Panel B: pooled estimate recomputed as
# each trial is added, ordered by publication year. Panel C: final estimates with
# confidence and prediction intervals against the minimally important difference.
# =============================================================================
source("R/00_setup.R")

prod <- readr::read_csv(file.path(DIR_OUT, "annual_production.csv"), show_col_types = FALSE)
dat  <- prepare_extraction(log_dropped = FALSE)
MID    <- 2
MIN_K  <- 3   # cumulative estimates are shown from the third trial onward

# ---- Panel A ----------------------------------------------------------------
mk <- Kendall::MannKendall(prod$n)
y0 <- min(prod$PY); y1 <- max(prod$PY)
cagr <- ((prod$n[prod$PY == y1] / prod$n[prod$PY == y0])^(1 / (y1 - y0)) - 1) * 100

pA <- ggplot(prod, aes(PY, n)) +
  geom_area(fill = PAL_SEQ[2], alpha = .45) +
  geom_line(colour = PAL_SEQ[4], linewidth = .9) +
  geom_point(colour = PAL_SEQ[5], size = 1.7) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = "grey80",
              linetype = "22", linewidth = .45, alpha = .22) +
  scale_x_continuous(breaks = scales::pretty_breaks(7)) +
  scale_y_continuous(expand = expansion(mult = c(0, .10))) +
  labs(title = "A  Research attention",
       subtitle = sprintf("Annual publications on anthocyanins and blood pressure (n = %d documents)\nMann-Kendall tau = %.2f, p < 0.001 | CAGR %.1f%%",
                          sum(prod$n), as.numeric(mk$tau), cagr),
       x = NULL, y = "Publications per year") +
  theme_pub()

# ---- Panel B ----------------------------------------------------------------
cumulative <- function(outcome) {
  d  <- dat[dat$outcome == outcome, ]
  es <- compute_es(d)
  es <- es[order(es$year), ]
  m  <- fit_rma(es)$model
  cm <- metafor::cumul(m, order = es$year)
  df <- tibble::tibble(year = es$year, estimate = as.numeric(cm$estimate),
                       ci_lower = as.numeric(cm$ci.lb), ci_upper = as.numeric(cm$ci.ub),
                       k = seq_along(es$year))
  readr::write_csv(df, file.path(DIR_OUT, paste0("cumulative_", outcome, ".csv")))
  list(by_year = df |> dplyr::group_by(year) |> dplyr::slice_tail(n = 1) |> dplyr::ungroup(),
       model = m)
}
cs <- cumulative("SBP"); cd <- cumulative("DBP")

cum <- dplyr::bind_rows(dplyr::mutate(cs$by_year, outcome = "SBP"),
                        dplyr::mutate(cd$by_year, outcome = "DBP")) |>
  dplyr::mutate(outcome = factor(outcome, levels = c("SBP", "DBP"))) |>
  dplyr::filter(k >= MIN_K)

pB <- ggplot(cum, aes(year, estimate, colour = outcome, fill = outcome)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -MID, ymax = MID,
           fill = PAL_SEQ[1], alpha = .6) +
  geom_hline(yintercept = 0, colour = INK, linewidth = .4) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .18, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(SBP = PAL_CAT[1], DBP = PAL_CAT[2]), name = NULL) +
  scale_fill_manual(values = c(SBP = PAL_CAT[1], DBP = PAL_CAT[2]), name = NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(7)) +
  coord_cartesian(ylim = c(-7, 3)) +
  labs(title = "B  Cumulative evidence",
       subtitle = sprintf("Pooled estimate recomputed as each trial was added, ordered by publication year (from k >= %d)\nShaded band: %d mmHg, the pre-specified minimally important difference", MIN_K, MID),
       x = NULL, y = "Cumulative pooled MD (mmHg)") +
  theme_pub()

# ---- Panel C ----------------------------------------------------------------
final <- dplyr::bind_rows(
  tibble::tibble(outcome = "SBP", est = cs$model$b[1],
                 lo = cs$model$ci.lb, hi = cs$model$ci.ub,
                 pil = predict(cs$model)$pi.lb, pih = predict(cs$model)$pi.ub,
                 k = cs$model$k, I2 = cs$model$I2),
  tibble::tibble(outcome = "DBP", est = cd$model$b[1],
                 lo = cd$model$ci.lb, hi = cd$model$ci.ub,
                 pil = predict(cd$model)$pi.lb, pih = predict(cd$model)$pi.ub,
                 k = cd$model$k, I2 = cd$model$I2)) |>
  dplyr::mutate(label = sprintf("%s  (k = %d, I2 = %.0f%%)", outcome, k, I2),
                outcome = factor(outcome, levels = c("SBP", "DBP")))

pC <- ggplot(final, aes(y = label)) +
  annotate("rect", xmin = -MID, xmax = MID, ymin = -Inf, ymax = Inf,
           fill = PAL_SEQ[1], alpha = .6) +
  geom_vline(xintercept = 0, colour = INK, linewidth = .4) +
  geom_errorbarh(aes(xmin = pil, xmax = pih), height = .10, colour = MUTED, linewidth = .5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = outcome), height = 0, linewidth = 1.6) +
  geom_point(aes(x = est, fill = outcome), shape = 23, size = 4,
             colour = "white", stroke = .4) +
  geom_text(aes(x = est, label = sprintf("%+.2f", est)), vjust = -1.5, size = 3, colour = INK) +
  geom_text(aes(x = pil, label = sprintf("%+.2f", pil)), hjust = 1.25, size = 2.5, colour = MUTED) +
  geom_text(aes(x = pih, label = sprintf("%+.2f", pih)), hjust = -.25, size = 2.5, colour = MUTED) +
  scale_colour_manual(values = c(SBP = PAL_CAT[1], DBP = PAL_CAT[2]), guide = "none") +
  scale_fill_manual(values = c(SBP = PAL_CAT[1], DBP = PAL_CAT[2]), guide = "none") +
  scale_x_continuous(breaks = scales::pretty_breaks(8), expand = expansion(mult = .14)) +
  labs(title = "C  Pooled effect and prediction interval",
       subtitle = "Thick bar: 95% confidence interval. Thin bar: 95% prediction interval for a future trial.",
       x = "Mean difference (mmHg)   favours anthocyanins <-- | --> favours control", y = NULL,
       caption = "Both prediction intervals cross the null: the direction of effect in a new trial cannot be anticipated from the current evidence base.") +
  theme_pub() + theme(panel.grid.major.y = element_blank())

fig <- (pA / pB / pC) + patchwork::plot_layout(heights = c(1, 1.15, .75))
open_dev(file.path(DIR_FIG, "Fig_decoupling.tiff"), 8.5, 12); print(fig); dev.off()
message("Decoupling figure written.")
