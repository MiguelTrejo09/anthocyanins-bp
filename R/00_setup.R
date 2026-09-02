# =============================================================================
# 00_setup.R — shared configuration
#
# Sourced by every analysis script. Defines paths, the graphics device, the
# visual identity used across all figures, and the effect-size functions shared
# by the meta-analytic scripts.
#
# Paths are relative to the project root. Set the working directory to the
# repository root before running anything.
# =============================================================================

PKGS <- c("bibliometrix", "metafor", "dplyr", "tidyr", "readr", "readxl",
          "stringr", "tibble", "ggplot2", "ggraph", "tidygraph", "ggrepel",
          "patchwork", "scales", "Kendall", "igraph", "irr", "robvis", "ragg")
miss <- PKGS[!PKGS %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss)
invisible(lapply(PKGS, library, character.only = TRUE))
set.seed(2026)

DIR_DATA <- "data"
DIR_OUT  <- "outputs"
DIR_FIG  <- "figures"
DIR_RAW  <- "data_raw"
for (d in c(DIR_OUT, DIR_FIG)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- Visual identity --------------------------------------------------------
# Palette derived from anthocyanin pigment hues (cyanidin to malvidin).
PAL_SEQ <- c("#F2E3EF", "#D9A7C7", "#B5679C", "#7B3294", "#4A1D5E")
PAL_CAT <- c("#5B2C6F", "#1F7A8C", "#BF6B04", "#A63446", "#2E5E4E", "#7D6608")
INK <- "#22202A"; MUTED <- "#6E6A78"; GRID <- "#E6E4EA"

theme_pub <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text          = ggplot2::element_text(colour = INK),
      plot.title    = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2),
                                            hjust = 0, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(colour = MUTED, size = ggplot2::rel(0.92),
                                            hjust = 0, margin = ggplot2::margin(b = 11),
                                            lineheight = 1.15),
      plot.caption  = ggplot2::element_text(colour = MUTED, size = ggplot2::rel(0.76),
                                            hjust = 0, margin = ggplot2::margin(t = 11)),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.title    = ggplot2::element_text(colour = MUTED, size = ggplot2::rel(0.9)),
      axis.text     = ggplot2::element_text(colour = INK, size = ggplot2::rel(0.84)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
      plot.margin   = ggplot2::margin(15, 17, 11, 15),
      legend.position = "top", legend.justification = "left",
      legend.title  = ggplot2::element_text(size = ggplot2::rel(0.8), colour = MUTED),
      legend.text   = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.key.height = ggplot2::unit(3, "mm"),
      legend.key.width  = ggplot2::unit(6.5, "mm"))
}
theme_net <- function(base_size = 11) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = INK),
      plot.title    = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0),
      plot.subtitle = ggplot2::element_text(colour = MUTED, size = ggplot2::rel(0.92),
                                            hjust = 0, margin = ggplot2::margin(b = 8)),
      plot.caption  = ggplot2::element_text(colour = MUTED, size = ggplot2::rel(0.76), hjust = 0),
      plot.title.position = "plot", legend.position = "bottom",
      plot.margin = ggplot2::margin(15, 17, 11, 15))
}

# ---- Graphics device --------------------------------------------------------
open_dev <- function(file, width = 8, height = 6, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(file, width = width, height = height, units = "in",
                   res = dpi, compression = "lzw")
  } else {
    grDevices::tiff(file, width = width, height = height, units = "in", res = dpi,
                    type = if (capabilities("aqua")) "quartz" else "cairo")
  }
}
save_fig <- function(plot_obj, file, width = 8, height = 6, dpi = 600) {
  open_dev(file.path(DIR_FIG, file), width, height, dpi)
  on.exit(dev.off())
  print(plot_obj)
}

# ---- Shared helpers ---------------------------------------------------------
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Assumed within-subject correlation for crossover trials that do not report it.
R_ASSUMED <- 0.5

# Normal quantile for a 95% interval, used when drawing intervals around
# individual effect estimates.
Z95 <- qnorm(0.975)

# Database indexing descriptors carrying no conceptual content. Removed before
# keyword clustering; the list is reported in the supplementary material.
INDEXING_TERMS <- c(
  "human","humans","male","female","adult","adults","article","aged",
  "middle aged","very elderly","young adult","aged, 80 and over","controlled study",
  "major clinical study","clinical article","priority journal","nonhuman","animal",
  "animals","review","procedures","chemistry","unclassified drug","physiology",
  "comparative study","questionnaire","pathophysiology","normal human","blood",
  "blood sampling","drug blood level","drug urine level","age","age distribution",
  "ethnicity","united states","new zealand","medline","web of science",
  "search engine","data extraction","outcome assessment","follow up",
  "drug effect","drug effects")

prepare_extraction <- function(path = file.path(DIR_DATA, "extraction_BP.csv"),
                               log_dropped = TRUE) {
  d <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    sh <- readxl::excel_sheets(path)
    as.data.frame(readxl::read_excel(
      path, sheet = if ("extraction_BP" %in% sh) "extraction_BP" else sh[1]))
  } else as.data.frame(readr::read_csv(path, show_col_types = FALSE))

  for (v in c("n_int","mean_int","sd_int","n_ctrl","mean_ctrl","sd_ctrl",
              "corr_r","duration_weeks","anthocyanin_dose_mg_day","year"))
    if (v %in% names(d)) d[[v]] <- as_num(d[[v]])
  if (!"corr_r" %in% names(d)) d$corr_r <- NA_real_
  if (!"shared_ctrl_id" %in% names(d)) d$shared_ctrl_id <- NA_character_
  d$shared_ctrl_id[d$shared_ctrl_id %in% c("", "NA")] <- NA

  # Categorical fields are normalised on loading so that comparisons downstream
  # are not sensitive to case or stray whitespace.
  for (v in c("design", "form", "bp_status", "duration_band", "outcome",
              "dose_quantified", "rob2_overall")) {
    if (v %in% names(d)) d[[v]] <- tolower(trimws(d[[v]]))
  }
  d$outcome <- toupper(d$outcome)

  valid_design <- c("parallel", "crossover")
  if (any(!d$design %in% valid_design)) {
    stop("Unrecognised values in 'design': ",
         paste(unique(d$design[!d$design %in% valid_design]), collapse = ", "))
  }

  need <- c("n_int","mean_int","sd_int","n_ctrl","mean_ctrl","sd_ctrl")
  dropped <- dplyr::filter(d, dplyr::if_any(dplyr::all_of(need), is.na))
  if (nrow(dropped) && log_dropped) {
    readr::write_csv(dropped, file.path(DIR_OUT, "log_dropped_rows.csv"))
    message("Dropped ", nrow(dropped), " comparisons with missing values ",
            "(logged; never imputed).")
  }
  d <- dplyr::filter(d, dplyr::if_all(dplyr::all_of(need), ~ !is.na(.)))

  # Multi-arm correction: arms sharing one control split its sample size.
  d <- d |>
    dplyr::group_by(shared_ctrl_id, outcome) |>
    dplyr::mutate(n_arms = ifelse(is.na(dplyr::first(shared_ctrl_id)), 1L, dplyr::n()),
                  n_ctrl_adj = n_ctrl / n_arms) |>
    dplyr::ungroup() |> as.data.frame()

  shared <- dplyr::filter(d, n_arms > 1)
  if (nrow(shared)) {
    readr::write_csv(shared[, c("study_id","author","outcome","shared_ctrl_id",
                                "n_arms","n_ctrl","n_ctrl_adj")],
                     file.path(DIR_OUT, "log_shared_control.csv"))
    message("Shared-control correction applied to ", nrow(shared), " comparisons.")
  }
  d
}

# Parallel comparisons use the independent-groups variance. Crossover
# comparisons are paired: both conditions are measured in the same participants,
# so the variance incorporates the within-subject correlation.
# mode = "independent" reproduces the approach of earlier syntheses and is used
# in script 08.
compute_es <- function(d, r_assumed = R_ASSUMED, mode = c("paired","independent")) {
  mode <- match.arg(mode)
  is_cross <- d$design == "crossover" & mode == "paired"
  par <- d[!is_cross, , drop = FALSE]
  cro <- d[ is_cross, , drop = FALSE]
  out <- NULL
  if (nrow(par)) {
    e <- metafor::escalc(measure = "MD",
                         m1i = mean_int,  sd1i = sd_int,  n1i = n_int,
                         m2i = mean_ctrl, sd2i = sd_ctrl, n2i = n_ctrl_adj, data = par)
    par$yi <- as.numeric(e$yi); par$vi <- as.numeric(e$vi); out <- rbind(out, par)
  }
  if (nrow(cro)) {
    r  <- ifelse(is.na(cro$corr_r), r_assumed, cro$corr_r)
    np <- pmin(cro$n_int, cro$n_ctrl_adj)
    sd_diff2 <- pmax(cro$sd_int^2 + cro$sd_ctrl^2 - 2 * r * cro$sd_int * cro$sd_ctrl,
                     .Machine$double.eps)
    cro$yi <- cro$mean_int - cro$mean_ctrl
    cro$vi <- sd_diff2 / np
    out <- rbind(out, cro)
  }
  out
}

# Random-effects model. Where the between-study variance cannot be estimated by
# restricted maximum likelihood, a moment-based estimator is used instead; the
# estimator applied is returned with the model.
fit_rma <- function(es) {
  if (is.null(es) || nrow(es) < 2) return(list(model = NULL, method = NA_character_))
  m <- tryCatch(metafor::rma(yi, vi, data = es, method = "REML"),
                error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(m)) return(list(model = m, method = "REML"))
  m <- tryCatch(metafor::rma(yi, vi, data = es, method = "DL"), error = function(e) NULL)
  if (!is.null(m)) return(list(model = m, method = "DerSimonian-Laird"))
  m <- tryCatch(metafor::rma(yi, vi, data = es, method = "FE"), error = function(e) NULL)
  if (!is.null(m)) return(list(model = m, method = "fixed effect"))
  list(model = NULL, method = "not estimable")
}

message("Setup loaded. R ", getRversion(),
        " | metafor ", as.character(packageVersion("metafor")),
        " | bibliometrix ", as.character(packageVersion("bibliometrix")))
