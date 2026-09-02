# =============================================================================
# run_all.R — runs the full pipeline in order
#
# Scripts 01 to 03 require the database exports in data_raw/. Scripts 04 to 08
# run from data/extraction_BP.csv alone, so the meta-analytic results can be
# reproduced without them.
# =============================================================================
bibliometric <- dir.exists("data_raw") && length(list.files("data_raw"))

if (bibliometric) {
  source("R/01_merge_dedup.R")
  source("R/02_screening_kappa.R")
  source("R/03_bibliometric.R")
} else {
  message("data_raw/ is empty; skipping the bibliometric arm. ",
          "See docs/search_strings.md to regenerate the exports.")
}

source("R/04_meta_analysis.R")
source("R/05_rob_figures.R")
source("R/06_forest_plots.R")
if (bibliometric) source("R/07_decoupling_figure.R")
source("R/08_literature_comparison.R")

message("\nPipeline complete. Tables in outputs/, figures in figures/.")
writeLines(capture.output(sessionInfo()), file.path("outputs", "sessionInfo.txt"))
