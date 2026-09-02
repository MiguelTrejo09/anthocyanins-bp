# =============================================================================
# 02_screening_kappa.R — inter-rater agreement and consensus corpora
#
# Input:  data/screening_decisions.csv       (bibliometric arm, 1523 records)
#         data/screening_decisions_meta.csv  (meta-analysis arm, 204 records)
#         outputs/corpus_merged.rds
# Output: outputs/interrater_kappa.csv, outputs/prisma_screening.csv,
#         outputs/corpus_included.rds
#
# Agreement is computed from the three reviewers' independent decisions.
# Disagreements were resolved by discussion.
# =============================================================================
source("R/00_setup.R")

agreement <- function(path, arm) {
  d <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(dplyr::across(c(reviewer1, reviewer2, reviewer3),
                                ~ tolower(stringr::str_trim(.))))
  r <- as.data.frame(d[, c("reviewer1", "reviewer2", "reviewer3")])
  fk <- irr::kappam.fleiss(r)
  tbl <- tibble::tibble(
    arm = arm,
    metric = c("Fleiss kappa (3 raters)", "Cohen kappa R1-R2",
               "Cohen kappa R1-R3", "Cohen kappa R2-R3", "Observed agreement"),
    value = round(c(fk$value,
                    irr::kappa2(r[, c(1, 2)])$value,
                    irr::kappa2(r[, c(1, 3)])$value,
                    irr::kappa2(r[, c(2, 3)])$value,
                    mean(apply(r, 1, function(x) length(unique(x)) == 1))), 3),
    n_records = nrow(d))
  list(table = tbl, decisions = d)
}

bib  <- agreement(file.path(DIR_DATA, "screening_decisions.csv"), "bibliometric")
meta <- agreement(file.path(DIR_DATA, "screening_decisions_meta.csv"), "meta-analysis")

kappa_tbl <- dplyr::bind_rows(bib$table, meta$table)
readr::write_csv(kappa_tbl, file.path(DIR_OUT, "interrater_kappa.csv"))
print(as.data.frame(kappa_tbl))

consensus <- function(d) {
  d |>
    dplyr::rowwise() |>
    dplyr::mutate(n_include = sum(c(reviewer1, reviewer2, reviewer3) == "include"),
                  decision = ifelse(n_include >= 2, "include", "exclude")) |>
    dplyr::ungroup()
}
bib_d  <- consensus(bib$decisions)
meta_d <- consensus(meta$decisions)

readr::write_csv(tibble::tibble(
  arm = c(rep("bibliometric", 3), rep("meta-analysis", 3)),
  stage = rep(c("Records screened", "Records excluded", "Records included"), 2),
  n = c(nrow(bib_d), sum(bib_d$decision == "exclude"), sum(bib_d$decision == "include"),
        nrow(meta_d), sum(meta_d$decision == "exclude"), sum(meta_d$decision == "include"))),
  file.path(DIR_OUT, "prisma_screening.csv"))

M <- readRDS(file.path(DIR_OUT, "corpus_merged.rds"))
M$record_id <- seq_len(nrow(M))
M <- M[M$record_id %in% bib_d$record_id[bib_d$decision == "include"], ]
class(M) <- c("bibliometrixDB", "data.frame")
saveRDS(M, file.path(DIR_OUT, "corpus_included.rds"))

message("Bibliometric corpus: ", nrow(M), " documents | ",
        "Meta-analysis: ", sum(meta_d$decision == "include"), " records to full text")
