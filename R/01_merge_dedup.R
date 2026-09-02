# =============================================================================
# 01_merge_dedup.R — bibliometric arm: import, merge, deduplicate
#
# Input:  data_raw/wos*.txt, data_raw/scopus*.csv, data_raw/pubmed*.txt
#         Database exports are not redistributed here; their licences do not
#         allow it. docs/search_strings.md contains the verbatim queries and
#         the search date, so they can be regenerated.
# Output: outputs/corpus_merged.rds, data/screening_index.csv,
#         outputs/log_duplicates.csv, outputs/prisma_identification.csv
# =============================================================================
source("R/00_setup.R")

import_db <- function(pattern, dbsource, format) {
  files <- list.files(DIR_RAW, pattern, full.names = TRUE)
  if (!length(files)) stop("No files matching ", pattern, " in ", DIR_RAW)
  M <- bibliometrix::convert2df(files, dbsource = dbsource, format = format)
  message(sprintf("  %-8s %5d records", dbsource, nrow(M)))
  M
}

message("Importing raw exports")
M_wos <- import_db("^wos.*\\.txt$",    "wos",    "plaintext")
M_sco <- import_db("^scopus.*\\.csv$", "scopus", "csv")
M_pub <- import_db("^pubmed.*\\.txt$", "pubmed", "pubmed")

n_wos <- nrow(M_wos); n_sco <- nrow(M_sco); n_pub <- nrow(M_pub)
n_total <- n_wos + n_sco + n_pub

M <- bibliometrix::mergeDbSources(M_wos, M_sco, M_pub, remove.duplicated = TRUE)
n_dedup <- nrow(M)
message(sprintf("Merged: %d records identified, %d after deduplication",
                n_total, n_dedup))

saveRDS(M, file.path(DIR_OUT, "corpus_merged.rds"))

readr::write_csv(
  tibble::tibble(
    stage = c("Records identified - Web of Science",
              "Records identified - Scopus",
              "Records identified - PubMed",
              "Total records identified",
              "Duplicates removed",
              "Records after deduplication"),
    n = c(n_wos, n_sco, n_pub, n_total, n_total - n_dedup, n_dedup)),
  file.path(DIR_OUT, "prisma_identification.csv"))

# Index handed to the three reviewers for independent screening.
idx <- tibble::tibble(
  record_id = seq_len(nrow(M)),
  author = sub(";.*", "", M$AU),
  year   = M$PY,
  title  = M$TI,
  doi    = tolower(trimws(M$DI)),
  source = M$SO,
  reviewer1 = "", reviewer2 = "", reviewer3 = "", exclusion_reason = "")
readr::write_csv(idx, file.path(DIR_DATA, "screening_index_TEMPLATE.csv"))

message("Screening template written. The completed file with the reviewers' ",
        "decisions is data/screening_decisions.csv")
