# =============================================================================
# 05_rob_figures.R — risk-of-bias figures
#
# Input:  data/rob2_assessments.csv
# Output: figures/Fig_rob_traffic_light_*.tiff, Fig_rob_summary_*.tiff,
#         outputs/rob_domain_counts.csv
#
# Parallel and crossover trials are plotted separately, since the crossover
# variant of RoB 2 includes an additional domain for period and carryover
# effects and the two designs therefore do not share a domain set.
# =============================================================================
source("R/00_setup.R")

rob <- readr::read_csv(file.path(DIR_DATA, "rob2_assessments.csv"),
                       show_col_types = FALSE) |> as.data.frame()
rob$Weight <- 1

par_d   <- rob[rob$D6 == "NA", setdiff(names(rob), "D6")]
cross_d <- rob[rob$D6 != "NA", ]

tl <- robvis::rob_traffic_light(data = par_d, tool = "ROB2", psize = 8)
save_fig(tl, "Fig_rob_traffic_light_parallel.tiff", 8, max(6, 0.18 * nrow(par_d) + 3))

sm <- robvis::rob_summary(data = par_d, tool = "ROB2", weighted = FALSE) +
  labs(title = sprintf("Risk of bias, parallel-group trials (n = %d)", nrow(par_d)))
save_fig(sm, "Fig_rob_summary_parallel.tiff", 8, 3.2)

if (nrow(cross_d) > 1) {
  names(cross_d) <- c("Study", paste0("D", 1:6), "Overall", "Weight")
  tlc <- robvis::rob_traffic_light(data = cross_d, tool = "Generic", psize = 8)
  save_fig(tlc, "Fig_rob_traffic_light_crossover.tiff", 8, 4.5)
  smc <- robvis::rob_summary(data = cross_d, tool = "Generic", weighted = FALSE) +
    labs(title = sprintf("Risk of bias, crossover trials (n = %d), including period and carryover",
                         nrow(cross_d)))
  save_fig(smc, "Fig_rob_summary_crossover.tiff", 8, 3.4)
}

counts <- do.call(rbind, lapply(c(paste0("D", 1:6), "Overall"), function(d) {
  if (!d %in% names(rob)) return(NULL)
  data.frame(domain = d, as.data.frame(table(rob[[d]])))
}))
names(counts) <- c("domain", "judgement", "n")
readr::write_csv(counts, file.path(DIR_OUT, "rob_domain_counts.csv"))
print(counts)
