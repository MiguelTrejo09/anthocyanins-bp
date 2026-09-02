# =============================================================================
# 03_bibliometric.R — descriptive indicators, laws, networks, thematic map
#
# Input:  outputs/corpus_included.rds
# Output: outputs/main_indicators.csv, outputs/annual_production.csv,
#         outputs/trend_tests.csv, outputs/bradford_zones.csv,
#         outputs/thematic_clusters.csv, outputs/dose_terms.csv,
#         figures/Fig_annual_production.tiff, Fig_keyword_network.tiff,
#         Fig_thematic_map.tiff, Fig_country_network.tiff
# =============================================================================
source("R/00_setup.R")

M <- readRDS(file.path(DIR_OUT, "corpus_included.rds"))
class(M) <- c("bibliometrixDB", "data.frame")
message("Corpus: ", nrow(M), " documents")

M <- bibliometrix::metaTagExtraction(M, Field = "AU_CO",  sep = ";")
M <- bibliometrix::metaTagExtraction(M, Field = "AU1_CO", sep = ";")

res <- bibliometrix::biblioAnalysis(M, sep = ";")

# ---- Descriptive indicators -------------------------------------------------
# safe() returns NA for fields not present in the installed bibliometrix version.
safe <- function(expr) {
  v <- tryCatch(expr, error = function(e) NA)
  if (is.null(v) || length(v) == 0) return(NA_character_) else as.character(v[1])
}
add <- function(df, k, v) rbind(df, data.frame(indicator = k, value = safe(v)))
ind <- data.frame(indicator = character(), value = character())

n_docs <- safe(res$Articles); n_auth <- safe(res$nAuthors); n_app <- safe(res$Appearances)
num <- function(x) suppressWarnings(as.numeric(x))

ind <- add(ind, "Documents", n_docs)
ind <- add(ind, "Source journals", length(unique(M$SO)))
ind <- add(ind, "Timespan", paste(min(M$PY, na.rm = TRUE), max(M$PY, na.rm = TRUE), sep = "-"))
ind <- add(ind, "Authors", n_auth)
ind <- add(ind, "Author appearances", n_app)
ind <- add(ind, "Authors per document", round(num(n_auth) / num(n_docs), 2))
ind <- add(ind, "Co-authors per document", round(num(n_app) / num(n_docs), 2))
ind <- add(ind, "International co-authorship (%)", {
  cc <- res$CountryCollaboration
  round(100 * sum(cc$MCP, na.rm = TRUE) / sum(cc$SCP + cc$MCP, na.rm = TRUE), 1) })
ind <- add(ind, "Mean citations per document",
           round(mean(as_num(M$TC), na.rm = TRUE), 2))
ind <- add(ind, "Documents never cited", sum(as_num(M$TC) == 0, na.rm = TRUE))
ind <- add(ind, "Documents with Keywords Plus", sum(!is.na(M$ID) & M$ID != ""))
ind <- add(ind, "Documents with author keywords", sum(!is.na(M$DE) & M$DE != ""))
ind <- add(ind, "Document types", paste(names(table(M$DT)), collapse = "; "))
readr::write_csv(ind, file.path(DIR_OUT, "main_indicators.csv"))
print(ind)

# ---- Annual production and trend --------------------------------------------
prod <- M |> dplyr::count(PY, name = "n") |> dplyr::filter(!is.na(PY)) |>
  dplyr::arrange(PY) |> dplyr::mutate(cumulative = cumsum(n))
readr::write_csv(prod, file.path(DIR_OUT, "annual_production.csv"))

mk <- Kendall::MannKendall(prod$n)
fit <- lm(n ~ PY, data = prod)
y0 <- min(prod$PY); y1 <- max(prod$PY)
cagr <- ((prod$n[prod$PY == y1] / prod$n[prod$PY == y0])^(1 / (y1 - y0)) - 1) * 100
readr::write_csv(tibble::tibble(
  statistic = c("Mann-Kendall tau","Mann-Kendall p","Linear slope (docs/year)",
                "Slope 95% CI lower","Slope 95% CI upper","R squared","CAGR (%)"),
  value = round(c(as.numeric(mk$tau), as.numeric(mk$sl), coef(fit)[["PY"]],
                  confint(fit)["PY", 1], confint(fit)["PY", 2],
                  summary(fit)$r.squared, cagr), 4)),
  file.path(DIR_OUT, "trend_tests.csv"))

p1 <- ggplot(prod, aes(PY, n)) +
  geom_area(fill = PAL_SEQ[2], alpha = .45) +
  geom_line(colour = PAL_SEQ[4], linewidth = .9) +
  geom_point(colour = PAL_SEQ[5], size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = "grey80",
              linetype = "22", linewidth = .5, alpha = .22) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .10))) +
  labs(title = "Research attention on anthocyanins and blood pressure",
       subtitle = sprintf("Annual output, %d-%d (n = %d documents), three databases combined\nMann-Kendall tau = %.2f, p < 0.001 | CAGR %.1f%%",
                          y0, y1, sum(prod$n), as.numeric(mk$tau), cagr),
       x = NULL, y = "Publications per year",
       caption = "Web of Science, Scopus and PubMed; deduplicated and screened by three independent reviewers.") +
  theme_pub()
save_fig(p1, "Fig_annual_production.tiff", 8, 5)

# ---- Bradford's law ---------------------------------------------------------
# Source names are normalised so that spelling variants of one journal are not
# ranked separately.
M2 <- M; M2$SO <- stringr::str_squish(toupper(M2$SO))
class(M2) <- c("bibliometrixDB", "data.frame")
B <- bibliometrix::bradford(M2)
readr::write_csv(B$table, file.path(DIR_OUT, "bradford_zones.csv"))
zones <- B$table |> dplyr::group_by(Zone) |>
  dplyr::summarise(journals = dplyr::n(), documents = sum(Freq), .groups = "drop")
readr::write_csv(zones, file.path(DIR_OUT, "bradford_summary.csv"))
print(as.data.frame(zones))

# ---- Networks ---------------------------------------------------------------
trim_net <- function(NM, drop = INDEXING_TERMS, top = 55) {
  NM <- as.matrix(NM)
  keep <- !(tolower(rownames(NM)) %in% tolower(drop))
  NM <- NM[keep, keep, drop = FALSE]
  sel <- names(sort(rowSums(NM), decreasing = TRUE))[seq_len(min(top, nrow(NM)))]
  NM[sel, sel, drop = FALSE]
}
gg_network <- function(NM, title, subtitle, caption, top = 55, seed = 2026) {
  NM <- trim_net(NM, top = top)
  g <- igraph::graph_from_adjacency_matrix(NM, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  g <- igraph::delete_vertices(g, igraph::degree(g) == 0)
  cl <- igraph::cluster_louvain(g)
  igraph::V(g)$cluster  <- as.factor(igraph::membership(cl))
  igraph::V(g)$strength <- igraph::strength(g)
  set.seed(seed)
  ggraph(tidygraph::as_tbl_graph(g), layout = "fr") +
    geom_edge_link(aes(width = weight, alpha = weight), colour = "grey62",
                   show.legend = FALSE) +
    scale_edge_width(range = c(.12, 1.4)) + scale_edge_alpha(range = c(.08, .45)) +
    geom_node_point(aes(size = strength, fill = cluster), shape = 21,
                    colour = "white", stroke = .5) +
    scale_size(range = c(2.5, 13), guide = "none") +
    scale_fill_manual(values = grDevices::colorRampPalette(PAL_CAT)(
      length(levels(igraph::V(g)$cluster))), name = "Cluster") +
    geom_node_text(aes(label = name), repel = TRUE, max.overlaps = 22,
                   colour = INK, size = 3, seed = seed,
                   bg.colour = "white", bg.r = .12) +
    labs(title = title, subtitle = subtitle, caption = caption) +
    theme_net() + guides(fill = guide_legend(override.aes = list(size = 4)))
}

NM_kw <- bibliometrix::biblioNetwork(M, analysis = "co-occurrences",
                                     network = "keywords", sep = ";")
save_fig(gg_network(NM_kw,
  "Conceptual structure of the field",
  "Keyword co-occurrence; node size = strength, colour = Louvain cluster",
  "Database indexing descriptors removed before clustering."),
  "Fig_keyword_network.tiff", 9, 7.5)

nc <- tryCatch({
  NM_co <- bibliometrix::biblioNetwork(M, analysis = "collaboration",
                                       network = "countries", sep = ";")
  gg_network(NM_co, "International collaboration network",
             "Country co-authorship; node size = strength, colour = Louvain cluster",
             "Derived from author affiliation.", top = 35)
}, error = function(e) { message("Country network skipped: ", conditionMessage(e)); NULL })
if (!is.null(nc)) save_fig(nc, "Fig_country_network.tiff", 8.5, 7)

# ---- Thematic map -----------------------------------------------------------
M3 <- M
strip <- function(x) vapply(strsplit(x, ";"), function(v) {
  v <- stringr::str_squish(v)
  paste(v[!(tolower(v) %in% tolower(INDEXING_TERMS))], collapse = ";")
}, character(1))
M3$ID <- strip(M3$ID)
class(M3) <- c("bibliometrixDB", "data.frame")

TM <- bibliometrix::thematicMap(M3, field = "ID", n = 250, minfreq = 5,
                                ngrams = 1, stemming = FALSE, size = 0.5,
                                n.labels = 3, repel = TRUE)
save_fig(TM$map +
  labs(title = "Thematic map of the field",
       subtitle = "Centrality (relevance) against density (development) of keyword clusters",
       caption = "Upper right: motor themes; lower right: basic; upper left: niche; lower left: emerging or declining.\nDatabase indexing descriptors removed before clustering.",
       size = "Cluster size (log frequency)") + theme_pub(),
  "Fig_thematic_map.tiff", 8.5, 7)
readr::write_csv(TM$words, file.path(DIR_OUT, "thematic_clusters.csv"))

# Frequency of dose-characterisation terms.
dose_terms <- TM$words |>
  dplyr::mutate(w = tolower(Words)) |>
  dplyr::filter(stringr::str_detect(w, "dose|bioavail|standard|pharmacokinet|quantif")) |>
  dplyr::select(Words, Occurrences, Cluster_Label) |>
  dplyr::arrange(dplyr::desc(Occurrences))
readr::write_csv(dose_terms, file.path(DIR_OUT, "dose_terms.csv"))
message("\nDose-characterisation terms in the thematic map:")
print(as.data.frame(dose_terms))
