# enrichment analysis
library(dplyr)
library(ggplot2)
library(readr)

categories <- c("cDC1", "cDC2", "pDC")   
OR_CUT     <- 1.4
FDR_CUT    <- 0.05
TOP_N      <- 10                          

dir = "TRS/"
file.list <- list.files(path = dir, pattern='*.xlsx',full.names=T)
df.list <- lapply(file.list, read_excel)
attr(df.list, "names") <- file.list
TRS <- bind_rows(df.list, .id = "id")
counts <- TRS %>%
  dplyr::filter(category %in% categories) %>%
  dplyr::group_by(id, category) %>%
  dplyr::summarise(a         = sum(true_cell_top_idx, na.rm = TRUE),
                   n_lineage = dplyr::n(),
                   TRS.mean  = mean(TRS, na.rm = TRUE),
                   .groups   = "drop") %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(enriched_total = sum(a),
                total_cells    = sum(n_lineage)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(b = n_lineage - a,
                c = enriched_total - a,
                d = (total_cells - n_lineage) - (enriched_total - a)) %>%
  dplyr::filter(enriched_total > 0)

enr <- counts %>%
  dplyr::rowwise() %>%
  dplyr::mutate(p_value    = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
                odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)),
                log2_OR    = log2(odds_ratio),
                expected   = n_lineage * (enriched_total / total_cells),
                obs_exp    = a / expected) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(FDR         = stats::p.adjust(p_value, method = "BH"),
                significant = FDR < FDR_CUT & odds_ratio > OR_CUT) %>%
  dplyr::rename(enriched_cell = a, n_cells_category = n_lineage) %>%
  dplyr::arrange(dplyr::desc(significant), FDR, dplyr::desc(odds_ratio))

readr::write_csv(
  enr %>% dplyr::select(id, category, enriched_cell, n_cells_category,
                        enriched_total, total_cells, expected, obs_exp,
                        TRS.mean, odds_ratio, log2_OR, p_value, FDR, significant),
  op("dc_enrichment_results.csv"))
