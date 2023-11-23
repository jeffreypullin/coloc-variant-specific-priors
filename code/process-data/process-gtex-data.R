source(here::here("renv/activate.R"))

library(readr)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(forcats)

gtex_paths <- snakemake@input[["gtex_paths"]]
plot_data_path <- snakemake@output[["plot_data_path"]]

process_file_gtex <- function(path) {

  raw_data <- read_tsv(path,
                       col_select = c(gene_id, pval_nominal, tss_distance),
                       show_col_types = FALSE, progress = FALSE)
  raw_data |>
    filter(pval_nominal < 5 * 10^-8) |>
    group_by(gene_id) |>
    arrange(pval_nominal) |>
    filter(row_number() == 1) |>
    mutate(abs_tss_distance = abs(tss_distance)) |>
    filter(abs_tss_distance > 0) |>
    ungroup() |>
    pull(abs_tss_distance)
}

gtex_data <- tibble(path = gtex_paths) |>
  rowwise() |>
  mutate(tissue = gsub("_", " ", unlist(strsplit(basename(path), "[.]"))[[1]])) |>
  mutate(abs_tss_distance = list(process_file_gtex(path))) |>
  ungroup() |>
  unnest(cols = abs_tss_distance)

saveRDS(gtex_data, here::here(plot_data_path))
