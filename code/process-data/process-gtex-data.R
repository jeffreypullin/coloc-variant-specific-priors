source(here::here("renv/activate.R"))

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(tidyr)

gtex_paths <- snakemake@input[["gtex_paths"]]
processed_data_path <- snakemake@output[["processed_data_path"]]

process_file_gtex <- function(path) {

  raw_data <- read_tsv(path,
                       col_select = c(gene_id, pval_nominal, tss_distance),
                       show_col_types = FALSE, progress = FALSE)
  raw_data |>
    filter(pval_nominal < 5 * 10^-8) |>
    group_by(gene_id) |>
    arrange(pval_nominal) |>
    filter(row_number() == 1) |>
    ungroup() |>
    pull(tss_distance)
}

gtex_data <- tibble(path = gtex_paths) |>
  rowwise() |>
  mutate(tissue = gsub("_", " ", unlist(strsplit(basename(path), "[.]"))[[1]])) |>
  mutate(tss_distance = list(process_file_gtex(path))) |>
  ungroup() |>
  unnest(cols = tss_distance) |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_tss_distance = log10(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

saveRDS(gtex_data, processed_data_path)
