source(here::here("renv/activate.R"))

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(tidyr)

data_paths <- snakemake@input[["adipos_express_marginal_paths"]]
processed_data_path <- snakemake@output[["processed_data_path"]]
tss_data_path <- snakemake@output[["tss_data_path"]]

process_file_adipo_express_marginal <- function(path, tss_data) {

  raw_data <- read_tsv(path, col_select = c(gene, pos, pval),
                       show_col_types = FALSE, progress = FALSE)

  gene_snp_data <- raw_data |>
    group_by(gene) |>
    arrange(pval) |>
    filter(row_number() == 1) |>
    ungroup()

  gene_snp_data |>
    left_join(tss_data, by = join_by(gene == ensembl_gene_id)) |>
    # NOTE: Remove genes that do not have a corresponding TSS.
    filter(!is.na(tss)) |>
    mutate(tss_distance = pos - tss) |>
    pull(tss_distance)
}

hg19_tss_data <- readRDS(tss_data_path)

adipos_express_marginal_data <- tibble(path = data_paths) |>
  mutate(file = gsub("(.+).txt", "\\1", basename(path))) |>
  rowwise() |>
  mutate(tss_distance = list(process_file_adipo_express_marginal(path, hg19_tss_data))) |>
  ungroup() |>
  unnest(cols = tss_distance) |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_tss_distance = log10(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

saveRDS(adipos_express_marginal_data, processed_data_path)
