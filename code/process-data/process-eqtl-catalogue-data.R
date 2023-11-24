source(here::here("renv/activate.R"))

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(tidyr)

metadata_path <- snakemake@input[["eqtl_catalogue_metadata"]]
data_paths <- snakemake@input[["eqtl_catalogue_paths"]]
tss_data_path <- snakemake@input[["tss_data_path"]]
processed_data_path <- snakemake@output[["processed_data_path"]]

process_file_eqtl_catalogue <- function(path, tss_data) {

  raw_data <- read_tsv(path, col_select = c(gene_id, position, pvalue),
                       show_col_types = FALSE, progress = FALSE)

  gene_snp_data <- raw_data |>
    group_by(gene_id) |>
    arrange(pvalue) |>
    filter(row_number() == 1) |>
    ungroup()

  gene_snp_data |>
    left_join(tss_data, by = join_by(gene_id == ensembl_gene_id)) |>
    # NOTE: Remove genes that do not have a corresponding TSS.
    filter(!is.na(tss)) |>
    mutate(tss_distance = position - tss) |>
    pull(tss_distance)
}

hg38_tss_data <- readRDS(tss_data_path)
eqtl_catalogue_metadata <- read_tsv(metadata_path,
                                    show_col_types = FALSE, progress = FALSE)

eqtl_catalogue_data <- tibble(path = data_paths) |>
  mutate(file = gsub("(.+).cc.tsv", "\\1", basename(path))) |>
  rowwise() |>
  mutate(tss_distance = list(process_file_eqtl_catalogue(path, hg38_tss_data))) |>
  ungroup() |>
  unnest(cols = tss_distance) |>
  left_join(eqtl_catalogue_metadata, join_by(file == dataset_id)) |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_tss_distance = log10(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

saveRDS(eqtl_catalogue_data, processed_data_path)
