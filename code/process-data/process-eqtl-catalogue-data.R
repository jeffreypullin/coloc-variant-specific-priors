source(here::here("renv/activate.R"))

library(readr)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(biomaRt)
library(glue)
library(tools)
library(forcats)

metadata_path <- snakemake@input[["eqtl_catalogue_metadata"]]
data_paths <- snakemake@input[["eqtl_catalogue_paths"]]
plot_data_path <- snakemake@output[["plot_data_path"]]

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
    mutate(abs_tss_distance = abs(position - tss)) |>
    filter(abs_tss_distance > 0) |>
    pull(abs_tss_distance)
}

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

tss_data <- getBM(
  attributes = c("ensembl_gene_id", "transcription_start_site"),
  filters = "ensembl_gene_id",
  values = "",
  mart = mart
) |>
  as_tibble() |>
  rename(tss = transcription_start_site) |>
  group_by(ensembl_gene_id) |>
  # NOTE: For genes with multiple TSSs we take the median of the listed TSSs.
  summarise(tss = median(tss))

eqtl_catalogue_metadata <- read_tsv(metadata_path)

eqtl_catalogue_data <- tibble(path = data_paths) |>
  mutate(file = gsub("(.+).cc.tsv", "\\1", basename(path))) |>
  rowwise() |>
  mutate(abs_tss_distance = list(process_file_eqtl_catalogue(path, tss_data))) |>
  ungroup() |>
  unnest(cols = abs_tss_distance) |>
  left_join(eqtl_catalogue_metadata, join_by(file == dataset_id))

saveRDS(eqtl_catalogue_data, plot_data_path)
