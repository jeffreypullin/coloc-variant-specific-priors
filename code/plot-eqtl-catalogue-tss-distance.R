source("/home/jp2045/coloc-estimated-eqtl-priors/renv/activate.R")

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
boxplot_by_dataset_path <- snakemake@output[["boxplot_by_dataset_path"]]

process_file_eqtl_catalogue <- function(path, mart) {

  raw_data <- read_tsv(path, col_select = c(gene_id, position, pvalue),
                       show_col_types = FALSE, progress = FALSE)

  gene_snp_data <- raw_data |>
    group_by(gene_id) |>
    arrange(pvalue) |>
    filter(row_number() == 1) |>
    ungroup()

  tss_data <- getBM(
    attributes = c("ensembl_gene_id", "transcription_start_site"),
    filters = "ensembl_gene_id",
    values = gene_snp_data$gene_id,
    mart = mart
  ) |>
    as_tibble() |>
    rename(tss = transcription_start_site) |>
    group_by(ensembl_gene_id) |>
    # NOTE: For genes with multiple TSSs we take the median of the listed TSSs.
    summarise(tss = median(tss))

  gene_snp_data |>
    left_join(tss_data, by = join_by(gene_id == ensembl_gene_id)) |>
    # NOTE: Remove genes that do not have a corresponding TSS.
    filter(!is.na(tss)) |>
    mutate(abs_tss_distance = abs(position - tss)) |>
    filter(abs_tss_distance > 0) |>
    pull(abs_tss_distance)
}

eqtl_catalogue_metadata <- read_tsv(metadata_path)

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

eqtl_catalogue_data <- tibble(path = data_paths) |>
  mutate(file = gsub("(.+).cc.tsv", "\\1", basename(path))) |>
  rowwise() |>
  mutate(abs_tss_distance = list(process_file_eqtl_catalogue(path, mart))) |>
  ungroup() |>
  unnest(cols = abs_tss_distance)

tss_boxplot_by_dataset <- eqtl_catalogue_data |>
  left_join(eqtl_catalogue_metadata, join_by(file == dataset_id)) |>
  ggplot(aes(study_label, abs_tss_distance, fill = study_label)) +
  geom_boxplot(outlier.alpha = 0.1) +
  coord_flip(ylim = c(0, 150000)) +
  labs(
    x = "Dataset",
    y = "Distance from the TSS"
  ) +
  theme_bw()
ggsave(boxplot_by_dataset_path, plot = tss_boxplot_by_dataset)

# eqtl_catalogue_data |>
#   left_join(eqtl_catalogue_metadata, join_by(file == dataset_id)) |>
#   mutate(y = 10^log10_abs_tss_distance) |>
#   summarise(y = median(y), .by = c(file, sample_size)) |>
#   ggplot(aes(sample_size, y)) +
#   geom_point(size = 2) +
#   labs(
#     x = "Study sample size",
#     y = "Median distance to TSS"
#   ) +
#   coord_cartesian(ylim = c(10000, 25000), xlim = c(0, 700)) +
#   theme_bw()
