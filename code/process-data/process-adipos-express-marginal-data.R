source(here::here("renv/activate.R"))

library(readr)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(biomaRt)
library(glue)
library(tools)
library(forcats)

data_paths <- snakemake@input[["adipos_express_marginal_paths"]]
plot_data_path <- snakemake@output[["plot_data_path"]]

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
    mutate(abs_tss_distance = abs(pos - tss)) |>
    filter(abs_tss_distance > 0) |>
    pull(abs_tss_distance)
}

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", "https://feb2014.archive.ensembl.org")

tss_data <- getBM(
  attributes = c("ensembl_gene_id", "transcript_start", "transcript_end", "strand"),
  filters = "ensembl_gene_id",
  values = "",
  mart = mart
) |>
  as_tibble() |>
  mutate(transcription_start_site = if_else(strand == 1, transcript_start, transcript_end)) |>
  dplyr::select(transcription_start_site, ensembl_gene_id) |>
  rename(tss = transcription_start_site) |>
  group_by(ensembl_gene_id) |>
  # NOTE: For genes with multiple TSSs we take the median of the listed TSSs.
  summarise(tss = median(tss))

adipos_express_marginal_data <- tibble(path = data_paths) |>
  mutate(file = gsub("(.+).txt", "\\1", basename(path))) |>
  rowwise() |>
  mutate(abs_tss_distance = list(process_file_adipo_express_marginal(path, tss_data))) |>
  ungroup() |>
  unnest(cols = abs_tss_distance)

saveRDS(adipos_express_marginal_data, plot_data_path)
