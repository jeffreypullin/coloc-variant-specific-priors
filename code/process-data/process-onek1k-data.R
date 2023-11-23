source(here::here("renv/activate.R"))

library(readr)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(biomaRt)
library(janitor)

onek1k_path <- snakemake@input[["onek1k_path"]]
plot_data_path <- snakemake@output[["plot_data_path"]]

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                   "https://feb2014.archive.ensembl.org")

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

raw_onek1k_data <- read_tsv(onek1k_path,
                            col_select = c(gene_id, pos, p_value, cell_type, round),
                            name_repair = make_clean_names)

onek1k_data <- raw_onek1k_data |>
  filter(p_value < 5 * 10^-8) |>
  group_by(gene_id, cell_type) |>
  arrange(p_value) |>
  filter(row_number() == 1) |>
  left_join(tss_data, by = join_by(gene_id == ensembl_gene_id)) |>
  mutate(abs_tss_distance = abs(pos - tss)) |>
  filter(abs_tss_distance > 0) |>
  ungroup() |>
  group_by(gene_id, abs_tss_distance, round) |>
  summarise(
    cell_type_collapsed = paste0(sort(unlist(cell_type)), collapse="-"),
    cell_type = list(cell_type),
    n_cell_types = n(),
    .groups = "drop"
  )

saveRDS(onek1k_data, plot_data_path)
