
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(janitor)
})

onek1k_path <- snakemake@input[["onek1k_path"]]
tss_data_path <- snakemake@input[["tss_data_path"]]
processed_data_path <- snakemake@output[["processed_data_path"]]

hg19_tss_data <- readRDS(tss_data_path)

raw_onek1k_data <- read_tsv(
  onek1k_path,
  col_select = c(gene_id, pos, p_value, cell_type, round),
  name_repair = make_clean_names,
  show_col_types = FALSE, progress = FALSE
)

onek1k_data <- raw_onek1k_data |>
  filter(p_value < 5e-8) |>
  group_by(gene_id, cell_type) |>
  arrange(p_value) |>
  filter(row_number() == 1) |>
  ungroup() |>
  left_join(hg19_tss_data, by = join_by(gene_id == ensembl_gene_id)) |>
  mutate(tss_distance = if_else(strand == 1, pos - tss, tss - pos))

saveRDS(onek1k_data, processed_data_path)
