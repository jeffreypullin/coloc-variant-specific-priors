source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
})

data_path <- snakemake@input[["eqtlgen_path"]]
hg19_tss_data_path <- snakemake@input[["hg19_tss_data_path"]]
processed_data_path <- snakemake@output[["processed_data_path"]]

hg19_tss_data <- read_rds(hg19_tss_data_path)

raw_eqtlgen_data <- read_tsv(
  data_path,
  col_select = c(gene_symbol, snp_pos, gene, gene_pos, pvalue),
  name_repair = make_clean_names,
  show_col_types = FALSE, progress = FALSE
)

eqtlgen_data <- raw_eqtlgen_data |>
  filter(pvalue < 5e-8) |>
  group_by(gene_symbol) |>
  arrange(pvalue) |>
  filter(row_number() == 1) |>
  ungroup() |>
  left_join(
    hg19_tss_data,
    by = join_by(gene == ensembl_gene_id),
    relationship = "many-to-many"
  ) |>
  mutate(tss_distance = if_else(strand == 1, snp_pos - tss, tss - snp_pos)) |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_tss_distance = log10(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

saveRDS(eqtlgen_data, processed_data_path)
