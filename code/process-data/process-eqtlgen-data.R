source(here::here("renv/activate.R"))

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(janitor)

data_path <- snakemake@input[["eqtlgen_path"]]
processed_data_path <- snakemake@output[["plot_data_path"]]

raw_eqtlgen_data <- read_tsv(data_path,
                             col_select = c(gene_symbol, snp_pos, gene_pos, pvalue),
                             name_repair = make_clean_names,
                             show_col_types = FALSE, progress = FALSE)

eqtlgen_data <- raw_eqtlgen_data |>
  filter(pvalue < 5 * 10^-8) |>
  group_by(gene_symbol) |>
  arrange(pvalue) |>
  filter(row_number() == 1) |>
  mutate(tss_distance = snp_pos - gene_pos) |>
  ungroup() |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_tss_distance = log10(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

saveRDS(eqtlgen_data, processed_data_path)
