source(here::here("renv/activate.R"))

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(janitor)

data_path <- snakemake@input[["eqtlgen_path"]]
plot_data_path <- snakemake@output[["plot_data_path"]]

raw_eqtlgen_data <- read_tsv(data_path,
                             col_select = c(gene_symbol, snp_pos, gene_pos, pvalue),
                             name_repair = make_clean_names)

eqtlgen_data <- raw_eqtlgen_data |>
  filter(pvalue < 5 * 10^-8) |>
  group_by(gene_symbol) |>
  arrange(pvalue) |>
  filter(row_number() == 1) |>
  mutate(abs_tss_distance = abs(snp_pos - gene_pos)) |>
  filter(abs_tss_distance > 0) |>
  ungroup()

saveRDS(eqtlgen_data, plot_data_path)
