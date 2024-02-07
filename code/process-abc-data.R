source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(janitor)
  library(dplyr)
})

input_path <- snakemake@input[[1]]
output_path <- snakemake@output[[1]]

raw_abc_score_data <- read_tsv(
  input_path,
  name_repair = make_clean_names,
  show_col_types = FALSE
)

abc_score_data <- raw_abc_score_data |>
  select(chr, start, end, abc_score, cell_type, target_gene)

write_tsv(abc_score_data, output_path)
