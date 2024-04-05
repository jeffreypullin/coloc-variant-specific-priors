source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

input_path <- snakemake@input[[1]]
output_path <- snakemake@output[[1]]

raw_gnocchi_data <- read_tsv(
  input_path,
  show_col_types = FALSE
)

gnocchi_data <- raw_gnocchi_data |>
  setNames(c("chromosome", "start_pos", "end_pos", "score")) |>
  mutate(across(c(start_pos, end_pos, score), as.numeric))

write_tsv(gnocchi_data, output_path)
