source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
})

input_path <- snakemake@input[[1]]
output_path <- snakemake@output[[1]]

raw_gnocchi_data <- read_tsv(
  input_path,
  show_col_types = FALSE
)

colnames(raw_gnocchi_data) <- c("chromosome", "start_pos", "end_pos", "score")
gnocchi_data <- raw_gnocchi_data

write_tsv(gnocchi_data, output_path)
