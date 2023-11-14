library(readr)
library(dplyr, warn.conflicts = FALSE)
library(glue)

eqtl_catalogue_metadata <- read_tsv(snakemake@input[["metadata_file"]],
                                    show_col_types = FALSE)
output_path <- snakemake@output[["dataset_file"]]

job_dataset_id <- gsub("\\..*", "", basename(output_path))

download_path <- eqtl_catalogue_metadata |>
  mutate(download_path = gsub("all", "cc", ftp_path)) |>
  filter(dataset_id == job_dataset_id) |>
  pull(download_path)

system(glue("wget -O {output_path} {download_path}"))
