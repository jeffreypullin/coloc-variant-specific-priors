source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(glue)
  library(tools)
})

eqtl_catalogue_metadata <- read_tsv(snakemake@input[["metadata_file"]],
                                    show_col_types = FALSE)
file_output_path <- snakemake@output[["dataset_file"]]
tabix_output_path <- snakemake@output[["tabix_file"]]
permuted_output_path <- snakemake@output[["permuted_file"]]

job_dataset_id <- gsub("\\..*", "", basename(file_output_path))

file_download_path <- eqtl_catalogue_metadata |>
  mutate(download_path = gsub("all", "cc", ftp_path)) |>
  filter(dataset_id == job_dataset_id) |>
  pull(download_path)

tabix_download_path <- paste0(file_download_path, ".tbi")

permuted_download_path <- file_download_path |>
  file_path_sans_ext(compression = TRUE) |>
  file_path_sans_ext() |>
  paste0(".permuted.tsv.gz")

system(glue("wget -O {file_output_path} {file_download_path}"))
system(glue("wget -O {tabix_output_path} {tabix_download_path}"))
system(glue("wget -O {permuted_output_path} {permuted_download_path}"))
