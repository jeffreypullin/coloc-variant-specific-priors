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

# Create tsv files with only p-values that are < 5 * 10^-8 (genome-wide
# significant).
# setwd(here::here("data", "eqtl-catalogue"))
# dir.create("gws-sumstats")
#
# setwd(here::here("data", "eqtl-catalogue", "gws-sumstats"))
# list.files(here::here("data", "eqtl-catalogue", "sumstats"), full.names = TRUE) |>
#   walk(function(x) {
#     out_filename <- paste0("gws-", basename(tools::file_path_sans_ext(x)))
#     glue("zcat {{x}} | awk 'NR == 1 { print } NR != 1 { if ($9 <= 5E-8) { print } }' > {{out_filename}}",
#          .open = "{{", .close = "}}") |>
#       system()
#   })
