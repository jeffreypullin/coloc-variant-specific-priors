
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(biomaRt)
  library(janitor)
  library(readr)
})

density_to_tibble <- function(dens) {
  tibble(x = dens$x, y = dens$y)
}

processed_eqtlgen_data_path <- snakemake@input[["processed_eqtlgen_data_path"]]
processed_onek1k_data_path <- snakemake@input[["processed_onek1k_data_path"]]

eqtlgen_density_path <- snakemake@output[["eqtlgen_density_path"]]
onek1k_r1_density_path <- snakemake@output[["onek1k_r1_density_path"]]
onek1k_r2_density_path <- snakemake@output[["onek1k_r2_density_path"]]
onek1k_r3_density_path <- snakemake@output[["onek1k_r3_density_path"]]

processed_eqtlgen_data <- read_rds(processed_eqtlgen_data_path)
processed_onek1k_data <- read_rds(processed_onek1k_data_path)

processed_onek1k_data <- processed_onek1k_data |>
  filter(cell_type == "CD4 Naive/Central memory T cell")

onek1k_r1_density <- processed_onek1k_data  |>
  filter(round == 1) |>
  filter(abs(tss_distance) <= 5e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 8) |>
  density_to_tibble()

onek1k_r2_density <- processed_onek1k_data |>
  filter(round == 2) |>
  filter(abs(tss_distance) <= 5e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 8) |>
  density_to_tibble()

onek1k_r3_density <- processed_onek1k_data |>
  filter(round >= 3) |>
  filter(abs(tss_distance) <= 5e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 8) |>
  density_to_tibble()

eqtlgen_density <- processed_eqtlgen_data |>
  filter(abs(tss_distance) <= 5e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 8) |>
  density_to_tibble()

write_rds(eqtlgen_density, eqtlgen_density_path)
write_rds(onek1k_r1_density, onek1k_r1_density_path)
write_rds(onek1k_r2_density, onek1k_r2_density_path)
write_rds(onek1k_r3_density, onek1k_r3_density_path)
