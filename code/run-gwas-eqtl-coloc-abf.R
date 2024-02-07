# The code in this file is heavily based on:
# https://github.com/ralf-tambets/coloc/blob/main/bin/coloc_v3_pqtl.R

source(here::here("renv/activate.R"))

set.seed(26012024)

suppressPackageStartupMessages({
  library(readr)
  library(seqminer)
  library(dplyr)
  library(data.table)
  library(dtplyr)
  library(janitor)
  library(arrow)
  devtools::load_all("~/coloc")
})

source("code/coloc-utils.R")

eqtl_file <- snakemake@input[["eqtl_data_file"]]
gwas_file <- snakemake@input[["gwas_data_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])

width <- 1e5
coloc_metadata <- eqtl_metadata_file |>
  read_tsv(show_col_types = FALSE) |>
  filter(gene_type == "protein_coding") |>
  mutate(tss = if_else(strand == 1, gene_start, gene_end)) |>
  rowwise() |>
  mutate(
    start_pos = max(tss - width, 1),
    end_pos = tss + width
  ) |>
  mutate(region = paste0(chromosome, ":", start_pos, "-", end_pos)) |>
  select(gene_name, region, gene_id, chromosome, start_pos, end_pos, tss, gene_name) |>
  ungroup()

all_eqtl_data <- tabix.read.table(eqtl_file, paste0(chr, ":1-2147483647")) |>
  as_tibble()
all_gwas_data <- tabix.read.table(gwas_file, paste0(chr, ":1-2147483647")) |>
  as_tibble()

all_gwas_data <- all_gwas_data |>
  rename(
    rsid = SNPID,
    chromosome = CHR38,
    position = BP38,
    beta = BETA,
    se = SE,
    maf = ALT_FREQ
  ) |>
  mutate(
    an = 2 * sample_size,
    molecular_trait_id = "t1d",
    variant = paste0("chr", chromosome, "_", position, "_", REF, "_", ALT)
  ) |>
  select(-c(sample_size, REF, ALT)) |>
  lazy_dt()

coloc_metadata <- coloc_metadata |>
  filter(chromosome == chr) |>
  filter(!is.na(gene_id))

n <- nrow(coloc_metadata)

gnocchi_data <- read_tsv("data/gnocchi-windows.bed",
                         col_names = FALSE, show_col_types = FALSE)
colnames(gnocchi_data) <- c("chromosome", "start_pos", "end_pos", "score")

abc_score_data <- read_tsv("data/abc-score-data.txt.gz", show_col_types = FALSE)

density_data_round_1 <- read_rds("output/densities/onek1k_cd4nc_round_1.rds")
density_data_round_2 <- read_rds("output/densities/onek1k_cd4nc_round_2.rds")
density_data_round_3 <- read_rds("output/densities/onek1k_cd4nc_round_3.rds")
eqtlgen_density_data <- read_rds("output/densities/eqtlgen.rds")

snp_var_data_1_7 <- read_parquet("data/snpvar_meta.chr1_7.parquet")
snp_var_data_8_22 <- read_parquet("data/snpvar_meta.chr8_22.parquet")

results <- list()
for (i in 1:n) {

  region <- paste0(coloc_metadata$chromosome[[i]], ":",
                   coloc_metadata$start_pos[[i]], "-",
                   coloc_metadata$end_pos[[i]])

  print(paste0("Region: ", region))
  print(paste0("Gene ID: ", coloc_metadata$gene_id[[i]]))
  print(paste0("Gene name: ", coloc_metadata$gene_name[[i]]))

  eqtl_data <- all_eqtl_data |>
    setNames(eqtl_catalouge_colnames) |>
    filter_qtl_dataset(
      trait_id = coloc_metadata$gene_id[[i]],
      chrom = coloc_metadata$chromosome[[i]],
      start_pos = coloc_metadata$start_pos[[i]],
      end_pos = coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble()

  gwas_data <- all_gwas_data |>
    filter(chromosome == coloc_metadata$chromosome[[i]]) |>
    filter(
      position >= coloc_metadata$start_pos[[i]] &
      position <= coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble()

  if (min(gwas_data$P) > 5e-6) {
    next
  }

  if (nrow(eqtl_data) == 0 || nrow(gwas_data) == 0) {
    next
  }

  gwas_data <- prepare_coloc_dataset(gwas_data)
  eqtl_data <- prepare_coloc_dataset(eqtl_data)

  n_before <- max(nrow(eqtl_data), nrow(gwas_data))
  if (n_before <= 300) {
    next
  }

  eqtl_data <- eqtl_data |>
    filter(variant %in% gwas_data$variant)
  gwas_data <- gwas_data |>
    filter(variant %in% eqtl_data$variant)

  n_after <- max(nrow(eqtl_data), nrow(gwas_data))
  if (n_after / n_before < 0.1) {
    next
  }

  is_oneover_na <- all(is.na(1 / eqtl_data$se^2))
  is_nvx_na <- all(is.na(2 * eqtl_data$an / 2 * eqtl_data$maf * (1 - eqtl_data$maf)))

  if (is_oneover_na || is_nvx_na) {
    next
  }

  eqtl_dataset <- list(
    varbeta = eqtl_data$se^2,
    N = eqtl_data$an / 2,
    MAF = eqtl_data$maf,
    type = "quant",
    beta = eqtl_data$beta,
    snp = eqtl_data$variant,
    position = eqtl_data$position
  )

  is_oneover_na <- all(is.na(1 / gwas_data$se^2))
  is_nvx_na <- all(is.na(2 * gwas_data$an / 2 * gwas_data$maf * (1 - gwas_data$maf)))

  if (is_oneover_na || is_nvx_na) {
    next
  }

  gwas_dataset <- list(
    varbeta = gwas_data$se^2,
    N = gwas_data$an / 2,
    MAF = gwas_data$maf,
    type = "quant",
    beta = gwas_data$beta,
    snp = gwas_data$variant
  )

  eqtl_prior_weights_eqtlgen <- compute_eqtl_tss_dist_prior_weights(
    eqtl_dataset$position, coloc_metadata$tss[[i]], eqtlgen_density_data
  )
  eqtl_prior_weights <- compute_eqtl_tss_dist_prior_weights(
    eqtl_dataset$position, coloc_metadata$tss[[i]], density_data_round_1
  )
  eqtl_prior_weights_round_2 <- compute_eqtl_tss_dist_prior_weights(
    eqtl_dataset$position, coloc_metadata$tss[[i]], density_data_round_2
  )
  eqtl_prior_weights_round_3 <- compute_eqtl_tss_dist_prior_weights(
    eqtl_dataset$position, coloc_metadata$tss[[i]], density_data_round_3
  )
  gnocchi_prior_weights <- compute_gnocchi_prior_weights(
    eqtl_dataset$position, chr, gnocchi_data
  )
  abc_score_prior_weights <- compute_abc_prior_weights(
    eqtl_dataset$position, chr, coloc_metadata$gene_name[[i]], abc_score_data
  )
  abc_score_primary_blood_prior_weights <- compute_abc_prior_weights(
    eqtl_dataset$position, chr, coloc_metadata$gene_name[[i]], abc_score_data, 
    biosamples = "primary_blood"
  )

  if (chr %in% 1:7) {
    polyfun_data <- snp_var_data_1_7
  } else {
    polyfun_data <- snp_var_data_8_22
  }

  polyfun_prior_weights <- compute_polyfun_prior_weights(
    eqtl_dataset$position, chr, polyfun_data
  )

  # Uniform priors.

  coloc_unif <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
  )

  # eQTLGen estimated density priors.

  coloc_eqtl_tss_eqtlgen <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtl_prior_weights_eqtlgen
  )

  coloc_gwas_tss_eqtlgen <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = eqtl_prior_weights_eqtlgen
  )

  coloc_eqtl_tss_gwas_tss_eqtlgen <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtl_prior_weights_eqtlgen,
    prior_weights2 = eqtl_prior_weights_eqtlgen
  )

  # OneK1K estimated density priors.

  coloc_eqtl_tss_onek1k_round_1 <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtl_prior_weights
  )

  coloc_eqtl_tss_onek1k_round_2 <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtl_prior_weights_round_2
  )

  coloc_eqtl_tss_onek1k_round_3 <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtl_prior_weights_round_3
  )

  coloc_gwas_tss_onek1k_round_1 <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = eqtl_prior_weights
  )

  coloc_eqtl_tss_gwas_tss_onek1k_round_1 <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtl_prior_weights,
    prior_weights2 = eqtl_prior_weights
  )

  # Polyfun priors.

  coloc_polyfun_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = polyfun_prior_weights
  )

  coloc_polyfun_gwas <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = polyfun_prior_weights
  )

  # Gnocchi priors.

  coloc_gnocchi_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = gnocchi_prior_weights
  )

  coloc_gnocchi_gwas <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = gnocchi_prior_weights
  )

  # ABC score priors.

  coloc_abc_score_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = abc_score_prior_weights
  )

  coloc_abc_score_gwas <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = abc_score_prior_weights
  )

  coloc_abc_score_gwas_primary_blood <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = abc_score_primary_blood_prior_weights
  )

  coloc_results <- bind_cols(
    coloc_to_tibble(coloc_unif, "unif"),
    # OneK1K estimated.
    coloc_to_tibble(coloc_eqtl_tss_onek1k_round_1, "eqtl_tss_onek1k_round_1"),
    coloc_to_tibble(coloc_eqtl_tss_onek1k_round_2, "eqtl_tss_onek1k_round_2"),
    coloc_to_tibble(coloc_eqtl_tss_onek1k_round_3, "eqtl_tss_onek1k_round_3"),
    coloc_to_tibble(coloc_gwas_tss_onek1k_round_1, "gwas_tss_onek1k_round_1"),
    coloc_to_tibble(coloc_eqtl_tss_gwas_tss_onek1k_round_1, "eqtl_tss_gwas_tss_onek1k_round_1"),
    # eQTLGen estimated.
    coloc_to_tibble(coloc_eqtl_tss_eqtlgen, "eqtl_tss_eqtlgen"),
    coloc_to_tibble(coloc_gwas_tss_eqtlgen, "gwas_tss_eqtlgen"),
    coloc_to_tibble(coloc_eqtl_tss_gwas_tss_eqtlgen, "eqtl_tss_gwas_tss_eqtlgen"),
    # Gnocchi.
    coloc_to_tibble(coloc_gnocchi_eqtl, "gnocchi_eqtl"),
    coloc_to_tibble(coloc_gnocchi_gwas, "gnocchi_gwas"),
    # ABC score.
    coloc_to_tibble(coloc_abc_score_eqtl, "abc_score_eqtl"),
    coloc_to_tibble(coloc_abc_score_gwas, "abc_score_gwas"),
    coloc_to_tibble(coloc_abc_score_gwas_primary_blood, "abc_score_gwas_primary_blood"),
    # polyfun.
    coloc_to_tibble(coloc_polyfun_eqtl, "polyfun_eqtl"),
    coloc_to_tibble(coloc_polyfun_gwas, "polyfun_gwas"),
    tibble(
      phenotype_id = coloc_metadata$phenotype_id[[i]],
      chromosome = coloc_metadata$chromosome[[i]],
      gene_id = coloc_metadata$gene_id[[i]],
      region = coloc_metadata$region[[i]],
      gene_name = coloc_metadata$gene_name[[i]],
      unif_result = list(coloc_unif$results),
      tss = coloc_metadata$tss[[i]]
  ))

  results[[i]] <- coloc_results
}

write_rds(
  bind_rows(!!!results),
  snakemake@output[["result_file"]]
)
