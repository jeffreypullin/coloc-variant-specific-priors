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
source("code/prior-probabilities-funs.R")

eqtl_file <- snakemake@input[["eqtl_data_file"]]
permutation_file <- snakemake@input[["permutation_file"]]
gwas_file <- snakemake@input[["gwas_data_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]
manifest_file <- snakemake@input[["manifest_file"]]
gwas_id <- snakemake@wildcards[["gwas_id"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])
coloc_results_file <- snakemake@output[["coloc_results_file"]]
finemapping_results_file <- snakemake@output[["finemapping_results_file"]]

gnocchi_data_path <- snakemake@input[["gnocchi_data_path"]]
polyfun_data_1_7_path <- snakemake@input[["polyfun_data_1_7_path"]]
polyfun_data_8_22_path <- snakemake@input[["polyfun_data_8_22_path"]]
abc_score_data_path <- snakemake@input[["abc_score_data_path"]]
eqtlgen_density_path <- snakemake@input[["eqtlgen_density_path"]]
onek1k_r1_density_path <- snakemake@input[["onek1k_r1_density_path"]]
onek1k_r2_density_path <- snakemake@input[["onek1k_r2_density_path"]]
onek1k_r3_density_path <- snakemake@input[["onek1k_r3_density_path"]]

eqtl_metadata <- read_tsv(eqtl_metadata_file, show_col_types = FALSE)
permutation_data <- read_tsv(permutation_file, show_col_types = FALSE)
manifest_data <- read_tsv(manifest_file, show_col_types = FALSE)

permutations <- permutation_data |>
  mutate(FDR = p.adjust(p = p_beta, method = "fdr")) |>
  filter(FDR < 0.01) |>
  select(molecular_trait_object_id, molecular_trait_id) |>
  distinct()

width <- 5e5
coloc_metadata <- eqtl_metadata |>
  filter(gene_type == "protein_coding") |>
  mutate(tss = if_else(strand == 1, gene_start, gene_end)) |>
  rowwise() |>
  mutate(
    start_pos = max(tss - width, 1),
    end_pos = tss + width
  ) |>
  mutate(region = paste0(chromosome, ":", start_pos, "-", end_pos)) |>
  select(gene_name, region, gene_id, chromosome, start_pos, end_pos, tss, gene_name) |>
  ungroup() |>
  filter(chromosome == chr) |>
  filter(!is.na(gene_id)) |>
  filter(gene_id %in% permutations$molecular_trait_id)

all_eqtl_data <- tabix.read.table(eqtl_file, paste0(chr, ":1-2147483647")) |>
  as_tibble() |>
  setNames(eqtl_catalouge_colnames) |>
  left_join(
    eqtl_metadata |>
      select(gene_id, gene_name, gene_type),
    by = join_by(molecular_trait_id == gene_id)
  ) |>
  filter(gene_type == "protein_coding") |>
  mutate(molecular_trait_id = gene_name) |>
  select(-c(gene_name, gene_type))

all_gwas_data <- tabix.read.table(gwas_file, paste0(chr, ":1-2147483647")) |>
  rename(
    rsid = rsids,
    chromosome = chrom,
    position = pos,
    beta = beta,
    se = sebeta,
    maf = af_alt
  ) |>
  mutate(
    molecular_trait_id = gwas_id,
    variant = paste0("chr", chromosome, "_", position, "_", ref, "_", alt)
  ) |>
  select(-c(af_alt_cases, af_alt_controls, ref, alt)) |>
  left_join(
    manifest_data |>
      select(phenocode, num_cases, num_controls),
    by = join_by(molecular_trait_id == phenocode)
  ) |>
  mutate(N = num_cases + num_controls) |>
  select(-c(num_cases, num_controls)) |>
  lazy_dt()

# Remove the HLA region.
if (chr == 6) {
  all_gwas_data <- all_gwas_data |>
    filter(!between(position, 29602228, 33410226))
  all_eqtl_data <- all_eqtl_data |>
    filter(!between(position, 29602228, 33410226))
  coloc_metadata <- coloc_metadata |>
    filter(!between(start_pos, 29602228, 33410226))
}

# Prior information.
gnocchi_data <- read_tsv(gnocchi_data_path, show_col_types = FALSE)
abc_score_data <- read_tsv(abc_score_data_path, show_col_types = FALSE)
density_data_round_1 <- read_rds(onek1k_r1_density_path)
density_data_round_2 <- read_rds(onek1k_r2_density_path)
density_data_round_3 <- read_rds(onek1k_r3_density_path)
eqtlgen_density_data <- read_rds(eqtlgen_density_path)
snp_var_data_1_7 <- read_parquet(polyfun_data_1_7_path)
snp_var_data_8_22 <- read_parquet(polyfun_data_8_22_path)

coloc_results <- list()
finemapping_results <- list()

for (i in seq_len(nrow(coloc_metadata))) {

  region <- paste0(coloc_metadata$chromosome[[i]], ":",
                   coloc_metadata$start_pos[[i]], "-",
                   coloc_metadata$end_pos[[i]])

  print(paste0("Region: ", region))
  print(paste0("Gene ID: ", coloc_metadata$gene_id[[i]]))
  print(paste0("Gene name: ", coloc_metadata$gene_name[[i]]))

  gwas_data <- all_gwas_data |>
    filter(
      position >= coloc_metadata$start_pos[[i]] &
      position <= coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble()

  if (nrow(gwas_data) <= 300) {
    next
  }

  if (min(gwas_data$pval) > 5e-8) {
    next
  }

  eqtl_data <- all_eqtl_data |>
    filter(molecular_trait_id == coloc_metadata$gene_name[[i]]) |>
    filter(
      position >= coloc_metadata$start_pos[[i]] &
      position <= coloc_metadata$end_pos[[i]]
    )

  if (nrow(eqtl_data) <= 300) {
    next
  }

  if (min(eqtl_data$pvalue) > 5e-8) {
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
  is_nvx_na <- all(is.na(2 * gwas_data$N * gwas_data$maf * (1 - gwas_data$maf)))

  if (is_oneover_na || is_nvx_na) {
    next
  }

  gwas_dataset <- list(
    varbeta = gwas_data$se^2,
    N = gwas_data$N,
    MAF = gwas_data$maf,
    type = "cc",
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

  rand_prior_weights <- compute_rand_prior_weights(eqtl_dataset$position)
  shuffled_polyfun_prior_weights <- sample(polyfun_prior_weights)

  # Colocalisation analysis.

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

  # Uncertainty estimation.
  N <- 100
  qs <- c(0, 0.05, 0.2, 0.25, 0.5, 0.75, 0.8, 0.95)

  pph4_rand_eqtl <- numeric(N)
  pph4_rand_gwas <- numeric(N)
  pph4_perm_eqtlgen_eqtl <- numeric(N)
  pph4_perm_eqtlgen_gwas <- numeric(N)
  for (j in seq_len(N)) {

    rand_prior_weights <- compute_rand_prior_weights(eqtl_dataset$position)
    pph4_rand_eqtl[[j]] <- coloc.abf(
      eqtl_dataset,
      gwas_dataset,
      prior_weights1 = rand_prior_weights
    )$summary["PP.H4.abf"]

    pph4_rand_gwas[[j]] <- coloc.abf(
      eqtl_dataset,
      gwas_dataset,
      prior_weights2 = rand_prior_weights
    )$summary["PP.H4.abf"]

    permuted_eqtl_prior_weights <- sample(eqtl_prior_weights)
    pph4_perm_eqtlgen_eqtl[[j]] <- coloc.abf(
      eqtl_dataset,
      gwas_dataset,
      prior_weights1 = permuted_eqtl_prior_weights
    )$summary["PP.H4.abf"]
    pph4_perm_eqtlgen_gwas[[j]] <- coloc.abf(
      eqtl_dataset,
      gwas_dataset,
      prior_weights2 = permuted_eqtl_prior_weights
    )$summary["PP.H4.abf"]
  }

  q_pph4_rand_eqtl <- quantile(pph4_rand_eqtl, qs)
  q_pph4_rand_gwas <- quantile(pph4_rand_gwas, qs)
  q_pph4_perm_eqtlgen_eqtl <- quantile(pph4_perm_eqtlgen_eqtl, qs)
  q_pph4_perm_eqtlgen_gwas <- quantile(pph4_perm_eqtlgen_eqtl, qs)

  colocs <- bind_cols(
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
    # Polyfun.
    coloc_to_tibble(coloc_polyfun_eqtl, "polyfun_eqtl"),
    coloc_to_tibble(coloc_polyfun_gwas, "polyfun_gwas"),
    tibble(
      q_pph4_rand_eqtl = list(q_pph4_rand_eqtl),
      q_pph4_rand_gwas = list(q_pph4_rand_gwas),
      q_pph4_perm_eqtlgen_eqtl = list(q_pph4_perm_eqtlgen_eqtl),
      q_pph4_perm_eqtlgen_gwas = list( q_pph4_perm_eqtlgen_gwas),
      gene_name = coloc_metadata$gene_name[[i]]
  ))
  colocs <- left_join(colocs, coloc_metadata, by = "gene_name")
  coloc_results <- c(coloc_results, list(colocs))

  # Finemapping analysis.

  finemapping_results <- c(finemapping_results, list(tibble(
    # Uniform.
    unif_eqtl = finemap.abf(eqtl_dataset)$SNP.PP,
    unif_gwas = finemap.abf(gwas_dataset)$SNP.PP,
    # eQTLGen.
    eqtlgen_dist_eqtl = finemap.abf(eqtl_dataset, prior_weights = eqtl_prior_weights_eqtlgen)$SNP.PP,
    eqtlgen_dist_gwas = finemap.abf(gwas_dataset, prior_weights = eqtl_prior_weights_eqtlgen)$SNP.PP,
    # OneK1K.
    onek1k_r1_dist_eqtl = finemap.abf(eqtl_dataset, prior_weights = eqtl_prior_weights)$SNP.PP,
    onek1k_r1_dist_gwas = finemap.abf(gwas_dataset, prior_weights = eqtl_prior_weights)$SNP.PP,
    # ABC score.
    abc_score_eqtl = finemap.abf(eqtl_dataset, prior_weights = abc_score_prior_weights)$SNP.PP,
    abc_score_gwas = finemap.abf(gwas_dataset, prior_weights = abc_score_prior_weights)$SNP.PP,
    # Gnocchi.
    gnocchi_eqtl = finemap.abf(eqtl_dataset, prior_weights = gnocchi_prior_weights)$SNP.PP,
    gnocchi_gwas = finemap.abf(gwas_dataset, prior_weights = gnocchi_prior_weights)$SNP.PP,
    # Polyfun.
    polyfun_eqtl = finemap.abf(eqtl_dataset, prior_weights = polyfun_prior_weights)$SNP.PP,
    polyfun_gwas = finemap.abf(gwas_dataset, prior_weights = polyfun_prior_weights)$SNP.PP,
    # Random.
    rand_eqtl = finemap.abf(eqtl_dataset, prior_weights = rand_prior_weights)$SNP.PP,
    rand_gwas = finemap.abf(gwas_dataset, prior_weights = rand_prior_weights)$SNP.PP,
    # Metadata.
    gene_name = coloc_metadata$gene_name[[i]],
  ) |>
    left_join(coloc_metadata, by = "gene_name")
  ))
}

write_rds(
  bind_rows(!!!coloc_results),
  coloc_results_file
)

write_rds(
  bind_rows(!!!finemapping_results),
  finemapping_results_file
)
