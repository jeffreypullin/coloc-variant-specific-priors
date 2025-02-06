
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(arrow)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(seqminer)
  library(janitor)
  library(purrr)
  library(data.table)
  library(dtplyr)
  devtools::load_all("~/coloc")
})

source("code/coloc-utils.R")
source("code/prior-probabilities-funs.R")

define_loci <- function(sig_pos, sig_chr, sig_pval) {

  stopifnot(length(sig_pos) == length(sig_pval))

  new_sig_pos <- sig_pos
  new_sig_pval <- sig_pval
  new_sig_chr <- sig_chr

  for (k in 1:5) {
    sig_pos <- new_sig_pos
    sig_pval <- new_sig_pval
    sig_chr <- new_sig_chr

    new_sig_pos <- numeric()
    new_sig_chr <- numeric()
    new_sig_pval <- numeric()
    for (i in seq_along(sig_pos)) {
      snp_pos <- sig_pos[[i]]
      snp_chr <- sig_chr[[i]]

      window_ind <- which(
        (sig_chr == snp_chr) &
        (sig_pos < snp_pos + 1e6) &
        (sig_pos > snp_pos - 1e6)
      )
      min_window_pval_ind <- which.min(sig_pval[window_ind])

      min_window_pval <- sig_pval[window_ind[min_window_pval_ind]]
      min_window_pval_pos <- sig_pos[window_ind[min_window_pval_ind]]
      min_window_pval_chr <- sig_chr[window_ind[min_window_pval_ind]]

      new_sig_chr <- c(new_sig_chr, min_window_pval_chr)
      new_sig_pos <- c(new_sig_pos, min_window_pval_pos)
      new_sig_pval <- c(new_sig_pval, min_window_pval)
    }
    variant <- paste0(new_sig_chr, "-", new_sig_pos, "-", new_sig_pval)
    unique_variant <- unique(variant)
    unique_ind <- numeric(length(unique_variant))
    for (j in seq_along(unique_variant)) {
      unique_ind[[j]] <- which.max(variant == unique_variant[[j]])
    }

    new_sig_chr <- new_sig_chr[unique_ind]
    new_sig_pos <- new_sig_pos[unique_ind]
    new_sig_pval <- new_sig_pval[unique_ind]
  }
  tibble(chrom = new_sig_chr, pos = new_sig_pos, pval = new_sig_pval)
}

ukbb_manifest <- tibble(
  icd_code = c("244", "250.2", "401"),
  num_cases = c(14871, 18945, 77977),
  num_controls = c(391429, 388756,  330366)
)

gwas_file <- snakemake@input[["gwas_data_file"]]
eqtl_file <- snakemake@input[["eqtl_data_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]

icd_code_val <- snakemake@wildcards[["icd_code"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])

polyfun_trait_specific_path <- snakemake@input[["polyfun_trait_specific_data_path"]]
polyfun_data_1_7_path <- snakemake@input[["polyfun_data_1_7_path"]]
polyfun_data_8_22_path <- snakemake@input[["polyfun_data_8_22_path"]]
abc_score_data_path <- snakemake@input[["abc_score_data_path"]]
eqtlgen_density_path <- snakemake@input[["eqtlgen_density_path"]]
onek1k_r1_density_path <- snakemake@input[["onek1k_r1_density_path"]]
onek1k_r2_density_path <- snakemake@input[["onek1k_r2_density_path"]]

eqtl_metadata <- read_tsv(eqtl_metadata_file, show_col_types = FALSE)

test <- tabix.read.table(gwas_file, paste0(chr, ":1-2147483647")) |>
  setNames(c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
             "consequence", "pval", "beta", "sebeta", "af", "ac", "tstat")) |>
  as_tibble()

all_gwas_data <- tabix.read.table(gwas_file, paste0(chr, ":1-2147483647")) |>
  setNames(c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
             "consequence", "pval", "beta", "sebeta", "af", "ac", "tstat")) |>
  rename(
    rsid = rsids,
    chromosome = chrom,
    position = pos,
    beta = beta,
    se = sebeta,
    maf = af
  ) |>
  mutate(
    icd_code = icd_code_val,
    variant = paste0("chr", chromosome, "_", position, "_", ref, "_", alt)
  ) |>
  select(-c(ref, alt)) |>
  left_join(
    ukbb_manifest |>
      select(icd_code, num_cases, num_controls),
    by = join_by(icd_code == icd_code)
  ) |>
  mutate(N = num_cases + num_controls) |>
  select(-c(num_cases, num_controls)) |>
  # Why does SAIGE produce p-values less than 0?
  filter(pval > 0) |>
  as_tibble()

sig_gwas_data <- all_gwas_data |>
  filter(pval < 5 * 10^-8)

gwas_loci <- define_loci(
  sig_gwas_data$position,
  sig_gwas_data$chromosome,
  sig_gwas_data$pval
) |>
  mutate(
    start_pos = pos - 1e6,
    end_pos = pos + 1e6,
    region = paste0(chrom, ":", start_pos, "-", end_pos),
  ) |>
  arrange(pval)

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

coloc_metadata <- gwas_loci |>
  left_join(
    all_eqtl_data |>
      select(position, molecular_trait_id),
    join_by(between(y$position, x$start_pos, x$end_pos))
  ) |>
  # Filter out loci with no eQTL signal.
  filter(!is.na(molecular_trait_id)) |>
  distinct(chrom, start_pos, end_pos, region, molecular_trait_id) |>
  left_join(
    eqtl_metadata |>
      select(gene_name, gene_start, gene_end, strand),
    join_by(molecular_trait_id == gene_name)
  ) |>
  mutate(tss = if_else(strand == 1, gene_start, gene_end))

all_gwas_data <- lazy_dt(all_gwas_data)

trait_specific_polyfun_data <- read_tsv(polyfun_trait_specific_path, show_col_types = FALSE)
density_data_round_1 <- read_rds(onek1k_r1_density_path)
eqtlgen_density_data <- read_rds(eqtlgen_density_path)
snp_var_data_1_7 <- read_parquet(polyfun_data_1_7_path)
snp_var_data_8_22 <- read_parquet(polyfun_data_8_22_path)

coloc_results <- list()
finemapping_results <- list()
for (i in seq_len(nrow(coloc_metadata))) {
  
  print(paste0("Region: ", coloc_metadata$region[[i]]))
  print(paste0("Gene: ", coloc_metadata$molecular_trait_id[[i]]))

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
    filter(molecular_trait_id == coloc_metadata$molecular_trait_id[[i]]) |>
    filter(
      position >= coloc_metadata$start_pos[[i]] &
      position <= coloc_metadata$end_pos[[i]]
    )

  if (nrow(eqtl_data) <= 300) {
    next
  }

  if (min(eqtl_data$pvalue) > 5e-6) {
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

  position <- eqtl_dataset$position
  tss <- coloc_metadata$tss[[i]]

  eqtlgen_prior_weights <- compute_eqtl_tss_dist_prior_weights(
    position, tss, eqtlgen_density_data
  )
  onek1k_r1_prior_weights <- compute_eqtl_tss_dist_prior_weights(
    position, tss, density_data_round_1
  )

  if (chr %in% 1:7) {
    precomputed_polyfun_data <- snp_var_data_1_7
  } else {
    precomputed_polyfun_data <- snp_var_data_8_22
  }

  precomputed_polyfun_prior_weights <- compute_polyfun_prior_weights(
    eqtl_dataset$position, chr, precomputed_polyfun_data
  )

  trait_specific_polyfun_prior_weights <- compute_polyfun_trait_specific_prior_weights(
    eqtl_dataset$position, chr, trait_specific_polyfun_data
  )

  # Colocalisation analysis.

  # Uniform.

  coloc_unif <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    p12 = 5e-6
  )

  # eQTLGen.

  coloc_eqtl_tss_eqtlgen <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtlgen_prior_weights,
    p12 = 5e-6
  )

  # Precomputed Polyfun.

  coloc_polyfun_precomputed <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = precomputed_polyfun_prior_weights,
    p12 = 5e-6
  )

  # Precomputed Polyfun and eQTLGen.
  coloc_polyfun_precomputed_eqtlgen <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtlgen_prior_weights,
    prior_weights2 = precomputed_polyfun_prior_weights,
    p12 = 5e-6
  )

  # Trait-specific Polyfun.

  coloc_polyfun_trait_specific <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights2 = trait_specific_polyfun_prior_weights,
    p12 = 5e-6
  )

  # Trait-specific Polyfun and eQTLGen.

  coloc_polyfun_trait_specific_eqtlgen <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    prior_weights1 = eqtlgen_prior_weights,
    prior_weights2 = trait_specific_polyfun_prior_weights,
    p12 = 5e-6
  )

  colocs <- bind_cols(
    coloc_to_tibble(coloc_unif, "unif"),
    coloc_to_tibble(coloc_eqtl_tss_eqtlgen, "eqtlgen"),
    coloc_to_tibble(coloc_polyfun_precomputed, "polyfun_precomputed"),
    coloc_to_tibble(coloc_polyfun_trait_specific, "polyfun_trait_specific"),
    coloc_to_tibble(coloc_polyfun_precomputed_eqtlgen, "polyfun_precomputed_eqtlgen"),
    coloc_to_tibble(coloc_polyfun_trait_specific_eqtlgen, "polyfun_trait_specific_eqtlgen"),
    tibble(molecular_trait_id = coloc_metadata$molecular_trait_id[[i]])
  )
  colocs <- left_join(colocs, coloc_metadata, by = "molecular_trait_id")
  coloc_results <- c(coloc_results, list(colocs))

  finemapping_results <- c(finemapping_results, list(tibble(
    unif = finemap.abf(gwas_dataset)$SNP.PP,
    eqtlgen = finemap.abf(gwas_dataset, prior_weights = eqtlgen_prior_weights)$SNP.PP,
    polyfun_precomputed = finemap.abf(gwas_dataset, prior_weights = precomputed_polyfun_prior_weights)$SNP.PP,
    polyfun_trait_specific = finemap.abf(gwas_dataset, prior_weights = trait_specific_polyfun_prior_weights)$SNP.PP,
    molecular_trait_id = coloc_metadata$molecular_trait_id[[i]],
    eqtlgen_max = max(eqtlgen_prior_weights),
    precomputed_polyfun_max = max(precomputed_polyfun_prior_weights),
    trait_specific_polyfun_max = max(trait_specific_polyfun_prior_weights)
  ) |>
    left_join(coloc_metadata, by = "molecular_trait_id")
  ))
}

write_rds(
  bind_rows(!!!coloc_results),
  snakemake@output[["coloc_results_file"]]
)

write_rds(
  bind_rows(!!!finemapping_results),
  snakemake@output[["finemapping_results_file"]]
)
