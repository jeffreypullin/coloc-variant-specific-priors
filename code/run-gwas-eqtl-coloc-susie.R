# The code in this file is heavily based on:
# https://github.com/ralf-tambets/coloc/blob/main/bin/coloc_v5_pqtl.R

source(here::here("renv/activate.R"))

set.seed(25012024)

suppressPackageStartupMessages({
  library(readr)
  library(seqminer)
  library(dplyr)
  library(janitor)
  library(arrow)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
  devtools::load_all("~/coloc")
})

source("code/coloc-utils.R")
source("code/prior-probabilities-funs.R")

eqtl_lbf_file <- snakemake@input[["eqtl_lbf_file"]]
eqtl_cs_file <- snakemake@input[["eqtl_cs_file"]]

gwas_lbf_file <- snakemake@input[["gwas_lbf_file"]]
gwas_cs_file <- snakemake@input[["gwas_cs_file"]]

permutation_file <- snakemake@input[["permutation_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])

gnocchi_data_path <- snakemake@input[["gnocchi_data_path"]]
polyfun_data_1_7_path <- snakemake@input[["polyfun_data_1_7_path"]]
polyfun_data_8_22_path <- snakemake@input[["polyfun_data_8_22_path"]]
abc_score_data_path <- snakemake@input[["abc_score_data_path"]]
eqtlgen_density_path <- snakemake@input[["eqtlgen_density_path"]]
onek1k_r1_density_path <- snakemake@input[["onek1k_r1_density_path"]]
onek1k_r2_density_path <- snakemake@input[["onek1k_r2_density_path"]]
onek1k_r3_density_path <- snakemake@input[["onek1k_r3_density_path"]]

calculate_region_overlap <- function(region_1, region_2) {

  split_region_1 <- str_extract_all(region_1, "\\d+") |>
    unlist() |>
    as.numeric()
  split_region_2 <- str_extract_all(region_2, "\\d+") |>
    unlist() |>
    as.numeric()

  chr_1 <- split_region_1[[1]]
  chr_2 <- split_region_2[[1]]

  if (chr_1 != chr_2) {
    out <- 0
  } else {
    max_start <- max(split_region_1[[2]], split_region_2[[2]])
    min_end <- min(split_region_1[[3]], split_region_2[[3]])
    overlap_length <- min_end - max_start
    if (overlap_length < 0) {
      out <- 0
    } else {
      out <- overlap_length
    }
  }

  out
}

eqtl_metadata <- read_tsv(eqtl_metadata_file, show_col_types = FALSE)
permutation_data <- read_tsv(permutation_file, show_col_types = FALSE)

permutations <- permutation_data |>
  mutate(FDR = p.adjust(p = p_beta, method = "fdr")) |>
  filter(FDR < 0.01) |>
  select(molecular_trait_object_id, molecular_trait_id) |>
  distinct()

eqtl_cs_data <- read_tsv(eqtl_cs_file, show_col_types = FALSE)
eqtl_lbf_data <- tabix.read.table(eqtl_lbf_file, paste0(chr, ":1-2147483647")) |>
  as_tibble() |>
  setNames(c("molecular_trait_id", "region", "variant", "chromosome", "position", paste0("lbf_variable", 1:10)))

gwas_cs_data <- read_tsv(gwas_cs_file, show_col_types = FALSE)
gwas_lbf_data <- tabix.read.table(gwas_lbf_file, paste0("chr", chr, ":1-2147483647")) |>
  as_tibble() |>
  setNames(c("trait", "region", "v", "rsid", "chromosome", "position", "allele1", "allele2",
             "maf", "beta", "se", "p", "mean", "sd", "prob", "cs", "cs_specific_prob",
             "low_purity", "lead_r2", "mean_99", "sd_99", "prob_99", "cs_99", "cs_specific_prob_99",
             "low_purity_99", "lead_r2_99", paste0("alpha", 1:10),  paste0("mean", 1:10),
             paste0("sd", 1:10), paste0("lbf_variable", 1:10)))

n_regions <- gwas_cs_data |>
  filter(!low_purity) |>
  mutate(chr = str_extract(region, "\\d+")) |>
  filter(chr == !!chr) |>
  nrow()

if (n_regions == 0) {
  coloc_metadata <- tibble()
} else {
  coloc_metadata <- gwas_cs_data |>
    filter(!low_purity) |>
    mutate(chr = str_extract(region, "\\d+")) |>
    filter(chr == !!chr) |>
    left_join(
      eqtl_cs_data |>
        distinct(region, gene_id) |>
        mutate(chr = str_extract(region, "\\d+")),
      by = "chr",
      relationship = "many-to-many"
    ) |>
    rowwise() |>
    mutate(overlap = calculate_region_overlap(region.x, region.y)) |>
    ungroup() |>
    filter(overlap > 0) |>
    rename(
      gwas_region = region.x,
      eqtl_region = region.y
    ) |>
    filter(gene_id %in% permutations$molecular_trait_id) |>
    left_join(
      eqtl_metadata |>
        select(gene_id, gene_name, chromosome, gene_type, gene_start,
              gene_end, strand),
      by = "gene_id"
    ) |>
    filter(gene_type == "protein_coding") |>
    mutate(tss = if_else(strand == 1, gene_start, gene_end))
}

gnocchi_data <- read_tsv(gnocchi_data_path, show_col_types = FALSE)
abc_score_data <- read_tsv(abc_score_data_path, show_col_types = FALSE)
density_data_round_1 <- read_rds(onek1k_r1_density_path)
density_data_round_2 <- read_rds(onek1k_r2_density_path)
density_data_round_3 <- read_rds(onek1k_r3_density_path)
eqtlgen_density_data <- read_rds(eqtlgen_density_path)
snp_var_data_1_7 <- read_parquet(polyfun_data_1_7_path)
snp_var_data_8_22 <- read_parquet(polyfun_data_8_22_path)

results <- list()
for (i in seq_len(nrow(coloc_metadata))) {

  print(paste0("GWAS Region: ", coloc_metadata$gwas_region[[i]]))
  print(paste0("eQTL Region: ", coloc_metadata$eqtl_region[[i]]))
  print(paste0("Gene ID: ", coloc_metadata$gene_id[[i]]))
  print(paste0("Gene name: ", coloc_metadata$gene_name[[i]]))

  region_eqtl_cs_data <- eqtl_cs_data |>
    filter(
      region == coloc_metadata$eqtl_region[[i]],
      molecular_trait_id == coloc_metadata$gene_id[[i]]
    )

  if (nrow(region_eqtl_cs_data) == 0) {
    next
  }

  region_gwas_cs_data <- gwas_cs_data |>
    filter(region == coloc_metadata$gwas_region[[i]])

  if (nrow(region_gwas_cs_data) == 0) {
    next
  }

  eqtl_cs_ids <- region_eqtl_cs_data |>
    pull(cs_id) |>
    unique()
  gwas_cs_ids <- region_gwas_cs_data |>
    pull(cs) |>
    unique()

  split_eqtl_cs_ids <- str_split(eqtl_cs_ids, "_")
  gene_ids <- unique(map_chr(split_eqtl_cs_ids, 1))
  eqtl_cs_inds <- map_chr(split_eqtl_cs_ids, \(x) str_sub(x[[2]], 2, 2)) |>
    as.numeric() |>
    order(decreasing = TRUE)
  eqtl_cs_col_names <- paste0("lbf_variable", eqtl_cs_inds)

  gwas_cs_col_names <- paste0("lbf_variable", gwas_cs_ids)

  eqtl_lbf_mat <- eqtl_lbf_data |>
    filter(
      region == coloc_metadata$eqtl_region[[i]],
      molecular_trait_id == coloc_metadata$gene_id[[i]]
    ) |>
    select(variant, !!eqtl_cs_col_names) |>
    column_to_rownames("variant") |>
    as.matrix() |>
    t()

  gwas_lbf_mat <- gwas_lbf_data |>
    filter(region == coloc_metadata$gwas_region[[i]]) |>
    select(rsid, !!gwas_cs_col_names) |>
    column_to_rownames("rsid") |>
    as.matrix() |>
    t()

  # FIXME: Which position should be used?
  position <- eqtl_lbf_data |>
    filter(molecular_trait_id %in% gene_ids) |>
    pull(position)

  eqtl_prior_weights_eqtlgen <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$tss[[i]], eqtlgen_density_data
  )
  eqtl_prior_weights <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$tss[[i]], density_data_round_1
  )
  eqtl_prior_weights_round_2 <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$tss[[i]], density_data_round_2
  )
  eqtl_prior_weights_round_3 <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$tss[[i]], density_data_round_3
  )
  gnocchi_prior_weights <- compute_gnocchi_prior_weights(
    position, chr, gnocchi_data
  )
  abc_score_prior_weights <- compute_abc_prior_weights(
    position, chr, coloc_metadata$gene_name[[i]], abc_score_data
  )
  abc_score_primary_blood_prior_weights <- compute_abc_prior_weights(
    position, chr, coloc_metadata$gene_name[[i]], abc_score_data,
    biosamples = "primary_blood"
  )

  if (chr %in% 1:7) {
    polyfun_data <- snp_var_data_1_7
  } else {
    polyfun_data <- snp_var_data_8_22
  }

  if (length(intersect(colnames(eqtl_lbf_mat), colnames(gwas_lbf_mat))) == 0) {
    next
  }

  polyfun_prior_weights <- compute_polyfun_prior_weights(
    position, chr, polyfun_data
  )

  # Uniform priors.

  coloc_unif <- coloc.bf_bf(eqtl_lbf_mat, gwas_lbf_mat)

  # eQTLGen estimated density priors.

  coloc_eqtl_tss_eqtlgen <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = eqtl_prior_weights_eqtlgen
  )

  coloc_gwas_tss_eqtlgen <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights2 = eqtl_prior_weights_eqtlgen
  )

  coloc_eqtl_tss_gwas_tss_eqtlgen <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = eqtl_prior_weights_eqtlgen,
    prior_weights2 = eqtl_prior_weights_eqtlgen
  )

  # OneK1K estimated density priors.

  coloc_eqtl_tss_onek1k_round_1 <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = eqtl_prior_weights
  )

  coloc_eqtl_tss_onek1k_round_2 <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = eqtl_prior_weights_round_2
  )

  coloc_eqtl_tss_onek1k_round_3 <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = eqtl_prior_weights_round_3
  )

  coloc_gwas_tss_onek1k_round_1 <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights2 = eqtl_prior_weights
  )

  coloc_eqtl_tss_gwas_tss_onek1k_round_1 <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = eqtl_prior_weights,
    prior_weights2 = eqtl_prior_weights
  )

  # Polyfun priors.

  coloc_polyfun_eqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = polyfun_prior_weights
  )

  coloc_polyfun_gwas <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights2 = polyfun_prior_weights
  )

  # Gnocchi priors.

  coloc_gnocchi_eqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = gnocchi_prior_weights
  )

  coloc_gnocchi_gwas <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights2 = gnocchi_prior_weights
  )

  # ABC score priors.

  coloc_abc_score_eqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights1 = abc_score_prior_weights
  )

  coloc_abc_score_gwas <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
    prior_weights2 = abc_score_prior_weights
  )

  coloc_abc_score_gwas_primary_blood <- coloc.bf_bf(
    eqtl_lbf_mat,
    gwas_lbf_mat,
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
    # PolyFun.
    coloc_to_tibble(coloc_polyfun_eqtl, "polyfun_eqtl"),
    coloc_to_tibble(coloc_polyfun_gwas, "polyfun_gwas"),
    tibble(
      chromosome = coloc_metadata$chromosome[[i]],
      gene_id = coloc_metadata$gene_id[[i]],
      gwas_region = coloc_metadata$gwas_region[[i]],
      eqtl_region = coloc_metadata$eqtl_region[[i]],
      gene_name = coloc_metadata$gene_name[[i]],
      tss = coloc_metadata$tss[[i]]
  ))

  results[[i]] <- coloc_results
}

all_results <- bind_rows(!!!results)
# NAs caused by 'snp overlap too small' issue.
if (nrow(all_results) != 0) {
  all_results <- all_results |>
    filter(!is.na(PP.H4.abf_unif))
}

write_rds(
  all_results,
  snakemake@output[["result_file"]]
)