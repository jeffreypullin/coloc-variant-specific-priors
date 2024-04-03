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
  devtools::load_all("~/coloc")
})

source("code/coloc-utils.R")

eqtl_lbf_file <- snakemake@input[["eqtl_lbf_file"]]
eqtl_cs_file <- snakemake@input[["eqtl_cs_file"]]

pqtl_lbf_file <- snakemake@input[["pqtl_lbf_file"]]
pqtl_cs_file <- snakemake@input[["pqtl_cs_file"]]

permutation_file <- snakemake@input[["permutation_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]
pqtl_metadata_file <- snakemake@input[["pqtl_metadata_file"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])

eqtl_metadata <- read_tsv(eqtl_metadata_file, show_col_types = FALSE)
pqtl_metadata <- read_tsv(pqtl_metadata_file, show_col_types = FALSE)
permutation_data <- read_tsv(permutation_file, show_col_types = FALSE)

permutations <- permutation_data |>
  mutate(FDR = p.adjust(p = p_beta, method = "fdr")) |>
  filter(FDR < 0.01) |>
  select(molecular_trait_object_id, molecular_trait_id) |>
  distinct()

width <- 5e5
coloc_metadata <- eqtl_metadata |>
  filter(gene_type == "protein_coding") |>
  mutate(gene_tss = if_else(strand == 1, gene_start, gene_end)) |>
  rowwise() |>
  mutate(
    start_pos = max(gene_tss - width, 1),
    end_pos = gene_tss + width
  ) |>
  mutate(region = paste0(chromosome, ":", start_pos, "-", end_pos)) |>
  select(gene_name, region, gene_id, chromosome, start_pos, end_pos, gene_tss) |>
  ungroup() |>
  left_join(
    pqtl_metadata |>
      mutate(protein_tss = if_else(strand == 1, gene_start, gene_end)) |>
      select(phenotype_id, chromosome, phenotype_pos, pqtl_gene_id = gene_id,
             phenotype_start = gene_start, phenotype_end = gene_end,
             protein_tss) |>
      mutate(chromosome = as.character(chromosome)),
    by = join_by(chromosome == chromosome,
                 start_pos <= phenotype_pos,
                 end_pos >= phenotype_pos)
  ) |>
  select(gene_name, gene_id, phenotype_id, region, chromosome,
         start_pos, end_pos, gene_tss, pqtl_gene_id, phenotype_pos,
         phenotype_start, phenotype_end, protein_tss) |>
  filter(chromosome == chr) |>
  filter(!is.na(phenotype_id)) |>
  filter(!is.na(gene_id)) |>
  filter(gene_id %in% permutations$molecular_trait_id)

eqtl_cs_data <- read_tsv(eqtl_cs_file, show_col_types = FALSE)
eqtl_lbf_data <- tabix.read.table(eqtl_lbf_file, paste0(chr, ":1-2147483647")) |>
  as_tibble() |>
  setNames(c("molecular_trait_id", "region", "variant", "chromosome", "position", paste0("lbf_variable", 1:10)))

pqtl_cs_data <- read_tsv(pqtl_cs_file, show_col_types = FALSE)
pqtl_lbf_data <- tabix.read.table(pqtl_lbf_file, paste0(chr, ":1-2147483647")) |>
  as_tibble() |>
  setNames(c("molecular_trait_id", "region", "variant", "chromosome", "position", paste0("lbf_variable", 1:10)))

gnocchi_data <- read_tsv("data/gnocchi-windows.bed",
                         col_names = FALSE, show_col_types = FALSE)
colnames(gnocchi_data) <- c("chromosome", "start_pos", "end_pos", "score")

abc_score_data <- read_tsv("data/abc-data.txt.gz", show_col_types = FALSE)

density_data_round_1 <- read_rds("output/densities/onek1k_cd4nc_round_1.rds")
density_data_round_2 <- read_rds("output/densities/onek1k_cd4nc_round_2.rds")
density_data_round_3 <- read_rds("output/densities/onek1k_cd4nc_round_3.rds")
eqtlgen_density_data <- read_rds("output/densities/eqtlgen.rds")

snp_var_data_1_7 <- read_parquet("data/snpvar_meta.chr1_7.parquet")
snp_var_data_8_22 <- read_parquet("data/snpvar_meta.chr8_22.parquet")

results <- list()
for (i in seq_len(nrow(coloc_metadata))) {

  region <- paste0(coloc_metadata$chromosome[[i]], ":",
                   coloc_metadata$start_pos[[i]], "-",
                   coloc_metadata$end_pos[[i]])

  print(paste0("Region: ", region))
  print(paste0("Gene ID: ", coloc_metadata$gene_id[[i]]))
  print(paste0("Gene name: ", coloc_metadata$gene_name[[i]]))
  print(paste0("Protein ID: ", coloc_metadata$phenotype_id[[i]]))

  gene_eqtl_cs_data <- eqtl_cs_data |>
    filter(molecular_trait_id == coloc_metadata$gene_id[[i]])

  if (nrow(gene_eqtl_cs_data) == 0) {
    next
  }

  gene_pqtl_cs_data <- pqtl_cs_data |>
    filter(molecular_trait_id == coloc_metadata$phenotype_id[[i]])

  if (nrow(gene_pqtl_cs_data) == 0) {
    next
  }

  eqtl_cs_ids <- gene_eqtl_cs_data |>
    pull(cs_id) |>
    unique()
  pqtl_cs_ids <- gene_pqtl_cs_data |>
    pull(cs_id) |>
    unique()

  split_pqtl_cs_ids <- str_split(pqtl_cs_ids, "_")
  protein_ids <- unique(map_chr(split_pqtl_cs_ids, 1))
  pqtl_cs_inds <- map_chr(split_pqtl_cs_ids, \(x) str_sub(x[[2]], 2, 2)) |>
    as.numeric() |>
    order(decreasing = TRUE)
  pqtl_cs_col_names <- paste0("lbf_variable", pqtl_cs_inds)

  split_eqtl_cs_ids <- str_split(eqtl_cs_ids, "_")
  gene_ids <- unique(map_chr(split_eqtl_cs_ids, 1))
  eqtl_cs_inds <- map_chr(split_eqtl_cs_ids, \(x) str_sub(x[[2]], 2, 2)) |>
    as.numeric() |>
    order(decreasing = TRUE)
  eqtl_cs_col_names <- paste0("lbf_variable", eqtl_cs_inds)

  pqtl_lbf_mat <- pqtl_lbf_data |>
    filter(molecular_trait_id %in% protein_ids) |>
    select(variant, !!pqtl_cs_col_names) |>
    column_to_rownames("variant") |>
    as.matrix() |>
    t()

  eqtl_lbf_mat <- eqtl_lbf_data |>
    filter(molecular_trait_id %in% gene_ids) |>
    select(variant, !!eqtl_cs_col_names) |>
    column_to_rownames("variant") |>
    as.matrix() |>
    t()

  # FIXME: Which position should be used?
  position <- eqtl_lbf_data |>
    filter(molecular_trait_id %in% gene_ids) |>
    pull(position)

  eqtl_prior_weights_eqtlgen <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$gene_tss[[i]], eqtlgen_density_data
  )
  eqtl_prior_weights <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$gene_tss[[i]], density_data_round_1
  )
  eqtl_prior_weights_round_2 <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$gene_tss[[i]], density_data_round_2
  )
  eqtl_prior_weights_round_3 <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$gene_tss[[i]], density_data_round_3
  )

  pqtl_prior_weights_eqtlgen <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$protein_tss[[i]], eqtlgen_density_data
  )
  pqtl_prior_weights <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$protein_tss[[i]], density_data_round_1
  )
  pqtl_prior_weights_round_2 <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$protein_tss[[i]], density_data_round_2
  )
  pqtl_prior_weights_round_3 <- compute_eqtl_tss_dist_prior_weights(
    position, coloc_metadata$protein_tss[[i]], density_data_round_3
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

  polyfun_prior_weights <- compute_polyfun_prior_weights(
    position, chr, polyfun_data
  )

  # Uniform priors.

  coloc_unif <- coloc.bf_bf(eqtl_lbf_mat, pqtl_lbf_mat)

  # eQTLGen estimated density priors.

  coloc_eqtl_tss_eqtlgen <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = eqtl_prior_weights_eqtlgen
  )

  coloc_pqtl_tss_eqtlgen <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights2 = pqtl_prior_weights_eqtlgen
  )

  coloc_eqtl_tss_pqtl_tss_eqtlgen <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = eqtl_prior_weights_eqtlgen,
    prior_weights2 = pqtl_prior_weights_eqtlgen
  )

  # OneK1K estimated density priors.

  coloc_eqtl_tss_onek1k_round_1 <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = eqtl_prior_weights
  )

  coloc_eqtl_tss_onek1k_round_2 <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = eqtl_prior_weights_round_2
  )

  coloc_eqtl_tss_onek1k_round_3 <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = eqtl_prior_weights_round_3
  )

  coloc_pqtl_tss_onek1k_round_1 <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights2 = pqtl_prior_weights
  )

  coloc_eqtl_tss_pqtl_tss_onek1k_round_1 <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = eqtl_prior_weights,
    prior_weights2 = pqtl_prior_weights
  )

  # Polyfun priors.

  coloc_polyfun_eqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = polyfun_prior_weights
  )

  coloc_polyfun_pqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights2 = polyfun_prior_weights
  )

  # Gnocchi priors.

  coloc_gnocchi_eqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = gnocchi_prior_weights
  )

  coloc_gnocchi_pqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights2 = gnocchi_prior_weights
  )

  # ABC score priors.

  coloc_abc_score_eqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights1 = abc_score_prior_weights
  )

  coloc_abc_score_pqtl <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights2 = abc_score_prior_weights
  )

  coloc_abc_score_pqtl_primary_blood <- coloc.bf_bf(
    eqtl_lbf_mat,
    pqtl_lbf_mat,
    prior_weights2 = abc_score_primary_blood_prior_weights
  )

  coloc_results <- bind_cols(
    coloc_to_tibble(coloc_unif, "unif"),
    # OneK1K estimated.
    coloc_to_tibble(coloc_eqtl_tss_onek1k_round_1, "eqtl_tss_onek1k_round_1"),
    coloc_to_tibble(coloc_eqtl_tss_onek1k_round_2, "eqtl_tss_onek1k_round_2"),
    coloc_to_tibble(coloc_eqtl_tss_onek1k_round_3, "eqtl_tss_onek1k_round_3"),
    coloc_to_tibble(coloc_pqtl_tss_onek1k_round_1, "pqtl_tss_onek1k_round_1"),
    coloc_to_tibble(coloc_eqtl_tss_pqtl_tss_onek1k_round_1, "eqtl_tss_pqtl_tss_onek1k_round_1"),
    # eQTLGen estimated.
    coloc_to_tibble(coloc_eqtl_tss_eqtlgen, "eqtl_tss_eqtlgen"),
    coloc_to_tibble(coloc_pqtl_tss_eqtlgen, "pqtl_tss_eqtlgen"),
    coloc_to_tibble(coloc_eqtl_tss_pqtl_tss_eqtlgen, "eqtl_tss_pqtl_tss_eqtlgen"),
    # Gnocchi.
    coloc_to_tibble(coloc_gnocchi_eqtl, "gnocchi_eqtl"),
    coloc_to_tibble(coloc_gnocchi_pqtl, "gnocchi_pqtl"),
    # ABC score.
    coloc_to_tibble(coloc_abc_score_eqtl, "abc_score_eqtl"),
    coloc_to_tibble(coloc_abc_score_pqtl, "abc_score_pqtl"),
    coloc_to_tibble(coloc_abc_score_pqtl_primary_blood, "abc_score_pqtl_primary_blood"),
    # polyfun.
    coloc_to_tibble(coloc_polyfun_eqtl, "polyfun_eqtl"),
    coloc_to_tibble(coloc_polyfun_pqtl, "polyfun_pqtl"),
    tibble(
      phenotype_id = coloc_metadata$phenotype_id[[i]],
      chromosome = coloc_metadata$chromosome[[i]],
      gene_id = coloc_metadata$gene_id[[i]],
      region = coloc_metadata$region[[i]],
      gene_name = coloc_metadata$gene_name[[i]],
      gene_tss = coloc_metadata$gene_tss[[i]],
      protein_tss = coloc_metadata$protein_tss[[i]]
  ))

  results[[i]] <- coloc_results
}

write_rds(
  bind_rows(!!!results),
  snakemake@output[["result_file"]]
)