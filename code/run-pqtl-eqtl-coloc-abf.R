# The code in this file is heavily based on:
# https://github.com/ralf-tambets/coloc/blob/main/bin/coloc_v3_pqtl.R

source(here::here("renv/activate.R"))

set.seed(25012024)

suppressPackageStartupMessages({
  library(readr)
  library(seqminer)
  library(dplyr)
  library(janitor)
  library(arrow)
  library(qvalue)
  devtools::load_all("~/coloc")
})

source("code/coloc-utils.R")
source("code/prior-probabilities-funs.R")

eqtl_file <- snakemake@input[["eqtl_data_file"]]
pqtl_file <- snakemake@input[["pqtl_data_file"]]
permutation_file <- snakemake@input[["permutation_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]
pqtl_metadata_file <- snakemake@input[["pqtl_metadata_file"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])

eqtl_metadata <- read_tsv(eqtl_metadata_file, show_col_types = FALSE)
pqtl_metadata <- read_tsv(pqtl_metadata_file, show_col_types = FALSE)
permutation_data <- read_tsv(permutation_file, show_col_types = FALSE)

gnocchi_data_path <- snakemake@input[["gnocchi_data_path"]]
polyfun_data_1_7_path <- snakemake@input[["polyfun_data_1_7_path"]]
polyfun_data_8_22_path <- snakemake@input[["polyfun_data_8_22_path"]]
abc_score_data_path <- snakemake@input[["abc_score_data_path"]]
eqtlgen_density_path <- snakemake@input[["eqtlgen_density_path"]]
onek1k_r1_density_path <- snakemake@input[["onek1k_r1_density_path"]]
onek1k_r2_density_path <- snakemake@input[["onek1k_r2_density_path"]]

permutations <- permutation_data |>
  mutate(qvalue = qvalue(p = p_beta)$qvalue) |>
  filter(qvalue < 0.05) |>
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
  select(gene_name, region, gene_id, chromosome, start_pos, end_pos, 
         gene_tss) |>
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

all_eqtl_data <- tabix.read.table(eqtl_file, paste0(chr, ":1-2147483647")) |>
  as_tibble()
all_pqtl_data <- tabix.read.table(pqtl_file, paste0(chr, ":1-2147483647")) |>
  as_tibble()

gnocchi_data <- read_tsv(gnocchi_data_path, show_col_types = FALSE)
abc_score_data <- read_tsv(abc_score_data_path, show_col_types = FALSE)
eqtlgen_density_data <- read_rds(eqtlgen_density_path)
onek1k_r1_density_data <- read_rds(onek1k_r1_density_path)
onek1k_r2_density_data <- read_rds(onek1k_r2_density_path)
if (chr %in% 1:7) {
  polyfun_data <- read_parquet(polyfun_data_1_7_path)
} else {
  polyfun_data <- read_parquet(polyfun_data_8_22_path)
}

results <- list()
for (i in seq_len(nrow(coloc_metadata))) {

  region <- paste0(coloc_metadata$chromosome[[i]], ":",
                   coloc_metadata$start_pos[[i]], "-",
                   coloc_metadata$end_pos[[i]])

  print(paste0("Region: ", region))
  print(paste0("Gene ID: ", coloc_metadata$gene_id[[i]]))
  print(paste0("Gene name: ", coloc_metadata$gene_name[[i]]))
  print(paste0("Protein ID: ", coloc_metadata$phenotype_id[[i]]))

  eqtl_data <- all_eqtl_data |>
    setNames(eqtl_catalouge_colnames) |>
    filter_qtl_dataset(
      trait_id = coloc_metadata$gene_id[[i]],
      chrom = coloc_metadata$chromosome[[i]],
      start_pos = coloc_metadata$start_pos[[i]],
      end_pos = coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble()

  pqtl_data <- all_pqtl_data |>
    setNames(eqtl_catalouge_colnames) |>
    filter_qtl_dataset(
      trait_id = coloc_metadata$phenotype_id[[i]],
      chrom = coloc_metadata$chromosome[[i]],
      start_pos = coloc_metadata$start_pos[[i]],
      end_pos = coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble()

  if (nrow(eqtl_data) == 0 || nrow(pqtl_data) == 0) {
    next
  }

  if (min(eqtl_data$pvalue) > 5e-6) {
    next
  }
  if (min(pqtl_data$pvalue) > 5e-8) {
    next
  }

  pqtl_data <- prepare_coloc_dataset(pqtl_data)
  eqtl_data <- prepare_coloc_dataset(eqtl_data)

  n_before <- max(nrow(eqtl_data), nrow(pqtl_data))
  if (n_before <= 300) {
    next
  }

  eqtl_data <- eqtl_data |>
    filter(variant %in% pqtl_data$variant)
  pqtl_data <- pqtl_data |>
    filter(variant %in% eqtl_data$variant)

  n_after <- max(nrow(eqtl_data), nrow(pqtl_data))
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

  is_oneover_na <- all(is.na(1 / pqtl_data$se^2))
  is_nvx_na <- all(is.na(2 * pqtl_data$an / 2 * pqtl_data$maf * (1 - pqtl_data$maf)))

  if (is_oneover_na || is_nvx_na) {
    next
  }

  pqtl_dataset <- list(
    varbeta = pqtl_data$se^2,
    N = pqtl_data$an / 2,
    MAF = pqtl_data$maf,
    type = "quant",
    beta = pqtl_data$beta,
    snp = pqtl_data$variant
  )

  gene_tss <- coloc_metadata$gene_tss[[i]]
  protein_tss <- coloc_metadata$protein_tss[[i]]
  position <- eqtl_dataset$position

  eqtlgen_prior_weights_eqtl <- compute_eqtl_tss_dist_prior_weights(
    position, gene_tss, eqtlgen_density_data
  )
  onek1k_r1_prior_weights_eqtl <- compute_eqtl_tss_dist_prior_weights(
    position, gene_tss, onek1k_r1_density_data
  )
  onek1k_r2_prior_weights_eqtl <- compute_eqtl_tss_dist_prior_weights(
    position, gene_tss, onek1k_r2_density_data
  )
  eqtlgen_prior_weights_pqtl <- compute_eqtl_tss_dist_prior_weights(
    position, protein_tss, eqtlgen_density_data
  )
  onek1k_r1_prior_weights_pqtl <- compute_eqtl_tss_dist_prior_weights(
    position, protein_tss, onek1k_r1_density_data
  )
  onek1k_r2_prior_weights_pqtl <- compute_eqtl_tss_dist_prior_weights(
    position, protein_tss, onek1k_r1_density_data
  )
  gnocchi_prior_weights <- compute_gnocchi_prior_weights(
    position, chr, gnocchi_data
  )
  abc_score_prior_weights <- compute_abc_prior_weights(
    position, chr, coloc_metadata$gene_name[[i]], abc_score_data
  )
  polyfun_prior_weights <- compute_polyfun_prior_weights(
    position, chr, polyfun_data
  )

  # Uniform priors.

  coloc_unif <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
  )

  # eQTLGen.

  coloc_eqtlgen_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = eqtlgen_prior_weights_eqtl
  )

  coloc_eqtlgen_pqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights2 = eqtlgen_prior_weights_pqtl
  )

  coloc_eqtlgen_pqtl_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = eqtlgen_prior_weights_eqtl,
    prior_weights2 = eqtlgen_prior_weights_pqtl
  )

  # OneK1K round 1.

  coloc_onek1k_r1_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = onek1k_r1_prior_weights_eqtl
  )

  coloc_onek1k_r1_pqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights2 = onek1k_r1_prior_weights_pqtl
  )

  coloc_onek1k_r1_pqtl_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = onek1k_r1_prior_weights_eqtl,
    prior_weights2 = onek1k_r1_prior_weights_pqtl
  )

  # OneK1k round 2.

  coloc_onek1k_r2_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = onek1k_r2_prior_weights_eqtl
  )

  coloc_onek1k_r2_pqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights2 = onek1k_r2_prior_weights_pqtl
  )

  coloc_onek1k_r2_pqtl_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = onek1k_r2_prior_weights_eqtl,
    prior_weights2 = onek1k_r2_prior_weights_pqtl
  )

  # Polyfun priors.

  coloc_polyfun_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = polyfun_prior_weights
  )

  coloc_polyfun_pqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights2 = polyfun_prior_weights
  )

  # Gnocchi priors.

  coloc_gnocchi_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = gnocchi_prior_weights
  )

  coloc_gnocchi_pqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights2 = gnocchi_prior_weights
  )

  # ABC score priors.

  coloc_abc_score_eqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights1 = abc_score_prior_weights
  )

  coloc_abc_score_pqtl <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    prior_weights2 = abc_score_prior_weights
  )

  coloc_results <- bind_cols(
    coloc_to_tibble(coloc_unif, "unif"),
    # eQTLGen.
    coloc_to_tibble(coloc_eqtlgen_eqtl, "eqtlgen-eqtl"),
    coloc_to_tibble(coloc_eqtlgen_pqtl, "eqtlgen-pqtl"),
    coloc_to_tibble(coloc_eqtlgen_pqtl_eqtl, "eqtlgen-pqtl_eqtl"),
    # OneK1K round 1.
    coloc_to_tibble(coloc_onek1k_r1_eqtl, "onek1k_r1-eqtl"),
    coloc_to_tibble(coloc_onek1k_r1_pqtl, "onek1k_r1-pqtl"),
    coloc_to_tibble(coloc_onek1k_r1_pqtl_eqtl, "onek1k_r1-pqtl_eqtl"),
    # OneK1k round 2+.
    coloc_to_tibble(coloc_onek1k_r2_eqtl, "onek1k_r2-eqtl"),
    coloc_to_tibble(coloc_onek1k_r2_pqtl, "onek1k_r2-pqtl"),
    coloc_to_tibble(coloc_onek1k_r2_pqtl_eqtl, "onek1k_r2-pqtl_eqtl"),
    # Gnocchi.
    coloc_to_tibble(coloc_gnocchi_eqtl, "gnocchi-eqtl"),
    coloc_to_tibble(coloc_gnocchi_pqtl, "gnocchi-pqtl"),
    # ABC score.
    coloc_to_tibble(coloc_abc_score_eqtl, "abc_score-eqtl"),
    coloc_to_tibble(coloc_abc_score_pqtl, "abc_score-pqtl"),
    # PolyFun.
    coloc_to_tibble(coloc_polyfun_eqtl, "polyfun-eqtl"),
    coloc_to_tibble(coloc_polyfun_pqtl, "polyfun-pqtl"),
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
