# The code in this file is heavily based on:
# https://github.com/ralf-tambets/coloc/blob/main/bin/coloc_v3_pqtl.R

source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(seqminer)
  library(dplyr)
  devtools::load_all("~/coloc")
})

eqtl_file <- snakemake@input[["eqtl_data_file"]]
pqtl_file <- snakemake@input[["pqtl_data_file"]]
eqtl_metadata_file <- snakemake@input[["eqtl_metadata_file"]]
pqtl_metadata_file <- snakemake@input[["pqtl_metadata_file"]]
chr <- as.numeric(snakemake@wildcards[["chr"]])

eqtl_metadata <- read_tsv(eqtl_metadata_file, show_col_types = FALSE)
pqtl_metadata <- read_tsv(pqtl_metadata_file, show_col_types = FALSE)

width <- 1e5
coloc_metadata <- eqtl_metadata |>
  mutate(tss = if_else(strand == 1, gene_start, gene_end)) |>
  rowwise() |>
  mutate(
    start_pos = max(tss - width, 1),
    end_pos = tss + width
  ) |>
  mutate(region = paste0(chromosome, ":", start_pos, "-", end_pos)) |>
  select(gene_name, region, gene_id, chromosome, start_pos, end_pos, tss) |>
  ungroup() |>
  left_join(
    pqtl_metadata |>
      select(phenotype_id, chromosome, phenotype_pos, pqtl_gene_id = gene_id) |>
      mutate(chromosome = as.character(chromosome)),
    by = join_by(chromosome == chromosome,
                 start_pos <= phenotype_pos,
                 end_pos >= phenotype_pos)
  ) |>
  select(gene_name, gene_id, phenotype_id, region, chromosome,
         start_pos, end_pos, tss, pqtl_gene_id, phenotype_pos)

all_eqtl_data <- tabix.read.table(eqtl_file, paste0(chr, ":1-2147483647")) |>
  as_tibble()
all_pqtl_data <- tabix.read.table(pqtl_file, paste0(chr, ":1-2147483647")) |>
  as_tibble()

coloc_metadata <- coloc_metadata |>
  filter(chromosome == chr) |>
  filter(!is.na(phenotype_id)) |>
  filter(!is.na(gene_id))

n <- nrow(coloc_metadata)

results <- list()
for (i in 1:n) {

  eqtl_data <- all_eqtl_data |>
    filter(V1 == coloc_metadata$gene_id[[i]]) |>
    filter(V2 == coloc_metadata$chromosome[[i]]) |>
    filter(
      V3 >= coloc_metadata$start_pos[[i]] &
      V3 <= coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble() |>
    setNames(c("molecular_trait_id", "chromosome", "position", "ref",
               "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
               "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id",
               "median_tpm", "rsid"))

  pqtl_data <- all_pqtl_data |>
    filter(V1 == coloc_metadata$phenotype_id[[i]]) |>
    filter(V2 == coloc_metadata$chromosome[[i]]) |>
    filter(
      V3 >= coloc_metadata$start_pos[[i]] &
      V3 <= coloc_metadata$end_pos[[i]]
    ) |>
    as_tibble() |>
    setNames(c("molecular_trait_id", "chromosome", "position", "ref",
               "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
               "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id",
               "median_tpm", "rsid"))

  if (nrow(eqtl_data) == 0 || nrow(pqtl_data) == 0) {
    next
  }

  pqtl_data <- pqtl_data |>
    select(-rsid) |>
    distinct() |>
    mutate(id = paste(chromosome, position, sep = ":")) |>
    group_by(id) |>
    mutate(row_count = n()) |>
    ungroup() |>
    filter(row_count == 1) |>
    filter(!is.nan(se)) |>
    filter(!is.na(se)) |>
    select(molecular_trait_id, variant, maf, beta, se, an) |>
    mutate(
      maf = as.numeric(maf),
      beta = as.numeric(beta),
      se = as.numeric(se),
      an = as.numeric(an)
    )

  eqtl_data <- eqtl_data |>
    select(-rsid) |>
    distinct() |>
    mutate(id = paste(chromosome, position, sep = ":")) |>
    group_by(id) |>
    mutate(row_count = n()) |>
    ungroup() |>
    filter(row_count == 1) |>
    filter(!is.nan(se)) |>
    filter(!is.na(se)) |>
    select(molecular_trait_id, variant, maf, beta, se, an, position) |>
    mutate(
      maf = as.numeric(maf),
      beta = as.numeric(beta),
      se = as.numeric(se),
      an = as.numeric(an)
    )

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

  varbeta <- eqtl_data$se^2
  n <- eqtl_data$an / 2
  maf <- eqtl_data$maf

  is_oneover_na <- all(is.na(1 / varbeta))
  is_nvx_na <- all(is.na(2 * n * maf * (1 - maf)))

  if (is_oneover_na || is_nvx_na) {
    next
  }

  eqtl_dataset <- list(
    varbeta = varbeta,
    N = n,
    MAF = maf,
    type = "quant",
    beta = eqtl_data$beta,
    snp = eqtl_data$variant,
    position = eqtl_data$position
  )

  varbeta <- pqtl_data$se^2
  n <- pqtl_data$an / 2
  maf <- pqtl_data$maf

  is_oneover_na <- all(is.na(1 / varbeta))
  is_nvx_na <- all(is.na(2 * n * maf * (1 - maf)))

  if (is_oneover_na || is_nvx_na) {
    next
  }

  pqtl_dataset <- list(
    varbeta = varbeta,
    N = n,
    MAF = maf,
    type = "quant",
    beta = pqtl_data$beta,
    snp = pqtl_data$variant
  )

  gnocchi_data <- read_tsv(
    "/home/jp2045/coloc-estimated-eqtl-priors/data/gnocchi-windows.bed",
    col_names = FALSE, show_col_types = FALSE
  )
  colnames(gnocchi_data) <- c("chromosome", "start_pos", "end_pos", "score")

  coloc_out_unif <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
  )

  coloc_out_eqtl_tss <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    tss1 = coloc_metadata$tss[[i]]
  )
  
  coloc_out_pqtl_tss <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    tss2 = coloc_metadata$phenotype_pos[[i]]
  )

  coloc_out_score <- coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = pqtl_dataset,
    score_data1 = gnocchi_data,
    chrom = chr
  )

  coloc_summary_to_tibble <- function(coloc_out) {
    t(as.data.frame(coloc_out$summary))
  }
  coloc_res_unif <- coloc_summary_to_tibble(coloc_out_unif)
  colnames(coloc_res_unif) <- paste0(colnames(coloc_res_unif), "_unif")

  coloc_res_eqtl_tss <- coloc_summary_to_tibble(coloc_out_eqtl_tss)
  colnames(coloc_res_eqtl_tss) <- paste0(colnames(coloc_res_eqtl_tss), "_eqtl_tss")

  coloc_res_pqtl_tss <- coloc_summary_to_tibble(coloc_out_pqtl_tss)
  colnames(coloc_res_pqtl_tss) <- paste0(colnames(coloc_res_pqtl_tss), "_pqtl_tss")

  coloc_res_score <- coloc_summary_to_tibble(coloc_out_score)
  colnames(coloc_res_score) <- paste0(colnames(coloc_res_score), "_score")

  coloc_results <- bind_cols(
    coloc_res_unif,
    coloc_res_eqtl_tss,
    coloc_res_pqtl_tss,
    coloc_res_score,
    tibble(
     phenotype_id = coloc_metadata$phenotype_id[[i]],
     chromosome = coloc_metadata$chromosome[[i]],
     gene_id = coloc_metadata$gene_id[[i]],
     region = coloc_metadata$region[[i]],
     gene_name = coloc_metadata$gene_name[[i]]
  ))
  results[[i]] <- coloc_results
}

write_tsv(
  bind_rows(!!!results),
  snakemake@output[["result_file"]]
)
