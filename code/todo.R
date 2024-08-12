
devtools::load_all("~/coloc")

suppressPackageStartupMessages({
  library(ggplot2)
  library(arrow)
  library(readr)
  library(tidyr)
  library(seqminer)
  library(janitor)
  library(purrr)
  library(EnsDb.Hsapiens.v86)
  library(locuszoomr)
  library(patchwork)
  library(glue)
  library(stringr)
  library(dplyr)
  library(tibble)
})

devtools::load_all("~/coloc")

source("code/coloc-utils.R")
source("code/plot-utils.R")
source("code/prior-probabilities-funs.R")

autoimmune_gwas_path <- "data/finngen/AUTOIMMUNE.gz"
gtex_thyroid_path <- "data/eqtl-catalogue/sumstats/QTD000341.cc.tsv.gz"

loci_range <- "9:123915444-124915444"
gwas_data <- tabix.read.table(autoimmune_gwas_path, loci_range)
eqtl_data <- tabix.read.table(gtex_thyroid_path, loci_range) |>
  setNames(eqtl_catalouge_colnames)

nek6_eqtl_data <- eqtl_data |>
  filter(molecular_trait_id == "ENSG00000119408")

psmb7_eqtl_data <- eqtl_data |>
  filter(molecular_trait_id == "ENSG00000136930")

gwas_data <- tabix.read.table(autoimmune_gwas_path, loci_range) |>
  rename(
    rsid = rsids,
    chromosome = chrom,
    position = pos,
    beta = beta,
    se = sebeta,
    maf = af_alt
  ) |>
  mutate(
    N = 105677 + 306504,
    variant = paste0("chr", chromosome, "_", position, "_", ref, "_", alt)
  ) |>
  as_tibble()

gwas_data <- prepare_coloc_dataset(gwas_data)
eqtl_data <- prepare_coloc_dataset(psmb7_eqtl_data)

eqtl_data_filtered <- eqtl_data |>
  filter(variant %in% gwas_data$variant)
gwas_data_filtered <- gwas_data |>
  filter(variant %in% eqtl_data$variant)

eqtl_dataset <- list(
  varbeta = eqtl_data$se^2,
  N = eqtl_data$an / 2,
  MAF = eqtl_data$maf,
  type = "quant",
  beta = eqtl_data$beta,
  snp = eqtl_data$variant,
  position = eqtl_data$position
)

gwas_dataset <- list(
 varbeta = gwas_data$se^2,
  N = gwas_data$N,
  MAF = gwas_data$maf,
  type = "cc",
  beta = gwas_data$beta,
  snp = gwas_data$variant
)

eqtl_dataset_filtered <- list(
  varbeta = eqtl_data_filtered$se^2,
  N = eqtl_data_filtered$an / 2,
  MAF = eqtl_data_filtered$maf,
  type = "quant",
  beta = eqtl_data_filtered$beta,
  snp = eqtl_data_filtered$variant,
  position = eqtl_data_filtered$position
)

gwas_dataset_filtered <- list(
  varbeta = gwas_data_filtered$se^2,
  N = gwas_data_filtered$N,
  MAF = gwas_data_filtered$maf,
  type = "cc",
  beta = gwas_data_filtered$beta,
  snp = gwas_data_filtered$variant
)

w_filtered <- compute_eqtl_tss_dist_prior_weights(
  eqtl_data_filtered$position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

w <- compute_eqtl_tss_dist_prior_weights(
  eqtl_data$position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

coloc.abf(gwas_dataset, eqtl_dataset)
coloc.abf(gwas_dataset_filtered, eqtl_dataset_filtered)

coloc.abf(gwas_dataset, eqtl_dataset, prior_weights2 = w)
coloc.abf(gwas_dataset_filtered, eqtl_dataset_filtered, prior_weights2 = w_filtered)

gwas_cs <- read_tsv("data/finngen/AUTOIMMUNE.SUSIE.cred.bgz", show_col_types = FALSE)
gwas_lbf <- tabix.read.table("data/finngen/AUTOIMMUNE.SUSIE.snp.bgz", "chr9:1-2147483647") |>
  as_tibble() |>
  setNames(c("trait", "region", "v", "rsid", "chromosome", "position", "allele1", "allele2",
             "maf", "beta", "se", "p", "mean", "sd", "prob", "cs", "cs_specific_prob",
             "low_purity", "lead_r2", "mean_99", "sd_99", "prob_99", "cs_99", "cs_specific_prob_99",
             "low_purity_99", "lead_r2_99", paste0("alpha", 1:10),  paste0("mean", 1:10),
             paste0("sd", 1:10), paste0("lbf_variable", 1:10)))

eqtl_cs <- read_tsv("data/eqtl-catalogue/susie/QTD000341.credible_sets.tsv.gz", show_col_types = FALSE)
eqtl_lbf <- tabix.read.table("data/eqtl-catalogue/susie/QTD000341.lbf_variable.txt.gz", "9:1-2147483647") |>
  as_tibble() |>
  setNames(c("molecular_trait_id", "region", "variant", "chromosome", "position", paste0("lbf_variable", 1:10)))

signal_lbf <- gwas_lbf |>
  filter(region == "chr9:122775517-125775517") |>
  select(variant = rsid, lbf_variable1)

signal_mat <- signal_lbf |>
  column_to_rownames("variant") |>
  t()

psmb7_lbf <- eqtl_lbf |>
  filter(molecular_trait_id == "ENSG00000136930") |>
  select(variant, lbf_variable1, position)

psmb7_mat <- psmb7_lbf |>
  select(-position) |>
  column_to_rownames("variant") |>
  t()

w_bf <- compute_eqtl_tss_dist_prior_weights(
  psmb7_lbf$position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

coloc.bf_bf(psmb7_mat, signal_mat)
coloc.bf_bf(psmb7_mat, signal_mat, prior_weights1 = w_bf)

psmb7_lbf_filt <- psmb7_lbf |>
  filter(variant %in% signal_lbf$variant)

signal_lbf_filt <- signal_lbf |>
  filter(variant %in% psmb7_lbf$variant)

psmb7_mat_filt <- psmb7_lbf_filt |>
  select(-position) |>
  column_to_rownames("variant") |>
  t()

signal_mat_filt <- signal_lbf_filt |>
  column_to_rownames("variant") |>
  t()

position <- psmb7_lbf$position
position_filt <- psmb7_lbf |>
  filter(variant %in% psmb7_lbf$variant) |>
  filter(variant %in% signal_lbf$variant) |>
  pull(position)

isnps <- intersect(colnames(psmb7_mat), colnames(signal_mat))
ind  <- match(isnps, colnames(psmb7_mat))

head(position_filt, n = 40)
head(position[ind], n = 40)

plot(position_filt, position[ind])

devtools::load_all("~/coloc")

filt_ind <- which(position_filt %in% position)

w_bf <- compute_eqtl_tss_dist_prior_weights(
  position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

w_bf_filt <- compute_eqtl_tss_dist_prior_weights(
  position[ind], 124415444, read_rds("output/densities/eqtlgen.rds")
)

devtools::load_all("~/coloc")

coloc.bf_bf(psmb7_mat, signal_mat, prior_weights1 = w_bf)$summary
coloc.bf_bf(psmb7_mat_filt, signal_mat_filt, prior_weights1 = w_bf_filt)$summary

isnps_test <- intersect(colnames(psmb7_mat), colnames(signal_mat))
filt_w_bf_test <- w_bf[match(isnps, colnames(psmb7_mat))]
summary(filt_w_bf / sum(filt_w_bf))

plot(filt_w_bf / sum(filt_w_bf), w_bf_filt)

a <- filt_w_bf / sum(filt_w_bf)
b <- w_bf_filt
plot(a / b)



a <- w_bf[filt_ind] / sum(w_bf[filt_ind])
b <- w_bf_filt

plot(position, w_bf)
plot(position_filt, w_bf_filt)

plot(w_bf[filt_ind] / w_bf_filt)
plot(a)
plot(b)

polyfun_data <- arrow::read_parquet("data/snpvar_meta.chr8_22.parquet")
w_p <- compute_polyfun_prior_weights(
  psmb7_lbf$position, 9, polyfun_data
)

w_p_filt <- compute_polyfun_prior_weights(
  position, 9, polyfun_data
)

a <- w_p[filt_ind] / sum(w_p[filt_ind])
b <- w_p_filt

plot(a / b)

compute_simple_weights <- function(pos) {
  dist <- abs(pos - mean(pos))
  dist / sum(dist)
}

w_s <- compute_simple_weights(psmb7_lbf$position)
w_s_filt <- compute_simple_weights(position)

a <- w_s[filt_ind] / sum(w_s[filt_ind])
b <- w_s_filt

summary(a / b)

devtools::load_all("~/coloc")

w_bf <- compute_eqtl_tss_dist_prior_weights(
  psmb7_lbf$position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

w_bf_filt <- compute_eqtl_tss_dist_prior_weights(
  position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

filt_ind <- which(colnames(psmb7_mat_filt) %in% colnames(psmb7_mat))

coloc.bf_bf(psmb7_mat, signal_mat, prior_weights1 = w_bf)$summary
coloc.bf_bf(psmb7_mat_filt, signal_mat_filt, prior_weights1 = w_bf[filt_ind])$summary


a <- w_bf[filt_ind] / sum(w_bf[filt_ind])
b <- w_bf_filt

w <- sample(1:100, 10)
w





filt_index <- which(position %in% psmb7_lbf$position)

head(w_bf[filt_index] / sum(w_bf[filt_index]))
head(w_bf_filt)

w_bf_filt / w_bf[filt_index]

density_data <- read_rds("output/densities/eqtlgen.rds"

rel <- pos - tss
x <- density_data$x
y <- density_data$y
closest <- numeric(length(pos))
for (i in seq_along(pos)) {
  z <- which(abs(x - rel[[i]]) == min(abs(x - rel[[i]])))
  if (length(z) > 1) {
    closest[[i]] <- z[[1]]
  } else {
    closest[[i]] <- z
  }
}
out <- y[closest] / sum(y[closest])

out


position

coloc.bf_bf(psmb7_mat, signal_mat, prior_weights1 = w_bf)
coloc.bf_bf(psmb7_mat_filt, signal_mat_filt, prior_weights1 = w_bf_filt)


all.equal(w_bf[filt_index] / sum(w_bf[filt_index]), w_bf_filt)

head(w_bf[filt_index])
head(w_bf_filt)

plot(w_bf[filt_index])
plot(w_bf_filt)

plot(w_bf)
plot(w_bf_filt)



filt_index <- which(psmb7_lbf_filt$variant %in% psmb7_lbf$variant)

a <- compute_eqtl_tss_dist_prior_weights(
  psmb7_lbf$position, 124415444, read_rds("output/densities/eqtlgen.rds")
)

b <- compute_eqtl_tss_dist_prior_weights(
  psmb7_lbf$position[filt_index], 124415444, read_rds("output/densities/eqtlgen.rds")
)

head(a[filt_index] / sum(a[filt_index]))
head(b)




density_data <- read_rds("output/densities/eqtlgen.rds")
rel <- psmb7_lbf$position - 124415444
x <- density_data$x
y <- density_data$y
z <- which(abs(x - rel[[i]]) == min(abs(x - rel[[i]])))
if (length(z) > 1) {
  closest[[i]] <- z[[1]]
} else {
  closest[[i]] <- z
}


summary(w_bf)
summary(w_bf_filt)

length(w_bf)
length(w_bf_filt)

plot(w_bf)
plot(w_bf_filt)

dim(psmb7_mat_filt)


test <- coloc:::process.dataset(eqtl_dataset, "1") |>
  as_tibble()

hist(test$lABF.1)
hist(psmb7_lbf$lbf_variable1)

gwas_eqtl_coloc_susie_data |>
  filter(gwas_id == "AUTOIMMUNE" & gene_name == "PSMB7") |>
  select(PP.H4.abf_unif, PP.H4.abf_eqtlgen, gwas_region, eqtl_region)
