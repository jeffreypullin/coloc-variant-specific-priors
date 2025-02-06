
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
  library(forcats)
  library(patchwork)
})

source("code/coloc-utils.R")
source("code/plot-utils.R")
source("code/prior-probabilities-funs.R")

sn <- snakemake
gnocchi_data <- read_tsv(sn@input[["gnocchi_data_path"]], show_col_types = FALSE)
abc_score_data <- read_tsv(sn@input[["abc_score_data_path"]], show_col_types = FALSE)
density_data_round_1 <- read_rds(sn@input[["density_data_r1_path"]])
density_data_round_2 <- read_rds(sn@input[["density_data_r2_path"]])
eqtlgen_density_data <- read_rds(sn@input[["eqtlgen_density_data_path"]])
snp_var_data_1_7 <- read_parquet(sn@input[["snp_var_data_1_7_path"]])

plot_path <- snakemake@output[["all_priors_plot_path"]]
polyfun_priors_plot_path <- snakemake@output[["polyfun_priors_plot_path"]]

chr <- 1
tss <- 107965180
width <- 5e5
gene_name <- "VAV3"
region <- paste0(chr, ":", tss - width, "-",  tss + width)

position <- tabix.read.table(
  sn@input[["autoimmune_gwas_path"]],
  region
) |>
  as_tibble() |>
  rename(
    rsid = rsids,
    chromosome = chrom,
    position = pos,
    se = sebeta,
    maf = af_alt
  ) |>
  mutate(
    # Sample size is not needed for prior calculation.
    sample_size = 100,
    molecular_trait_id = "t1d",
    variant = paste0("chr", chromosome, "_", position, "_", ref, "_", alt)
  ) |>
  select(-c(sample_size, ref, alt)) |>
  prepare_coloc_dataset() |>
  pull(position)

eqtl_prior_weights_eqtlgen <- compute_eqtl_tss_dist_prior_weights(
  position, tss, eqtlgen_density_data
)
eqtl_prior_weights <- compute_eqtl_tss_dist_prior_weights(
  position, tss, density_data_round_1
)
eqtl_prior_weights_round_2 <- compute_eqtl_tss_dist_prior_weights(
  position, tss, density_data_round_2
)
gnocchi_prior_weights <- compute_gnocchi_prior_weights(
  position, chr, gnocchi_data
)
abc_score_prior_weights <- compute_abc_prior_weights(
  position, chr, gene_name, abc_score_data
)
polyfun_prior_weights <- compute_polyfun_prior_weights(
  position, chr, snp_var_data_1_7
)

prior_type_lookup <- c(
  eqtlgen = "eQTLGen",
  onek1k_round_1 = "OneK1K (R1)",
  onek1k_round_2 = "OneK1K (R2+)",
  abc_score_all = "ABC Score",
  gnocchi = "Gnocchi",
  polyfun = "PolyFun"
)
prior_type_lookup <- factor(prior_type_lookup, levels = prior_type_lookup)

priors_plot <- tibble(
  position = position,
  eqtlgen = eqtl_prior_weights_eqtlgen,
  onek1k_round_1 = eqtl_prior_weights,
  onek1k_round_2 = eqtl_prior_weights_round_2,
  abc_score_all = abc_score_prior_weights,
  gnocchi = gnocchi_prior_weights,
  polyfun = polyfun_prior_weights
) |>
  mutate(position = position / 1e6) |>
  pivot_longer(
    cols = -position,
    names_to = "prior_type",
    values_to = "weight"
  ) |>
  mutate(prior_type = prior_type_lookup[prior_type]) |>
  ggplot(aes(position, weight)) +
  geom_point(col = "darkgrey") +
  geom_vline(xintercept = tss / 1e6, linetype = "dashed", color = "black") +
  labs(
    x = "Position (Mb)",
    y = "Prior probability",
  ) +
  scale_x_continuous(breaks = seq(107.5, 108.3, by = 0.4)) +
  facet_wrap(~prior_type) +
  theme_jp() +
  theme(axis.text = element_text(size = 14))

ggsave(
  plot_path,
  plot = priors_plot,
  width = 8,
  height = 6
)

# Plot Polyfun UKBB priors.

ukbb_244_1 <- read_tsv("data/polyfun-output/244/hg38-polyfun-out.1.snpvar_constrained.gz", show_col_types = FALSE)
ukbb_250_2_1 <- read_tsv("data/polyfun-output/250.2/hg38-polyfun-out.1.snpvar_constrained.gz", show_col_types = FALSE)
ukbb_401_1 <- read_tsv("data/polyfun-output/401/hg38-polyfun-out.1.snpvar_constrained.gz",  show_col_types = FALSE)

normalise <- function(x) {
  x / sum(x)
}

ukbb_polyfun_weights_plot <- bind_rows(
  ukbb_244_1 |>
    filter(BP > tss - width, BP < tss + width) |>
    mutate(trait = "Hypothyroidism") |>
    select(pos = BP, weight = "SNPVAR", trait) |>
    mutate(weight = normalise(weight)),
  ukbb_250_2_1 |>
    filter(BP > tss - width, BP < tss + width) |>
    mutate(trait = "Type 2 diabetes") |>
    select(pos = BP, weight = "SNPVAR", trait) |>
    mutate(weight = normalise(weight)),
  ukbb_401_1 |>
    filter(BP > tss - width, BP < tss + width) |>
    mutate(trait = "Hypertension") |>
    select(pos = BP, weight = "SNPVAR", trait) |>
    mutate(weight = normalise(weight)),
  tibble(
    pos = position,
    weight = polyfun_prior_weights,
    trait = "PolyFun pre-computed"
  )
) |>
  mutate(pos = pos / 1e6) |>
  ggplot(aes(pos, weight)) +
  geom_point(col = "darkgrey") +
  geom_vline(xintercept = tss / 1e6, linetype = "dashed", color = "black") +
  labs(
    x = "Position (Mb)",
    y = "Prior probability",
  ) +
  scale_x_continuous(breaks = seq(107.5, 108.3, by = 0.4)) +
  facet_wrap(~trait) +
  theme_jp() +
  theme(axis.text = element_text(size = 14))

ggsave(
  polyfun_priors_plot_path,
  plot = ukbb_polyfun_weights_plot,
  width = 14,
  height = 8
)