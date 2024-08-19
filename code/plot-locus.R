
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
  library(EnsDb.Hsapiens.v86)
  library(locuszoomr)
  library(patchwork)
  library(glue)
  library(stringr)
})

source("code/coloc-utils.R")
source("code/plot-utils.R")
source("code/prior-probabilities-funs.R")

# For debugging.
config <- yaml::read_yaml("config.yaml")
coloc_abf_paths <- glue(
  "data/output/gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = length(config$gwas_eqtl_coloc_ids))
)

coloc_abf_paths <- snakemake@input[["coloc_abf_paths"]]
autoimmune_gwas_path <- snakemake@input[["autoimmune_gwas_path"]]
gtex_thyroid_path <- snakemake@input[["gtex_thyroid_path"]]

loci_range <- "9:123257606-125257606"
gwas_data <- tabix.read.table(autoimmune_gwas_path, loci_range)
eqtl_data <- tabix.read.table(gtex_thyroid_path, loci_range) |>
  setNames(eqtl_catalouge_colnames)

nek6_eqtl_data <- eqtl_data |>
  dplyr::filter(molecular_trait_id == "ENSG00000119408")

psmb7_eqtl_data <- eqtl_data |>
  dplyr::filter(molecular_trait_id == "ENSG00000136930")

gwas_loc <- locus(
  data = gwas_data,
  gene = "NEK6",
  flank = 1.1e5,
  ens_db = "EnsDb.Hsapiens.v86"
)

nek6_eqtl_loc <- locus(
  data = nek6_eqtl_data,
  gene = "NEK6",
  flank = 1.1e5,
  ens_db = "EnsDb.Hsapiens.v86"
)

psmb7_eqtl_loc <- locus(
  data = psmb7_eqtl_data,
  gene = "NEK6",
  flank = 1.1e5,
  ens_db = "EnsDb.Hsapiens.v86"
)

nek6_prior_data <- nek6_eqtl_data
nek6_prior_data$pvalue <- 0.1^(compute_eqtl_tss_dist_prior_weights(
  nek6_eqtl_data$pos,
  124257606,
  read_rds("output/densities/eqtlgen.rds")
))
nek6_prior_loc <- locus(
  data = nek6_prior_data,
  gene = "NEK6",
  flank = 1.1e5,
  ens_db = "EnsDb.Hsapiens.v86"
)

psmb7_prior_data <- psmb7_eqtl_data
psmb7_prior_data$pvalue <- 0.1^(compute_eqtl_tss_dist_prior_weights(
  psmb7_eqtl_data$pos,
  124415444,
  read_rds("output/densities/eqtlgen.rds")
))
psmb7_prior_loc <- locus(
  data = psmb7_prior_data,
  gene = "NEK6",
  flank = 1.1e5,
  ens_db = "EnsDb.Hsapiens.v86"
)

gwas_eqtl_coloc_abf_data <- tibble(path = coloc_abf_paths) |>
  rowwise() |>
  mutate(split_path = str_split(basename(path), "-")) |>
  mutate(
    gwas_id = split_path[[5]],
    eqtl_id = split_path[[6]],
  ) |>
  dplyr::select(-split_path) |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

pp_h4_plot <- gwas_eqtl_coloc_abf_data |>
  dplyr::filter(
    gene_name %in% c("NEK6", "PSMB7"),
    gwas_id == "AUTOIMMUNE",
    eqtl_id == "QTD000341"
 ) |>
 dplyr::select(PP.H4.abf_unif, PP.H4.abf_eqtlgen, gene_name) |>
 pivot_longer(cols = -gene_name) |>
 mutate(name = if_else(name == "PP.H4.abf_unif", "Uniform", "eQTLGen")) |>
 mutate(name = factor(name, levels = c("Uniform", "eQTLGen"))) |>
 ggplot(aes(gene_name, value)) +
 geom_col(fill = "grey") +
 geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red") +
 facet_wrap(~name) +
 labs(
    x = "Gene",
    y = "Pr(H4)"
 ) +
 theme_jp()

nek6_psmb7_example_plot <- gg_scatter(gwas_loc, xticks = FALSE) +
  ggtitle("FinnGen autoimmune disease GWAS") +
  gg_scatter(nek6_eqtl_loc, xticks = FALSE) +
  ggtitle("NEK6 eQTL (GTEx thyroid)") +
  gg_scatter(psmb7_eqtl_loc, xticks = FALSE) +
  ggtitle("PSMB7 eQTL (GTEx thyroid)") +
  gg_genetracks(gwas_loc, filter_gene_name = c("NEK6", "PSMB7")) +
  pp_h4_plot +
  gg_scatter(nek6_prior_loc, yzero = FALSE, ylab = "Prior probability", index_snp = NULL, xticks = FALSE) +
  gg_scatter(psmb7_prior_loc, yzero = FALSE, ylab = "Prior probability", index_snp = NULL, xticks = FALSE) +
  gg_genetracks(gwas_loc, filter_gene_name = c("NEK6", "PSMB7")) +
  plot_layout(ncol = 2, nrow = 4, byrow = FALSE, heights = c(1, 1, 1, 0.3)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "", "d", "e", "f")))

ggsave(
  snakemake@output[["nek6_psmb7_plot_path"]],
  nek6_psmb7_example_plot,
  width = 8,
  height = 10
)
