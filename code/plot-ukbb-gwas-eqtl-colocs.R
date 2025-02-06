
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(glue)
  library(patchwork)
  library(stringr)
  library(janitor)
  library(xtable)
})

source("code/coloc-utils.R")
source("code/plot-utils.R")

# Debugging.
config <- yaml::read_yaml("config.yaml")
coloc_paths <- glue(
  "data/output/ukbb-gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$ukbb_gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = length(config$ukbb_gwas_eqtl_coloc_ids))
)
fm_paths <- glue(
  "data/output/ukbb-gwas-eqtl-finemapping-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$ukbb_gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = length(config$ukbb_gwas_eqtl_coloc_ids))
)

coloc_paths <- snakemake@input[["coloc_paths"]]
fm_paths <- snakemake@input[["fm_paths"]]
ukbb_gwas_eqtl_plot_path <- snakemake@output[["ukbb_gwas_eqtl_plot_path"]]

ukbb_gwas_eqtl_coloc_abf_data <- tibble(path = coloc_paths) |>
  rowwise() |>
  mutate(split_path = str_split(basename(path), "-")) |>
  mutate(
    gwas_id = split_path[[6]],
    eqtl_id = split_path[[7]],
  ) |>
  select(-split_path) |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

change_coloc_data <- ukbb_gwas_eqtl_coloc_abf_data |>
  rename(gene_name = molecular_trait_id) |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior",
    values_to = "pp_h4"
  ) |>
  rename(pp_h4_unif = PP.H4.abf_unif) |>
  separate_wider_delim(prior, "_", names = c("junk", "prior"), too_many = "merge") |>
  select(-junk) |>
  mutate(type = case_when(
    pp_h4_unif > 0.8 & pp_h4 > 0.8 ~ "Unchanged significant",
    pp_h4_unif < 0.8 & pp_h4 < 0.8 ~ "Unchanged not significant",
    pp_h4_unif < 0.8 & pp_h4 > 0.8 ~ "Newly signifcant",
    pp_h4_unif > 0.8 & pp_h4 < 0.8 ~ "Newly non-significant",
  )) |>
  filter(pp_h4_unif > 0.5)

change_coloc_plot <- change_coloc_data |>
  filter(prior %in% c("polyfun_precomputed", "polyfun_trait_specific")) |>
  count(prior, gwas_id, type)  |>
  mutate(prop = n / sum(n), .by = c(prior, gwas_id)) |>
  mutate(gwas_id = ukbb_lookup[gwas_id]) |>
  mutate(prior = ukbb_prior_lookup[prior]) |>
  ggplot(aes(x = prior, y = prop, fill = type)) +
  geom_col() +
  scale_fill_manual(values = c("#CC3311", "#009988", "grey", "grey35")) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(
    vars(gwas_id),
    scales = "free"
  ) +
  coord_flip() +
  labs(
    x = "Prior",
    y = "Percentage of loci",
    fill = "Effect"
  ) +
  theme_jp_vgrid() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Finemapping.

cs_size <- function(p, alpha = 0.95) {
  which.min(cumsum(sort(p, decreasing = TRUE)) < 0.95)
}

ukbb_gwas_eqtl_fm_data <- tibble(path = fm_paths) |>
  rowwise() |>
  mutate(split_path = str_split(basename(path), "-")) |>
  mutate(
    gwas_id = split_path[[5]],
    eqtl_id = split_path[[6]],
  ) |>
  select(-split_path) |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

cs_size_data <- ukbb_gwas_eqtl_fm_data |>
  group_by(region, gwas_id) |>
  filter(molecular_trait_id == first(molecular_trait_id)) |>
  summarise(across(c("unif", "eqtlgen", "polyfun_precomputed", "polyfun_trait_specific"), cs_size), .groups = "drop") |>
  ungroup()

cs_size_plot <- cs_size_data |>
  pivot_longer(
    -c(region, gwas_id),
    names_to = "prior",
    values_to = "cs_size"
  ) |>
  filter(prior != "eqtlgen") |>
  mutate(gwas_id = ukbb_lookup[gwas_id]) |>
  mutate(prior = ukbb_prior_lookup[prior]) |>
  mutate(prior = fct_reorder(factor(prior), cs_size)) |>
  ggplot(aes(prior, cs_size)) +
  geom_boxplot(fill = "grey") +
  coord_flip() +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1, 10, 100),
  ) +
  labs(
    y = "Size of credible set",
    x = "Prior",
  ) +
  facet_wrap(~gwas_id) +
  theme_jp_vgrid()

ukbb_gwas_eqtl_plot <- cs_size_plot / change_coloc_plot +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  ukbb_gwas_eqtl_plot_path,
  ukbb_gwas_eqtl_plot,
  width = 12,
  height = 12
)
