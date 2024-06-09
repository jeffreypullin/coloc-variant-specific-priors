
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(glue)
  library(patchwork)
  library(openxlsx)
  library(stringr)
  library(janitor)
})

source("code/coloc-utils.R")
source("code/plot-utils.R")
devtools::load_all("~/coloc")

coloc_abf_paths <- snakemake@input[["coloc_abf_paths"]]
coloc_susie_paths <- snakemake@input[["coloc_susie_paths"]]

abf_n_colocs_plot_path <- snakemake@output[["abf_n_colocs_plot_path"]]
abf_bootstrap_scatter_plot_path <- snakemake@output[["abf_bootstrap_scatter_plot_path"]]
abf_prob_sig_scatter_plot_path <- snakemake@output[["abf_prob_sig_scatter_plot_path"]]
abf_coloc_results_table_path <- snakemake@output[["abf_coloc_results_table_path"]]
susie_coloc_results_table_path <- snakemake@output[["susie_coloc_results_table_path"]]
susie_prior_effect_plot_path <- snakemake@output[["susie_prior_effect_plot_path"]]
abf_prior_effect_plot_path <- snakemake@output[["abf_prior_effect_plot_path"]]

# Debugging.
config <- yaml::read_yaml("config.yaml")
coloc_susie_paths <- glue(
 "data/output/gwas-eqtl-coloc-susie-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = 6)
)
coloc_abf_paths <- glue(
 "data/output/gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = 6)
)

gwas_eqtl_coloc_abf_data <- tibble(path = coloc_abf_paths) |>
  rowwise() |>
  mutate(split_path = str_split(basename(path), "-")) |>
  mutate(
    gwas_id = split_path[[5]],
    eqtl_id = split_path[[6]],
  ) |>
  select(-split_path) |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

gwas_eqtl_coloc_susie_data <- tibble(path = coloc_susie_paths) |>
  rowwise() |>
  mutate(split_path = strsplit(basename(path), "-")) |>
  mutate(
    gwas_id = split_path[[5]],
    eqtl_id = split_path[[6]],
  ) |>
  select(-split_path) |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

gwas_id_lookup <- c(
  "AUTOIMMUNE" = "Autoimmune",
  "I9_HYPTENS" = "Hypertension",
  "T2D_WIDE" = "Type 2 diabetes"
)

abf_change_coloc_data <- gwas_eqtl_coloc_abf_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior",
    values_to = "pp_h4"
  ) |>
  rename(pp_h4_unif = PP.H4.abf_unif) |>
  separate_wider_delim(prior, "_", names = c("junk", "prior"), too_many = "merge") |>
  select(-junk) |>
  filter(prior %in% c("onek1k_r1", "eqtlgen")) |>
  mutate(prior = prior_method_lookup[prior]) |>
  mutate(type = case_when(
    pp_h4_unif > 0.8 & pp_h4 > 0.8 ~ "no_change",
    pp_h4_unif < 0.8 & pp_h4 < 0.8 ~ "no_change",
    pp_h4_unif < 0.8 & pp_h4 > 0.8 ~ "Newly signifcant",
    pp_h4_unif > 0.8 & pp_h4 < 0.8 ~ "Newly non-significant",
  )) |>
  filter(type != "no_change")

susie_change_coloc_data <- gwas_eqtl_coloc_susie_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior",
    values_to = "pp_h4"
  ) |>
  rename(pp_h4_unif = PP.H4.abf_unif) |>
  separate_wider_delim(prior, "_", names = c("junk", "prior"), too_many = "merge") |>
  select(-junk) |>
  filter(prior %in% c("onek1k_r1", "eqtlgen")) |>
  mutate(prior = prior_method_lookup[prior]) |>
  mutate(type = case_when(
    pp_h4_unif > 0.8 & pp_h4 > 0.8 ~ "no_change",
    pp_h4_unif < 0.8 & pp_h4 < 0.8 ~ "no_change",
    pp_h4_unif < 0.8 & pp_h4 > 0.8 ~ "Newly signifcant",
    pp_h4_unif > 0.8 & pp_h4 < 0.8 ~ "Newly non-significant",
  )) |>
  filter(type != "no_change")

abf_n_colocs_plot <- bind_rows(
  abf_change_coloc_data |>
    count(prior, gwas_id, type) |>
    mutate(method = "abf"),
  susie_change_coloc_data |>
    count(prior, gwas_id, type) |>
    mutate(method = "susie")
) |>
  mutate(n = if_else(type == "Newly non-significant", -n, n)) |>
  ggplot(aes(prior, n, fill = type)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(-5, 5), labels = abs) +
  facet_grid(
    vars(gwas_id),
    vars(method),
    labeller = labeller(gwas_id = gwas_id_lookup),
  ) +
  labs(
    x = "Prior specifcation method",
    y = "Number of unique genes"
  ) +
  theme_jp_vgrid() +
  theme(panel.spacing = unit(1, "lines"))

ggsave(
  abf_n_colocs_plot_path,
  abf_n_colocs_plot,
  width = 12,
  height = 8
)

boostrap_scatter_plot <- gwas_eqtl_coloc_abf_data |>
  rowwise() |>
  mutate(
    q25 = q_pph4_rand_eqtl[[3]],
    q50 = q_pph4_rand_eqtl[[4]],
    q75 = q_pph4_rand_eqtl[[5]],
  ) |>
  ungroup() |>
  ggplot(aes(PP.H4.abf_unif, q50)) +
  geom_point() +
  facet_wrap(~gwas_id) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.005) +
  geom_vline(xintercept = 0.8, linetype = "dotted") +
  theme_jp()

ggsave(
  abf_bootstrap_scatter_plot_path,
  boostrap_scatter_plot,
  width = 12,
  height = 8
)

prob_sig_scatter_plot <- gwas_eqtl_coloc_abf_data |>
  ggplot(aes(PP.H4.abf_unif, prop_sig_pph4_rand_eqtl)) +
  geom_point() +
  geom_vline(xintercept = 0.8) +
  theme_jp()

ggsave(
  abf_prob_sig_scatter_plot_path,
  prob_sig_scatter_plot,
  width = 12,
  height = 8
)

# GWAS-eQTL coloc results

wb <- createWorkbook()

addWorksheet(wb, sheetName = "eQTLGen dist (eQTL)")
gwas_eqtl_coloc_abf_data |>
  mutate(
    sig_eqtlgen_eqtl = PP.H4.abf_eqtlgen > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  arrange(gwas_id) |>
  filter(sig_eqtlgen_eqtl != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_eqtlgen) |>
  writeDataTable(wb, sheet = 1, x = _)

saveWorkbook(
  wb,
  abf_coloc_results_table_path,
  overwrite = TRUE
)

rm(wb)


# coloc.susie()

wb <- createWorkbook()

addWorksheet(wb, sheetName = "eQTLGen dist (eQTL)")
gwas_eqtl_coloc_susie_data |>
  mutate(
    sig_eqtlgen_eqtl = PP.H4.abf_eqtlgen > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  arrange(gwas_id) |>
  filter(sig_eqtlgen_eqtl != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_eqtlgen) |>
  writeDataTable(wb, sheet = 1, x = _)

saveWorkbook(
  wb,
  susie_coloc_results_table_path,
  overwrite = TRUE
)

rm(wb)

# Effect of priors on Pr(H4).

susie_prior_effect_plot <- gwas_eqtl_coloc_susie_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior_type",
    values_to = "pp_h4"
  ) |>
  mutate(diff = pp_h4 - PP.H4.abf_unif) |>
  mutate(prior_type = fct_reorder(factor(prior_type), diff, .fun = var)) |>
  ggplot(aes(prior_type, diff)) +
  geom_boxplot() +
  labs(
    y = "Difference in Pr(H4) values",
    x = "Prior type"
  ) +
  coord_flip() +
  theme_jp_vgrid()

ggsave(
  susie_prior_effect_plot_path,
  susie_prior_effect_plot,
  width = 8,
  height = 8
)

abf_prior_effect_plot <- gwas_eqtl_coloc_abf_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior_type",
    values_to = "pp_h4"
  ) |>
  mutate(diff = pp_h4 - PP.H4.abf_unif) |>
  mutate(prior_type = fct_reorder(factor(prior_type), diff, .fun = var)) |>
  ggplot(aes(prior_type, diff)) +
  geom_boxplot() +
  labs(
    y = "Difference in Pr(H4) values",
    x = "Prior type"
  ) +
  coord_flip() +
  theme_jp_vgrid()

ggsave(
  abf_prior_effect_plot_path,
  abf_prior_effect_plot,
  width = 8,
  height = 8
)
