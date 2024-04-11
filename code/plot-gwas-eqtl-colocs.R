
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
#config <- yaml::read_yaml("config.yaml")
#coloc_susie_paths <- glue(
# "data/output/gwas-eqtl-coloc-susie-{gwas_id_eqtl_id}-{chr}.rds",
#  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
#  chr = rep(1:22, each = 6)
#)
#coloc_abf_paths <- glue(
# "data/output/gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds",
#  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
#  chr = rep(1:22, each = 6)
#)

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

abf_n_colocs_plot <- gwas_eqtl_coloc_abf_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id),
    names_to = "prior",
    values_to = "pp_h4"
  ) |>
  mutate(sig_coloc = pp_h4 > 0.8) |>
  summarise(n_sig = sum(sig_coloc), .by = c(prior, gwas_id)) |>
  mutate(is_unif = prior == "PP.H4.abf_unif") |>
  mutate(prior = fct_reorder(factor(prior), n_sig, mean)) |>
  ggplot(aes(prior, n_sig, fill = is_unif)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("#66CCEE", "#EE6677")) +
  labs(
    y = "Number of signifcant (Pr(H4) > 0.8) colocalisations",
    x = "Prior type"
  ) +
  facet_wrap(~gwas_id, scales = "free_x") +
  theme_jp_vgrid() +
  theme(legend.position = "none")

ggsave(
  abf_n_colocs_plot_path,
  abf_n_colocs_plot,
  width = 12,
  height = 8
)

boostrap_scatter_plot <- gwas_eqtl_coloc_abf_data |>
  rowwise() |>
  mutate(
    q25 = q_pph4_rand_gwas[[3]],
    q50 = q_pph4_rand_gwas[[4]],
    q75 = q_pph4_rand_gwas[[5]],
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
  ggplot(aes(PP.H4.abf_unif, prop_sig_pph4_rand_gwas)) +
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
    sig_eqtlgen_eqtl = PP.H4.abf_eqtl_tss_eqtlgen > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  arrange(gwas_id) |>
  filter(sig_eqtlgen_eqtl != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_eqtl_tss_eqtlgen) |>
  writeDataTable(wb, sheet = 1, x = _)


addWorksheet(wb, sheetName = "Gnocchi (GWAS)")
gwas_eqtl_coloc_abf_data |>
  mutate(
    sig_gnocchi_gwas = PP.H4.abf_gnocchi_gwas > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  arrange(gwas_id) |>
  filter(sig_gnocchi_gwas != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_gnocchi_gwas) |>
  writeDataTable(wb, sheet = 2, x = _)

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
    sig_eqtlgen_eqtl = PP.H4.abf_eqtl_tss_eqtlgen > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  arrange(gwas_id) |>
  filter(sig_eqtlgen_eqtl != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_eqtl_tss_eqtlgen) |>
  writeDataTable(wb, sheet = 1, x = _)

addWorksheet(wb, sheetName = "Gnocchi (GWAS)")
gwas_eqtl_coloc_susie_data |>
  mutate(
    sig_gnocchi_gwas = PP.H4.abf_gnocchi_gwas > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  arrange(gwas_id) |>
  filter(sig_gnocchi_gwas != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_gnocchi_gwas) |>
  writeDataTable(wb, sheet = 2, x = _)

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
