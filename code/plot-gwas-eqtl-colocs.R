
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
  library(latex2exp)
  library(openxlsx)
})

source("code/coloc-utils.R")
source("code/plot-utils.R")
devtools::load_all("~/coloc")

# Debugging.
config <- yaml::read_yaml("config.yaml")
coloc_susie_paths <- glue(
  "data/output/gwas-eqtl-coloc-susie-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = length(config$gwas_eqtl_coloc_ids))
)
coloc_abf_paths <- glue(
  "data/output/gwas-eqtl-coloc-abf-{gwas_id_eqtl_id}-{chr}.rds",
  gwas_id_eqtl_id = rep(config$gwas_eqtl_coloc_ids, 22),
  chr = rep(1:22, each = length(config$gwas_eqtl_coloc_ids))
)

coloc_abf_paths <- snakemake@input[["coloc_abf_paths"]]
coloc_susie_paths <- snakemake@input[["coloc_susie_paths"]]

overall_impact_plot_path <- snakemake@output[["overall_impact_plot_path"]]
coloc_results_excel_path <- snakemake@output[["coloc_results_excel_path"]]

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
  "AUTOIMMUNE" = "Autoimmune\ndisease",
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
    pp_h4_unif > 0.8 & pp_h4 > 0.8 ~ "Unchanged significant",
    pp_h4_unif < 0.8 & pp_h4 < 0.8 ~ "Unchanged not significant",
    pp_h4_unif < 0.8 & pp_h4 > 0.8 ~ "Newly signifcant",
    pp_h4_unif > 0.8 & pp_h4 < 0.8 ~ "Newly non-significant",
  )) |>
  filter(pp_h4_unif > 0.5)

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
    pp_h4_unif > 0.8 & pp_h4 > 0.8 ~ "Unchanged significant",
    pp_h4_unif < 0.8 & pp_h4 < 0.8 ~ "Unchanged not significant",
    pp_h4_unif < 0.8 & pp_h4 > 0.8 ~ "Newly signifcant",
    pp_h4_unif > 0.8 & pp_h4 < 0.8 ~ "Newly non-significant",
  )) |>
  filter(pp_h4_unif > 0.5)

n_colocs_plot <- bind_rows(
  abf_change_coloc_data |>
    count(prior, gwas_id, type) |>
    mutate(method = "coloc-single"),
  susie_change_coloc_data |>
    count(prior, gwas_id, type) |>
    mutate(method = "coloc-susie")
) |>
  mutate(gwas_id = gwas_id_lookup[gwas_id]) |>
  ggplot(aes(x = gwas_id, y = n, fill = type)) +
  geom_col() +
  facet_grid(
    vars(method),
    vars(prior),
    scales = "free"
  ) +
  labs(
    x = "Trait",
    y = "Number of loci",
    fill = "Effect"
  ) +
  theme_jp() +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

abf_eqtlgen_scatter_plot <- gwas_eqtl_coloc_abf_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior_type",
    values_to = "pp_h4"
  ) |>
  filter(gwas_id == "AUTOIMMUNE", prior_type == "PP.H4.abf_eqtlgen") |>
  rename(pp_h4_unif = PP.H4.abf_unif) |>
  ggplot(aes(pp_h4_unif, pp_h4)) +
  geom_point() +
  geom_vline(xintercept = 0.8, linetype = "dashed", col = "red") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  labs(
    x = TeX("Uniform prior $\\  \\Pr(H_4)$"),
    y = TeX("eQTLGen prior $\\  \\Pr(H_4)$")
  ) +
  theme_jp()

susie_eqtlgen_scatter_plot <- gwas_eqtl_coloc_susie_data |>
  select(gene_name, gwas_id, starts_with("PP.H4.abf")) |>
  pivot_longer(
    -c(gene_name, gwas_id, PP.H4.abf_unif),
    names_to = "prior_type",
    values_to = "pp_h4"
  ) |>
  filter(gwas_id == "AUTOIMMUNE", prior_type == "PP.H4.abf_eqtlgen") |>
  rename(pp_h4_unif = PP.H4.abf_unif) |>
  ggplot(aes(pp_h4_unif, pp_h4)) +
  geom_point() +
  geom_vline(xintercept = 0.8, linetype = "dashed", col = "red") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  labs(
    x = TeX("Uniform prior $\\  \\Pr(H_4)$"),
    y = TeX("eQTLGen prior $\\  \\Pr(H_4)$")
  ) +
  theme_jp()

overall_impact_plot <- n_colocs_plot / (abf_eqtlgen_scatter_plot + susie_eqtlgen_scatter_plot) +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  overall_impact_plot_path,
  overall_impact_plot,
  width = 12,
  height = 10
)

# GWAS-eQTL coloc results

wb <- createWorkbook()
addWorksheet(wb, sheetName = "coloc-single")
addWorksheet(wb, sheetName = "coloc-susie")

gwas_eqtl_coloc_abf_data |>
  mutate(
    sig_eqtlgen = PP.H4.abf_eqtlgen > 0.8,
    sig_onek1k = PP.H4.abf_onek1k_r1 > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  mutate(gwas_id = gwas_id_lookup[gwas_id]) |>
  arrange(gwas_id) |>
  filter(sig_eqtlgen != sig_unif | sig_onek1k != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_eqtlgen, PP.H4.abf_onek1k_r1) |>
  writeDataTable(wb, sheet = 1, x = _)

gwas_eqtl_coloc_susie_data |>
  mutate(
    sig_eqtlgen = PP.H4.abf_eqtlgen > 0.8,
    sig_onek1k = PP.H4.abf_onek1k_r1 > 0.8,
    sig_unif = PP.H4.abf_unif > 0.8
  ) |>
  mutate(gwas_id = gwas_id_lookup[gwas_id]) |>
  arrange(gwas_id) |>
  filter(sig_eqtlgen != sig_unif | sig_onek1k != sig_unif) |>
  select(gwas_id, eqtl_id, gene_name,
         PP.H4.abf_unif, PP.H4.abf_eqtlgen, PP.H4.abf_onek1k_r1) |>
  writeDataTable(wb, sheet = 2, x = _)

saveWorkbook(
  wb,
  coloc_results_excel_path,
  overwrite = TRUE
)

bind_rows(
  abf_change_coloc_data |>
    count(prior, gwas_id, type) |>
    mutate(method = "coloc-single"),
  susie_change_coloc_data |>
    count(prior, gwas_id, type) |>
    mutate(method = "coloc-susie")
) |>
  mutate(gwas_id = gwas_id_lookup[gwas_id]) |>
  mutate(unchanged = if_else(substr(type, 1, 9) == "Unchanged", "unchanged", "changed")) |>
  summarise(n = sum(n), .by = c(gwas_id, prior, unchanged)) |>
  mutate(prop = 100 * n / sum(n), .by = c(gwas_id, prior)) |>
  write_csv("output/tables/gwas-eqtl-overall-results.csv")
