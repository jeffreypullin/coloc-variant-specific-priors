
source(here::here("renv/activate.R"))

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(yaml)
library(forcats)
library(glue)
library(seqminer)

source(here::here("code/coloc-utils.R"))
source(here::here("code/plot-utils.R"))

config <- read_yaml(here::here("config.yaml"))
eqtl_data_ids <- config$pqtl_eqtl_coloc_dataset_ids
eqtl_data_ids_lookup <- c(
  "QTD000373" = "Lepik et. al. 2017",
  "QTD000341" = "GTEx thyroid",
  "QTD000116" = "GTEx adipose (subcutaneous)",
  "QTD000021" = "BLUEPRINT Monocytes",
  "QTD000539" = "TwinsUK LCL"
)

pqtl_eqtl_coloc_abf_data <- expand_grid(
  eqtl_data_id = eqtl_data_ids,
  chr = 1:22
) |>
  rowwise() |>
  mutate(file = paste0(
    here::here("data/output/pqtl-eqtl-coloc-abf"),
    "-", eqtl_data_id, "-", chr, ".rds")
  ) |>
  mutate(data = list(read_rds(file))) |>
  unnest(data)

pqtl_eqtl_coloc_susie_data <- expand_grid(
  eqtl_data_id = eqtl_data_ids,
  chr = 1:22
) |>
  rowwise() |>
  mutate(file = paste0(
    here::here("data/output/pqtl-eqtl-coloc-susie"),
    "-", eqtl_data_id, "-", chr, ".rds")
  ) |>
  mutate(data = list(read_rds(file))) |>
  unnest(data)

protein_metadata <- read_tsv(
  here::here("data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"),
  show_col_types = FALSE
)

# coloc.abf()

coloc_abf_perf_data <- left_join(
  pqtl_eqtl_coloc_abf_data |>
    select(eqtl_data_id, starts_with("PP.H4.abf"), phenotype_id,
           coloc_gene_id = gene_id, gene_name),
  protein_metadata |>
    select(coding_gene_id = gene_id, phenotype_id),
  by = "phenotype_id"
) |>
  rename(protein_id = phenotype_id) |>
  pivot_longer(-c(protein_id, coding_gene_id, eqtl_data_id,
                  coloc_gene_id, gene_name),
               values_to = "pp_h4",
               names_to = c("junk", "prior_type"),
               names_pattern = "(.*?)_(.*)") |>
  select(-junk) |>
  mutate(
    sig_coloc = pp_h4 >= 0.8,
    gene_match = coding_gene_id == coloc_gene_id
  ) |>
  select(-c(coding_gene_id, coloc_gene_id, pp_h4)) |>
  mutate(
    tp = sig_coloc & gene_match,
    fp = sig_coloc & !gene_match,
    fn = !sig_coloc & gene_match,
  ) |>
  summarise(
    n_tp = sum(tp),
    n_fp = sum(fp),
    n_fn = sum(fn),
    .by = c(eqtl_data_id, prior_type)
  ) |>
  mutate(
    recall = n_tp / (n_tp + n_fn),
    precision = n_tp / (n_tp + n_fp)
  )

pqtl_eqtl_abf_perf_by_eqtl_data_plot <- coloc_abf_perf_data |>
  mutate(is_unif = prior_type == "unif") |>
  mutate(eqtl_data_name = eqtl_data_ids_lookup[eqtl_data_id]) |>
  ggplot(aes(recall, precision, col = is_unif, label = prior_type)) +
  geom_point() +
  scale_colour_manual(values = c("#66CCEE", "#EE6677")) +
  geom_text_repel() +
  labs(
    x = "Recall",
    y = "Precision",
  ) +
  facet_wrap(~eqtl_data_name, scales = "free") +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  "output/figures/pqtl-eqtl-abf-perf-by-eqtl-data-plot.pdf",
  plot = pqtl_eqtl_abf_perf_by_eqtl_data_plot,
  width = 12,
  height = 10
)

pqtl_eqtl_abf_perf_median_plot <- coloc_abf_perf_data |>
  summarise(
    recall = median(recall),
    precision = median(precision),
    .by = prior_type
  ) |>
  mutate(is_unif = prior_type == "unif") |>
  ggplot(aes(recall, precision, col = is_unif, label = prior_type)) +
  geom_point() +
  geom_text_repel() +
  scale_colour_manual(values = c("#66CCEE", "#EE6677")) +
  labs(
    x = "Recall",
    y = "Precision",
    title = "coloc.abf(), median over eQTL datasets"
  ) +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  "output/figures/pqtl-eqtl-abf-perf-median-data-plot.pdf",
  plot = pqtl_eqtl_abf_perf_median_plot,
  width = 8,
  height = 6
)

pqtl_eqtl_abf_n_colocs_plot <- coloc_abf_perf_data |>
  mutate(eqtl_data_name = eqtl_data_ids_lookup[eqtl_data_id]) |>
  mutate(n_colocalised = n_tp + n_fp) |>
  mutate(prior_type = fct_reorder(
    factor(prior_type),
    n_colocalised
  )) |>
  mutate(eqtl_data_name = fct_reorder(
    factor(eqtl_data_name),
    n_colocalised
  )) |>
  ggplot(aes(prior_type, n_colocalised, fill = eqtl_data_name)) + 
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    y = "Number of signigicantly colocalised pQTL-eQTL pairs",
    x = "Prior"
  ) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

ggsave(
  "output/figures/pqtl-eqtl-abf-n-colocs-plot.pdf",
  plot = pqtl_eqtl_abf_n_colocs_plot,
  width = 8,
  height = 6
)

pph4_scatter_plot <- pqtl_eqtl_coloc_abf_data |>
  ggplot(aes(PP.H4.abf_unif, PP.H4.abf_eqtl_tss_eqtlgen)) + 
  geom_point() +
  geom_vline(xintercept = 0.8, colour = "red") +
  geom_hline(yintercept = 0.8, colour = "red") +
  theme_jp()

ggsave(
  "output/figures/pqtl-eqtl-abf-pph4-scatter-plot.pdf",
  plot = pph4_scatter_plot,
  width = 8,
  height = 6
)

pqtl_eqtl_diff_trend_plot <- pqtl_eqtl_coloc_abf_data |>
  mutate(diff = PP.H4.abf_unif - PP.H4.abf_eqtl_tss_eqtlgen) |>
  ggplot(aes(PP.H4.abf_unif, diff)) +
  geom_point() +
  coord_cartesian(ylim = c(-0.75, 0.75)) +
  theme_jp()

ggsave(
  "output/figures/pqtl-eqtl-abf-diff-trend-plot.pdf",
  plot = pqtl_eqtl_diff_trend_plot,
  width = 8,
  height = 6
)

pqtl_eqtl_diff_mag_plot <- pqtl_eqtl_coloc_abf_data |>
  mutate(diff = abs(PP.H4.abf_unif - PP.H4.abf_eqtl_tss_eqtlgen)) |>
  ggplot(aes(diff)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Magnitudie of difference",
    y = "Count"
  ) +
  theme_jp()

ggsave(
  "output/figures/pqtl-eqtl-abf-diff-mag-plot.pdf",
  plot = pqtl_eqtl_diff_mag_plot,
  width = 8,
  height = 6
)

# coloc.susie()

coloc_susie_perf_data <- left_join(
  pqtl_eqtl_coloc_susie_data |>
    select(eqtl_data_id, starts_with("PP.H4.abf"), phenotype_id, 
           coloc_gene_id = gene_id, gene_name), 
  protein_metadata |>
    select(coding_gene_id = gene_id, phenotype_id),
  by = "phenotype_id"
) |>
  rename(protein_id = phenotype_id) |>
  pivot_longer(-c(protein_id, coding_gene_id, eqtl_data_id, 
                  coloc_gene_id, gene_name),
               values_to = "pp_h4", 
               names_to = c("junk", "prior_type"),
               names_pattern = "(.*?)_(.*)") |>
  select(-junk) |>
  mutate(sig_coloc = pp_h4 >= 0.8) |>
  mutate(gene_match = coding_gene_id == coloc_gene_id) |>
  select(-c(coding_gene_id, coloc_gene_id, pp_h4)) |>
  mutate(
    tp = sig_coloc & gene_match,
    fp = sig_coloc & !gene_match,
    fn = !sig_coloc & gene_match,
  ) |>
  summarise(
    n_tp = sum(tp),
    n_fp = sum(fp),
    n_fn = sum(fn),
    .by = c(eqtl_data_id, prior_type)
  ) |>
  mutate(
    recall = n_tp / (n_tp + n_fn),
    precision = n_tp / (n_tp + n_fp)
  )

pqtl_eqtl_susie_perf_by_eqtl_data_plot <- coloc_susie_perf_data |>
  mutate(is_unif = prior_type == "unif") |>
  mutate(eqtl_data_name = eqtl_data_ids_lookup[eqtl_data_id]) |>
  ggplot(aes(recall, precision, col = is_unif, label = prior_type)) +
  geom_point() +
  scale_colour_manual(values = c("#66CCEE", "#EE6677")) +
  geom_text_repel() +
  labs(
    x = "Recall",
    y = "Precision",
  ) +
  facet_wrap(~eqtl_data_name, scales = "free") +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  "output/figures/pqtl-eqtl-susue-perf-by-eqtl-data-plot.pdf",
  plot = pqtl_eqtl_susie_perf_by_eqtl_data_plot,
  width = 12,
  height = 10
)

pqtl_eqtl_susie_perf_median_plot <- coloc_susie_perf_data |>
  summarise(
    recall = median(recall),
    precision = median(precision),
    .by = prior_type
  ) |>
  mutate(is_unif = prior_type == "unif") |>
  ggplot(aes(recall, precision, col = is_unif, label = prior_type)) +
  geom_point() +
  geom_text_repel() +
  scale_colour_manual(values = c("#66CCEE", "#EE6677")) +
  labs(
    x = "Recall",
    y = "Precision",
    title = "coloc.susie(), median over eQTL datasets"
  ) +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  "output/figures/pqtl-eqtl-susie-perf-median-data-plot.pdf",
  plot = pqtl_eqtl_susie_perf_median_plot,
  width = 8,
  height = 6
)

pph4_scatter_plot <- pqtl_eqtl_coloc_susie_data |>
  ggplot(aes(PP.H4.abf_unif, PP.H4.abf_eqtl_tss_pqtl_tss_eqtlgen)) + 
  geom_point() +
  geom_vline(xintercept = 0.8, colour = "red") +
  geom_hline(yintercept = 0.8, colour = "red") +
  theme_jp()

ggsave(
  "output/figures/pqtl-eqtl-susie-pph4-scatter-plot.pdf",
  plot = pph4_scatter_plot,
  width = 8,
  height = 6
)

