
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(forcats)
  library(glue)
  library(seqminer)
  library(stringr)
  library(yaml)
})

source(here::here("code/coloc-utils.R"))
source(here::here("code/plot-utils.R"))

eqtl_data_ids_lookup <- c(
  "QTD000373" = "Lepik et. al. 2017",
  "QTD000341" = "GTEx thyroid",
  "QTD000116" = "GTEx adipose (subcutaneous)",
  "QTD000021" = "BLUEPRINT Monocytes",
  "QTD000539" = "TwinsUK LCL"
)

coloc_abf_paths <- snakemake@input[["coloc_abf_paths"]]
coloc_susie_paths <- snakemake@input[["coloc_susie_paths"]]
protein_metadata_path <- snakemake@input[["protein_metadata_path"]]

abf_perf_by_dataset_plot_path <- snakemake@output[["abf_perf_by_dataset_plot_path"]]
abf_perf_median_plot_path <- snakemake@output[["abf_perf_median_plot_path"]]
abf_perf_max_plot_path <- snakemake@output[["abf_perf_max_plot_path"]]
abf_n_colocs_plot_path <- snakemake@output[["abf_n_colocs_plot_path"]]
abf_pph4_scatter_plot_path <- snakemake@output[["abf_pph4_scatter_plot_path"]]
susie_perf_by_dataset_plot_path <- snakemake@output[["susie_perf_by_dataset_plot_path"]]
susie_perf_median_plot_path <- snakemake@output[["susie_perf_median_plot_path"]]
susie_perf_max_plot_path <- snakemake@output[["susie_perf_max_plot_path"]]
susie_pph4_scatter_plot_path <- snakemake@output[["susie_pph4_scatter_plot_path"]]
abf_perf_max_curve_plot_path <- snakemake@output[["abf_perf_max_curve_plot_path"]]
pqtl_eqtl_perf_plot_path <- snakemake@output[["pqtl_eqtl_perf_plot_path"]]

# Debugging.
#config <- read_yaml("config.yaml")
#coloc_abf_paths <- glue(
#  "data/output/pqtl-eqtl-coloc-abf-{eqtl_id}-{chr}.rds",
#  eqtl_id = rep(config$pqtl_eqtl_coloc_dataset_ids, 22),
#  chr = rep(1:22, each = 5)
#)
#coloc_susie_paths <- glue(
#  "data/output/pqtl-eqtl-coloc-susie-{eqtl_id}-{chr}.rds",
#  eqtl_id = rep(config$pqtl_eqtl_coloc_dataset_ids, 22),
#  chr = rep(1:22, each = 5)
#)
#protein_metadata <- read_tsv(
#  "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz", 
#  show_col_types = FALSE
#)

pqtl_eqtl_coloc_abf_data <- tibble(path = coloc_abf_paths) |>
  mutate(eqtl_data_id = str_extract(path, "QTD[0-9]{6}")) |>
  rowwise() |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

pqtl_eqtl_coloc_susie_data <- tibble(path = coloc_susie_paths) |>
  mutate(eqtl_data_id = str_extract(path, "QTD[0-9]{6}")) |>
  rowwise() |>
  mutate(data = list(read_rds(path))) |>
  unnest(data)

protein_metadata <- read_tsv(protein_metadata_path, show_col_types = FALSE)

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

abf_perf_by_dataset_plot <- coloc_abf_perf_data |>
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
  abf_perf_by_dataset_plot_path,
  plot = abf_perf_by_dataset_plot,
  width = 12,
  height = 10
)

abf_perf_median_plot <- coloc_abf_perf_data |>
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
  abf_perf_median_plot_path,
  plot = abf_perf_median_plot,
  width = 8,
  height = 6
)

abf_n_colocs_plot <- coloc_abf_perf_data |>
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
  abf_n_colocs_plot_path,
  plot = abf_n_colocs_plot,
  width = 8,
  height = 6
)

coloc_abf_perf_data_max <- left_join(
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
  summarise(
    tp = any(sig_coloc) & all(gene_match),
    fp = any(sig_coloc) & !all(gene_match),
    fn = !any(sig_coloc) & all(gene_match),
    .by = c(prior_type, gene_name, protein_id)
  ) |>
  summarise(
    n_tp = sum(tp),
    n_fp = sum(fp),
    n_fn = sum(fn),
    .by = c(prior_type)
  ) |>
  mutate(
    recall = n_tp / (n_tp + n_fn),
    precision = n_tp / (n_tp + n_fp)
  ) |>
  separate_wider_delim(
    prior_type,
    delim = "-",
    names = c("method", "dataset"),
    too_few = "align_start"
  ) |>
  mutate(dataset = if_else(is.na(dataset), "pqtl_eqtl", dataset))

dataset_lookup <- c(
  "eqtl" = "eQTL",
  "pqtl" = "pQTL",
  "pqtl_eqtl" = "Both"
)

abf_perf_max_plot <- coloc_abf_perf_data_max |>
  mutate(
    method = prior_method_lookup[method],
    dataset = dataset_lookup[dataset]
  ) |>
  ggplot(aes(recall, precision, col = method, shape = dataset)) +
  geom_point(size = 5, alpha = 0.7, position = position_jitter(h = 3e-04, w = 3e-04)) +
  labs(
    x = "Recall",
    y = "Precision",
  ) +
  theme_jp() +
  labs(
    col = "Method",
    shape = "Dataset"
  ) +
  theme(
    legend.title = element_text(family = "Helvetica", size = 14, color = "#222222"),
    legend.position = "none"
  )

ggsave(
  abf_perf_max_plot_path,
  plot = abf_perf_max_plot,
  width = 8,
  height = 6
)

sig_levels <- function(pp_h4) {
  level <- c(
    seq(0.01, 0.09, by = 0.01),
    seq(0.1, 0.9, by = 0.1),
    seq(0.9, 0.99, by = 0.01)
  )
  tibble(
    level = level,
    sig_coloc = pp_h4 > level
  )
}

abf_perf_max_curve_data <- left_join(
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
  rowwise() |>
  mutate(
    sig_coloc = list(sig_levels(pp_h4)),
    gene_match = coding_gene_id == coloc_gene_id
  ) |>
  ungroup() |>
  select(-c(coding_gene_id, coloc_gene_id, pp_h4)) |>
  unnest(sig_coloc) |>
  summarise(
    tp = any(sig_coloc) & all(gene_match),
    fp = any(sig_coloc) & !all(gene_match),
    fn = !any(sig_coloc) & all(gene_match),
    tn = !any(sig_coloc) & !all(gene_match),
    .by = c(prior_type, gene_name, protein_id, level)
  ) |>
  summarise(
    n_tp = sum(tp),
    n_fp = sum(fp),
    n_fn = sum(fn),
    n_tn = sum(tn),
    .by = c(prior_type, level)
  ) |>
  mutate(
    tpr = n_tp / (n_tp + n_fn),
    fpr = n_fp / (n_fp + n_tn)
  )

abf_perf_max_curve_plot <- abf_perf_max_curve_data |>
  separate_wider_delim(
    prior_type,
    delim = "-",
    names = c("method", "dataset"),
    too_few = "align_start"
  ) |>
  mutate(dataset = if_else(is.na(dataset), "pqtl_eqtl", dataset)) |>
  mutate(
    method = prior_method_lookup[method],
    dataset = dataset_lookup[dataset]
  ) |>
  ggplot(aes(fpr, tpr, color = method, shape = dataset)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    col = "Method",
    shape = "Dataset"
  ) +
  theme_jp() +
  theme(
    legend.title = element_text(family = "Helvetica", size = 14, color = "#222222"),
    legend.position = "right"
  )

ggsave(
  abf_perf_max_curve_plot_path,
  plot = abf_perf_max_curve_plot,
  width = 12,
  height = 10
)

pqtl_eqtl_perf_plot <- abf_perf_max_plot + abf_perf_max_curve_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  pqtl_eqtl_perf_plot_path,
  plot = pqtl_eqtl_perf_plot,
  width = 14,
  height = 8
)

pph4_scatter_plot <- pqtl_eqtl_coloc_abf_data |>
  ggplot(aes(PP.H4.abf_unif, `PP.H4.abf_eqtlgen-eqtl`)) +
  geom_point() +
  geom_vline(xintercept = 0.8, colour = "red") +
  geom_hline(yintercept = 0.8, colour = "red") +
  theme_jp()

ggsave(
  abf_pph4_scatter_plot_path,
  plot = pph4_scatter_plot,
  width = 8,
  height = 6
)

# FIXME: Should these be made with the GWAS-eQTL colocs instead?
#pqtl_eqtl_diff_trend_plot <- pqtl_eqtl_coloc_abf_data |>
#  mutate(diff = PP.H4.abf_unif - PP.H4.abf_eqtl_tss_eqtlgen) |>
#  ggplot(aes(PP.H4.abf_unif, diff)) +
#  geom_point() +
#  coord_cartesian(ylim = c(-0.75, 0.75)) +
#  theme_jp()

#ggsave(
#  "output/figures/pqtl-eqtl-abf-diff-trend-plot.pdf",
#  plot = pqtl_eqtl_diff_trend_plot,
#  width = 8,
#  height = 6
#)

#pqtl_eqtl_diff_mag_plot <- pqtl_eqtl_coloc_abf_data |>
#  mutate(diff = abs(PP.H4.abf_unif - PP.H4.abf_eqtl_tss_eqtlgen)) |>
#  ggplot(aes(diff)) +
#  geom_histogram(binwidth = 0.05) +
#  labs(
#    x = "Magnitudie of difference",
#    y = "Count"
#  ) +
#  theme_jp()

#ggsave(
#  "output/figures/pqtl-eqtl-abf-diff-mag-plot.pdf",
#  plot = pqtl_eqtl_diff_mag_plot,
#  width = 8,
#  height = 6
#)

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

susie_perf_by_dataset_plot <- coloc_susie_perf_data |>
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
  susie_perf_by_dataset_plot_path,
  plot = susie_perf_by_dataset_plot,
  width = 12,
  height = 10
)

susie_perf_median_plot <- coloc_susie_perf_data |>
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
  susie_perf_median_plot_path,
  plot = susie_perf_median_plot,
  width = 8,
  height = 6
)

pqtl_eqtl_coloc_susie_data |>
  select(eqtl_data_id, phenotype_id, gene_name, hit1_unif, hit2_unif)

coloc_susie_perf_data_max <- left_join(
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
  mutate(
    sig_coloc = pp_h4 >= 0.8,
    gene_match = coding_gene_id == coloc_gene_id
  ) |>
  select(-c(coding_gene_id, coloc_gene_id, pp_h4)) |>
  summarise(
    tp = any(sig_coloc) & all(gene_match),
    fp = any(sig_coloc) & !all(gene_match),
    fn = !any(sig_coloc) & all(gene_match),
    .by = c(prior_type, gene_name, protein_id)
  ) |>
  summarise(
    n_tp = sum(tp),
    n_fp = sum(fp),
    n_fn = sum(fn),
    .by = c(prior_type)
  ) |>
  mutate(
    recall = n_tp / (n_tp + n_fn),
    precision = n_tp / (n_tp + n_fp)
  )

susie_perf_max_plot <- coloc_susie_perf_data_max |>
  mutate(is_unif = prior_type == "unif") |>
  ggplot(aes(recall, precision, col = is_unif, label = prior_type)) +
  geom_point() +
  geom_text_repel() +
  scale_colour_manual(values = c("#66CCEE", "#EE6677")) +
  labs(
    x = "Recall",
    y = "Precision",
    title = "coloc.susie(), max over eQTL datasets"
  ) +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  susie_perf_max_plot_path,
  plot = susie_perf_max_plot,
  width = 8,
  height = 6
)

susie_pph4_scatter_plot <- pqtl_eqtl_coloc_susie_data |>
  ggplot(aes(PP.H4.abf_unif, `PP.H4.abf_eqtlgen-pqtl_eqtl`)) +
  geom_point() +
  geom_vline(xintercept = 0.8, colour = "red") +
  geom_hline(yintercept = 0.8, colour = "red") +
  theme_jp()

ggsave(
  susie_pph4_scatter_plot_path,
  plot = susie_pph4_scatter_plot,
  width = 8,
  height = 6
)
