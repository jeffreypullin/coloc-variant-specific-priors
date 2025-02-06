
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
  library(latex2exp)
  library(ggrepel)
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

# Debugging.
config <- read_yaml("config.yaml")
coloc_abf_paths <- glue(
  "data/output/pqtl-eqtl-coloc-abf-{eqtl_id}-{chr}.rds",
  eqtl_id = rep(config$pqtl_eqtl_coloc_dataset_ids, 22),
  chr = rep(1:22, each = 5)
)
coloc_susie_paths <- glue(
  "data/output/pqtl-eqtl-coloc-susie-{eqtl_id}-{chr}.rds",
  eqtl_id = rep(config$pqtl_eqtl_coloc_dataset_ids, 22),
  chr = rep(1:22, each = 5)
)
protein_metadata <- read_tsv(
  "data/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz",
  show_col_types = FALSE
)

coloc_abf_paths <- snakemake@input[["coloc_abf_paths"]]
coloc_susie_paths <- snakemake@input[["coloc_susie_paths"]]
protein_metadata_path <- snakemake@input[["protein_metadata_path"]]

pqtl_eqtl_perf_both_plot_path <- snakemake@output[["pqtl_eqtl_perf_both_plot_path"]]
prior_effect_plot_path <- snakemake@output[["prior_effect_plot_path"]]

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

to_label <- c(
  "Uniform-Both",
  "eQTLGen-Both",
  "eQTLGen-eQTL",
  "eQTLGen-pQTL",
  "OneK1K (R1)-Both",
  "OneK1K (R1)-eQTL",
  "Gnocchi-eQTL",
  "ABC score-eQTL",
  "ABC score-pQTL"
)

abf_perf_max_plot <- coloc_abf_perf_data_max |>
  mutate(
    method = prior_method_lookup[method],
    dataset = dataset_lookup[dataset]
  ) |>
  mutate(label = paste0(method, "-", dataset)) |>
  mutate(label = if_else(label %in% to_label, label, "")) |>
  ggplot(aes(recall, precision, col = method, shape = dataset, label = label)) +
  geom_point(size = 5, position = position_jitter(h = 2e-04, w = 2e-04)) +
  geom_label_repel(box.padding = 0.4) +
  scale_colour_manual(values = prior_cols) +
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

sig_levels <- function(pp_h4) {
  level <- c(
    seq(0.5, 0.95, by = 0.05)
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
  scale_colour_manual(values = prior_cols) +
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

abf_perf_both_plot <- (abf_perf_max_plot + ggtitle("coloc-single")) / abf_perf_max_curve_plot + 
  plot_layout(guides = "collect")

# coloc.susie()

coloc_susie_perf_data_max <- left_join(
  pqtl_eqtl_coloc_susie_data |>
    select(eqtl_data_id, starts_with("PP.H4.abf"), phenotype_id,
           coloc_gene_id = gene_id, gene_name, idx1_unif, idx2_unif),
  protein_metadata |>
    select(coding_gene_id = gene_id, phenotype_id),
  by = "phenotype_id"
) |>
  rename(
    protein_id = phenotype_id,
    cs_ind1 = idx1_unif,
    cs_ind2 = idx2_unif
  ) |>
  pivot_longer(-c(protein_id, coding_gene_id, eqtl_data_id,
                  coloc_gene_id, gene_name, cs_ind1, cs_ind2),
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

susie_perf_max_plot <- coloc_susie_perf_data_max |>
  mutate(
    method = prior_method_lookup[method],
    dataset = dataset_lookup[dataset]
  ) |>
  mutate(label = paste0(method, "-", dataset)) |>
  mutate(label = if_else(label %in% to_label, label, "")) |>
  ggplot(aes(recall, precision, col = method, shape = dataset, label = label)) +
  geom_point(size = 5, position = position_jitter(h = 2e-04, w = 2e-04)) +
  geom_label_repel(box.padding = 0.4) +
  scale_colour_manual(values = prior_cols) +
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

susie_perf_max_curve_data <- left_join(
  pqtl_eqtl_coloc_susie_data |>
    select(eqtl_data_id, starts_with("PP.H4.abf"), phenotype_id,
           coloc_gene_id = gene_id, gene_name, idx1_unif, idx2_unif),
  protein_metadata |>
    select(coding_gene_id = gene_id, phenotype_id),
  by = "phenotype_id") |>
  rename(
    protein_id = phenotype_id,
    cs_ind1 = idx1_unif,
    cs_ind2 = idx2_unif
  ) |>
  pivot_longer(-c(protein_id, coding_gene_id, eqtl_data_id,
                  coloc_gene_id, gene_name, cs_ind1, cs_ind2),
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

susie_perf_max_curve_plot <- susie_perf_max_curve_data  |>
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
  scale_colour_manual(values = prior_cols) +
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

a <- ((abf_perf_max_plot + ggtitle("coloc-single")) / abf_perf_max_curve_plot)
b <- ((susie_perf_max_plot + ggtitle("coloc-susie")) / susie_perf_max_curve_plot)

pqtl_eqtl_perf_plot <- (a | b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  pqtl_eqtl_perf_both_plot_path,
  plot = pqtl_eqtl_perf_plot,
  width = 14,
  height = 8
)

# Prior impact figure (supplementary).

scatter_plot <- pqtl_eqtl_coloc_abf_data |>
  ggplot(aes(PP.H4.abf_unif, `PP.H4.abf_eqtlgen-eqtl`)) +
  geom_point() +
  geom_vline(xintercept = 0.8, colour = "red") +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(
    x = TeX("Uniform $\\ \\Pr(H_4)$"),
    y = TeX("eQTLGen density $\\ \\Pr(H_4)$")
  ) +
  theme_jp()

diff_scatter_plot <- pqtl_eqtl_coloc_abf_data |>
  mutate(diff = `PP.H4.abf_eqtlgen-eqtl` - PP.H4.abf_unif) |>
  ggplot(aes(PP.H4.abf_unif, diff)) +
  geom_point() +
  geom_abline(slope = -1, linetype = "dotted") +
  geom_abline(slope = -1, intercept = 1, linetype = "dotted") +
  coord_cartesian(ylim = c(-0.75, 0.75)) +
  labs(
    y = TeX("eQTLGen density $\\ \\Pr(H_4) \\ $ - uniform $\\ \\Pr(H_4)$"),
    x = TeX("Uniform $\\ \\Pr(H_4)$")
  ) +
  theme_jp()

coloc_abf_diff_data <- pqtl_eqtl_coloc_abf_data |>
  select(gene_name, starts_with("PP.H4.abf"), max_lbf) |>
  pivot_longer(
    -c(gene_name, PP.H4.abf_unif, max_lbf),
    names_to = "prior_type",
    values_to = "pp_h4"
  ) |>
  separate_wider_delim(
    prior_type,
    delim = "-",
    names = c("method", "dataset"),
    too_few = "align_start"
  ) |>
  separate_wider_delim(
    method,
    delim = "_",
    names = c("junk", "method"),
    too_many = "merge"
  ) |>
  select(-junk) |>
  mutate(dataset = if_else(is.na(dataset), "pqtl_eqtl", dataset)) |>
  mutate(
    method = prior_method_lookup[method],
    dataset = dataset_lookup[dataset]
  ) |>
  mutate(prior_type = paste0(method, "-", dataset)) |>
  mutate(diff = abs(pp_h4 - PP.H4.abf_unif)) |>
  mutate(method = "coloc-single")

coloc_susie_diff_data <- pqtl_eqtl_coloc_susie_data |>
  select(gene_name, starts_with("PP.H4.abf"), max_lbf) |>
  pivot_longer(
    -c(gene_name, PP.H4.abf_unif, max_lbf),
    names_to = "prior_type",
    values_to = "pp_h4"
  ) |>
  separate_wider_delim(
    prior_type,
    delim = "-",
    names = c("method", "dataset"),
    too_few = "align_start"
  ) |>
  separate_wider_delim(
    method,
    delim = "_",
    names = c("junk", "method"),
    too_many = "merge"
  ) |>
  select(-junk) |>
  mutate(dataset = if_else(is.na(dataset), "pqtl_eqtl", dataset)) |>
  mutate(
    method = prior_method_lookup[method],
    dataset = dataset_lookup[dataset]
  ) |>
  mutate(prior_type = paste0(method, "-", dataset)) |>
  mutate(diff = abs(pp_h4 - PP.H4.abf_unif)) |>
  mutate(method = "coloc-susie")

prop_small_diff_plot <- bind_rows(
  coloc_abf_diff_data,
  coloc_susie_diff_data
) |>
  summarise(
    prop_small_diff = sum(diff < 0.01) / n(),
    .by = c(prior_type, method)
  ) |>
  mutate(prior_type = fct_reorder(factor(prior_type), prop_small_diff)) |>
  ggplot(aes(prior_type, prop_small_diff, fill = method)) +
  geom_col(position = "dodge") +
  labs(
    y = TeX("Prop. abs($\\ \\Pr(H_4)\\ $ difference) < 0.01"),
    x = "Prior type"
  ) +
  coord_flip(ylim = c(0.75, 1)) +
  theme_jp_vgrid()

max_lbf_by_method_plot <- bind_rows(
  coloc_abf_diff_data,
  coloc_susie_diff_data
  ) |>
  summarise(diff = median(diff), .by = c(PP.H4.abf_unif, max_lbf, method)) |>
  ggplot(aes(method, max_lbf)) +
  geom_boxplot() +
  labs(
    x = "Coloc method",
    y = "Maximum log Bayes factor"
  ) + 
  theme_jp()

max_lbf_diff_plot <- coloc_abf_diff_data |>
  summarise(diff = median(diff), .by = c(PP.H4.abf_unif, max_lbf, method)) |>
  ggplot(aes(max_lbf, diff)) +
  geom_point() +
  labs(
    x = "Maximum log Bayes factor",
    y = "Unifrom vs variant-specific difference"
  ) +
  theme_jp()

prior_effect_plot <- (prop_small_diff_plot |
  (max_lbf_by_method_plot / max_lbf_diff_plot) |
  (scatter_plot / diff_scatter_plot)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  prior_effect_plot_path,
  plot = prior_effect_plot,
  width = 14,
  height = 10
)
