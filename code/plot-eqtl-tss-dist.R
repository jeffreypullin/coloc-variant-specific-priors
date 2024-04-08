
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(tools)
  library(forcats)
  library(janitor)
  library(ggupset)
  library(purrr)
  library(patchwork)
  library(dplyr, warn.conflicts = FALSE)
})

source(here::here("code", "plot-utils.R"))

eqtl_catalogue_data <- read_rds(snakemake@input[["eqtl_catalogue_data_file"]])
gtex_data <- read_rds(snakemake@input[["gtex_data_file"]])
eqtlgen_data <- read_rds(snakemake@input[["eqtlgen_data_file"]])
onek1k_data <- read_rds(snakemake@input[["onek1k_data_file"]])

dist_plot_file <- snakemake@output[["dist_plot_file"]]
onek1k_plot_file <- snakemake@output[["onek1k_plot_file"]]
dataset_plot_file <- snakemake@output[["dataset_plot_file"]]

# Distance plot.

eqtl_catalogue_dist_plot <- eqtl_catalogue_data |>
  filter(abs_tss_distance > 0) |>
  mutate(abs_tss_distance = log10(abs_tss_distance)) |>
  ggplot(aes(file, abs_tss_distance, fill = study_label)) +
  geom_boxplot(outlier.alpha = 0.1) +
  coord_flip() +
  labs(
    x = "eQTL catalogue datasets",
    y = "TSS distance (log10)",
    fill = "Study",
    title = "eQTL catalogue"
  ) +
  theme_jp_vgrid() +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
  )

gtex_dist_plot <- gtex_data |>
  mutate(tissue = fct_rev(factor(tissue))) |>
  filter(abs_tss_distance > 0) |>
  mutate(abs_tss_distance = log10(abs_tss_distance)) |>
  ggplot(aes(tissue, abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  labs(
    x = "Tissue",
    y = "TSS distance (log10)",
    title = "GTEx v8"
  ) +
  coord_flip() +
  theme_jp_vgrid() +
  theme(
    axis.text.y = element_text(size = 8)
  )

dist_plot <- gtex_dist_plot + eqtl_catalogue_dist_plot +
  plot_layout(axis_titles = "collect")

ggsave(
  dist_plot_file,
  dist_plot,
  width = 10,
  height = 6
)

# OneK1K analysis.

onek1k_data <- onek1k_data |>
  group_by(gene_id, tss_distance, round) |>
  summarise(
    cell_type_collapsed = paste0(sort(unlist(cell_type)), collapse = "-"),
    cell_type = list(cell_type),
    n_cell_types = n(),
    .groups = "drop"
  ) |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_tss_distance = log10(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

print(onek1k_data)

n_cell_types_plot <- onek1k_data |>
  mutate(abs_tss_distance = log10(abs_tss_distance)) |>
  ggplot(aes(factor(n_cell_types), abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  labs(
    x = "Number of cell types with eQTL",
    y = "TSS distance (log10)",
  ) +
  theme_jp()

cond_round_plot <- onek1k_data |>
  mutate(abs_tss_distance = log10(abs_tss_distance)) |>
  ggplot(aes(factor(round), abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  labs(
    x = "Round of conditional analysis",
    y = "TSS distance (log10)",
  ) +
  theme_jp()

cell_type_plot <- onek1k_data |>
  filter(n_cell_types == 1) |>
  mutate(cell_type = unlist(cell_type)) |>
  mutate(abs_tss_distance = log10(abs_tss_distance)) |>
  mutate(cell_type = fct_reorder(factor(cell_type), abs_tss_distance)) |>
  ggplot(aes(cell_type, abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  coord_flip() +
  labs(
    y = "Distance",
    x = "Cell type"
  ) +
  theme_jp_vgrid() +
  theme(axis.text = element_text(size = 10))

onek1k_plot <- (cond_round_plot + cell_type_plot) / n_cell_types_plot +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(
  onek1k_plot_file,
  onek1k_plot,
  width = 10,
  height = 8
)

# Dataset plot.

all_studies_data <- bind_rows(
  gtex_data |>
    mutate(study = "GTEx v8"),
  eqtlgen_data |>
    mutate(study = "eQTLGen"),
  eqtl_catalogue_data |>
    mutate(study = "eQTL Catalogue"),
  onek1k_data |>
    mutate(study = "OneK1K")
)

all_studies_dist_plot <- all_studies_data |>
  mutate(log10_abs_tss_distance = log10(abs_tss_distance)) |>
  ggplot(aes(study, log10_abs_tss_distance)) +
  geom_violin(fill = "grey") +
  geom_boxplot(outlier.alpha = 0.1) +
  labs(
    x = "Dataset",
    y = "TSS distance (log10)",
  ) +
  coord_flip() +
  theme_jp_vgrid()

all_studies_density_plot <- all_studies_data |>
  filter(abs(tss_distance) <= 2e5) |>
  mutate(kb = tss_distance / 1000) |>
  ggplot(aes(kb)) +
  geom_density(fill = "grey", bw = "SJ") +
  facet_wrap(~ study, scales = "free_y") +
  theme_jp() +
  labs(
    y = "",
    x = "Distance to TSS (kb)"
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 20, hjust = 0)
  )

dataset_plot <- all_studies_dist_plot + all_studies_density_plot

ggsave(
  dataset_plot_file,
  dataset_plot,
  width = 10,
  height = 8
)
