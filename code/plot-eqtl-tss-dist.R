
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
eqtlgen_data <- read_rds(snakemake@input[["eqtlgen_data_file"]])
onek1k_data <- read_rds(snakemake@input[["onek1k_data_file"]])

dist_plot_file <- snakemake@output[["dist_plot_file"]]
onek1k_plot_file <- snakemake@output[["onek1k_plot_file"]]

# Debugging.
eqtl_catalogue_data <- read_rds("data/processed-data/eqtl-catalogue.rds")
onek1k_data <- read_rds("data/processed-data/onek1k.rds")
eqtlgen_data <- read_rds("data/processed-data/eqtlgen.rds")

# Distance plot.

all_studies_data <- bind_rows(
  eqtlgen_data |>
    mutate(
      study = "eQTLGen",
      study_label = "eQTLGen",
      file = "eQTLGen"
    ),
  eqtl_catalogue_data |>
    mutate(study = "eQTL Catalogue"),
  onek1k_data |>
    mutate(
      study = "OneK1K",
      study_label = "OneK1K",
      file = "OneK1K"
    )
) |>
  filter(tss_distance != 0) |>
  mutate(
    abs_tss_distance = abs(tss_distance),
    log10_abs_tss_distance = log10(abs_tss_distance)
  )

dist_plot <- all_studies_data |>
  # Remove pQTL dataset.
  filter(study_label != "Sun_2018") |>
  filter(abs_tss_distance > 0) |>
  mutate(study = if_else(study_label == "GTEx", "GTEx", study)) |>
  mutate(file = fct_reorder(factor(file), abs_tss_distance)) |>
  ggplot(aes(file, abs_tss_distance, fill = study)) +
  geom_boxplot(outlier.alpha = 0.05, lwd = 0.3, colour = "grey50") +
  coord_flip() +
  scale_y_log10() +
  labs(
    x = "Dataset",
    y = "Distance to TSS (bp)",
    fill = "Study",
  ) +
  theme_jp_vgrid() +
  theme(
    legend.title = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
  )

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
  )

n_cell_types_plot <- onek1k_data |>
  ggplot(aes(factor(n_cell_types), abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  scale_y_log10() +
  labs(
    x = "Number of cell types with eQTL",
    y = "Distance to TSS (bp)",
  ) +
  theme_jp()

cond_round_plot <- onek1k_data |>
  ggplot(aes(factor(round), abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  scale_y_log10() +
  labs(
    x = "Round of conditional analysis",
    y = "Distance to TSS (bp)",
  ) +
  theme_jp()

cell_type_plot <- onek1k_data |>
  filter(n_cell_types == 1) |>
  mutate(cell_type = unlist(cell_type)) |>
  mutate(cell_type = fct_reorder(factor(cell_type), abs_tss_distance)) |>
  ggplot(aes(cell_type, abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  scale_y_log10() +
  coord_flip() +
  labs(
    y = "Distance to TSS (bp)",
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
