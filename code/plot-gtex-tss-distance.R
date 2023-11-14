library(readr)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(forcats)

gtex_paths <- snakemake@input[["gtex_paths"]]
boxplot_by_tissue_path <- snakemake@output[["boxplot_by_tissue_path"]]

process_file_gtex <- function(path) {

  raw_data <- read_tsv(path,
                       col_select = c(gene_id, pval_nominal, tss_distance),
                       show_col_types = FALSE, progress = FALSE)
  raw_data |>
    filter(pval_nominal < 5 * 10^-8) |>
    group_by(gene_id) |>
    arrange(pval_nominal) |>
    filter(row_number() == 1) |>
    mutate(abs_tss_distance = abs(tss_distance)) |>
    filter(abs_tss_distance > 0) |>
    ungroup() |>
    pull(abs_tss_distance)
}

gtex_data <- tibble(path = gtex_paths) |>
  rowwise() |>
  mutate(tissue = gsub("_", " ", unlist(strsplit(basename(path), "[.]"))[[1]])) |>
  mutate(abs_tss_distance = list(process_file_gtex(path))) |>
  ungroup() |>
  unnest(cols = abs_tss_distance)

tss_boxplot_by_tissue <- gtex_data |>
  mutate(tissue = fct_rev(factor(tissue))) |>
  ggplot(aes(tissue, abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  labs(
    x = "Tissue",
    y = "Distance to TSS"
  ) +
  coord_flip() +
  theme_bw()
ggsave(boxplot_by_tissue_path, plot = tss_boxplot_by_tissue)

# gtex_data |>
#   ggplot(aes(abs_tss_distance)) +
#   geom_histogram(binwidth = 1000, colour = "black", fill = "grey") +
#   coord_cartesian(xlim = c(0, 500000)) +
#   labs(
#     y = "Count",
#     x = "Absolute distance from variant to TSS",
#     title = "eQTL distance distribution in GTEx v8"
#   ) +
#   theme_bw()
