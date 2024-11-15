
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(glue)
  library(arrow)
  library(latex2exp)
})

source(here::here("code/plot-utils.R"))

otg_plot_path <- snakemake@output[["otg_plot_path"]]

cs_data <- open_dataset("/home/jp2045/rds/hpc-work/otg-data/v2d_credset") |>
  slice_sample(n = 1e6) |>
  select(postprob, lead_variant_id) |>
  filter(postprob > 0.5) |>
  distinct(postprob, lead_variant_id) |>
  collect()

coloc_data <- open_dataset("/home/jp2045/rds/hpc-work/otg-data/v2d_coloc") |>
  slice_sample(n = 1e6) |>
  select(coloc_h4) |>
  filter(coloc_h4 > 0.5) |>
  collect()

otg_plot <- cs_data |>
  ggplot(aes(postprob)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Posterior causal probability",
    y = "Number of variants"
  ) +
  theme_jp() +
  coloc_data |>
  ggplot(aes(coloc_h4)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = TeX("$\\Pr(H_4)$"),
    y = "Number of loci"
  ) +
  theme_jp() +
  plot_annotation(tag_levels = "a")  &
  theme(plot.tag = element_text(size = 18))

ggsave(
  otg_plot_path,
  otg_plot,
  width = 12,
  height = 10
)
