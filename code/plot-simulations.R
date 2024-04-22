
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(latex2exp)
  library(stringr)
  library(tools)
  library(glue)
  library(yaml)
  library(tidyr)
  library(patchwork)
})

source("code/plot-utils.R")

devtools::load_all("~/coloc")

# Debugging.
# config <- yaml::read_yaml("config.yaml")
# simulation_paths <- glue(
# "data/output/sim-result-{gene}.rds",
#  gene = config$simulation_genes
# )

simulation_paths <- snakemake@input

simulation_data <- lapply(simulation_paths, function(x) {
  gene_name <- file_path_sans_ext(x) |>
    basename() |>
    str_extract("([^-]+$)")

  data <- read_rds(x) |>
    mutate(gene = gene_name)
}) |>
  bind_rows()

p1 <- simulation_data |>
  pivot_wider(
    names_from = "prior",
    values_from = "pp_h4",
    values_fn = list
  ) |>
  unnest(cols = c(unif, non_unif)) |>
  mutate(diff = non_unif - unif) |>
  ggplot(aes(dist, diff)) +
  geom_point() +
  facet_wrap(~hyp) +
  theme_jp()

p2 <- simulation_data |>
  pivot_wider(
    names_from = "prior",
    values_from = "pp_h4",
    values_fn = list
  ) |>
  unnest(cols = c(unif, non_unif)) |>
  mutate(diff = non_unif - unif) |>
  ggplot(aes(unif, diff)) +
  geom_point() +
  facet_wrap(~hyp) +
  theme_jp()

gene_plot <- simulation_data |>
  ggplot(aes(hyp, pp_h4, fill = prior)) +
  geom_boxplot(position = "dodge") +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "True hypothesis",
    y = TeX("$\\Pr(H_4)$"),
    fill = "Prior type"
  ) +
  scale_fill_manual(
    labels = c("Variant-specific prior", "Uniform prior"),
    values = c("#66CCEE", "#EE6677")
  ) +
  scale_x_discrete(labels = c(TeX("$H_3$"), TeX("$H_4$"))) +
  facet_wrap(~gene) +
  theme_jp()

simulation_plot <- gene_plot + (p1 / p2)

ggsave(
  snakemake@output[[1]],
  plot = simulation_plot,
  width = 10,
  height = 6
)

