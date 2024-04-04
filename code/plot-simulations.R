
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(latex2exp)
  library(stringr)
  library(tools)
})

source("code/plot-utils.R")

devtools::load_all("~/coloc")

simulation_data <- lapply(snakemake@input, function(x) {
  gene_name <- file_path_sans_ext(x) |>
    basename() |>
    str_extract("([^-]+$)")

  data <- read_rds(x) |>
    mutate(gene = gene_name)
}) |>
  bind_rows()

simulation_plot <- simulation_data |>
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

ggsave(
  snakemake@output[[1]],
  plot = simulation_plot,
  width = 8,
  height = 6
)

