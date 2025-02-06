
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(readr)
  library(bench)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(ggbeeswarm)
  library(patchwork)
})

source(here::here("code/plot-utils.R"))
devtools::load_all("~/coloc")

benchmark_path <- snakemake@output[["benchmark_plot_path"]]

data(coloc_test_data)
attach(coloc_test_data)

Q <- length(D1$beta)
w <- rpois(Q, lambda = 10) + 1
w <- w / sum(w)

bench_results <- bench::mark(
  coloc.abf(D1, D2),
  coloc.abf(D1, D2, prior_weights1 = w),
  check = FALSE
)

time <- bench_results$time
name_vec <- c(rep("Uniform", length(time[[1]])), rep("Variant-specific", length(time[[2]])))
time_vec <- unlist(lapply(time, function(x) as.numeric(str_sub(x, 1, -3))))

benchmark_plot_regular <- tibble(prior = name_vec, time = time_vec) |>
  ggplot(aes(prior, time, colour = prior)) +
  geom_quasirandom(method = "pseudorandom", size = 4) +
  coord_flip(ylim = c(3, 15)) +
  labs(
    y = "Time (ms)",
    x = "Method",
    title = "500 variants"
  ) +
  theme_jp_vgrid()

D1_large <- list(
  beta = rep(D1$beta, 4),
  varbeta = rep(D1$varbeta, 4),
  snp = paste0("SNP_", 1:(500 * 4)),
  position = as.character(1:(500 * 4)),
  type = "quant",
  sdY = 1.1
)

D2_large <- list(
  beta = rep(D2$beta, 4),
  varbeta = rep(D2$varbeta, 4),
  snp = paste0("SNP_", 1:(500 * 4)),
  position = as.character(1:(500 * 4)),
  type = "quant",
  sdY = 1.1
)

Q_large <- length(D1_large$beta)
w_large <- rpois(Q_large, lambda = 10) + 1
w_large <- w_large / sum(w_large)

bench_results <- bench::mark(
  coloc.abf(D1_large, D2_large),
  coloc.abf(D1_large, D2_large, prior_weights1 = w_large),
  check = FALSE
)

time <- bench_results$time
name_vec <- c(rep("Uniform", length(time[[1]])), rep("Variant-specific", length(time[[2]])))
time_vec <- unlist(lapply(time, function(x) as.numeric(str_sub(x, 1, -3))))

benchmark_plot_large <- tibble(prior = name_vec, time = time_vec) |>
  ggplot(aes(prior, time, colour = prior)) +
  geom_quasirandom(method = "pseudorandom", size = 4) +
  coord_flip(ylim = c(3, 15)) +
  labs(
    y = "Time (ms)",
    x = "Method",
    title = "2000 variants"
  ) +
  theme_jp_vgrid()

benchmark_plot <- benchmark_plot_regular / benchmark_plot_large +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.position = "top"
  )

ggsave(
  benchmark_path,
  plot = benchmark_plot,
  width = 14,
  height = 8
)
