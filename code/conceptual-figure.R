
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("code/plot-utils.R")

conceptual_figure <- tibble(
    hyp = c("H4", "H3", "H4", "H3"),
    val = c(0.6, 0.35, 0.87, 0.1),
    type = c("Uniform", "Uniform", "Variant-specific", "Variant-specific")
  ) |>
  ggplot(aes(hyp, val, fill = factor(hyp))) +
  geom_col(col = "black", alpha = 0.8) +
  geom_hline(yintercept = 0.8) +
  labs(
    x = "Hypothesis",
    y = "Probability",
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c("deeppink1", "darkorange1")) +
  facet_wrap(~type) +
  theme_jp() +
  theme(
    legend.position = "none",
    axis.text = element_text(family = "Helvetica", size = 36),
    axis.title = element_text(family = "Helvetica", size = 45),
    strip.text.x = element_text(family = "Helvetica", size = 50),
  )

ggsave(
  "output/figures/conceptual-figure-panel.pdf",
  plot = conceptual_figure,
  width = 20,
  height = 6
)
