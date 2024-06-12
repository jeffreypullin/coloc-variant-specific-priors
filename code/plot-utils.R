
theme_jp <- function() {

  font <- "Helvetica"

  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = font, size = 24, color = "#222222"
    ),
    plot.subtitle = ggplot2::element_text(
      family = font, size = 20, margin = ggplot2::margin(0, 0, 10, 0)
    ),
    plot.caption = ggplot2::element_blank(),
    plot.title.position = "plot",
    legend.position = "top",
    legend.box.margin = ggplot2::margin(t = -5),
    legend.text.align = 0,
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family = font, size = 14, color = "#222222"),
    axis.title = ggplot2::element_text(family = font, size = 14, color = "#222222"),
    axis.text = ggplot2::element_text(family = font, size = 14, color = "#222222"),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"),
    panel.grid.minor.y = ggplot2::element_line(color = "#cdcdcd"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(family = font, size = 18)
  )
}

theme_jp_vgrid <- function() {
  theme_jp() %+replace%
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "#cbcbcb"),
      panel.grid.minor.x = ggplot2::element_line(color = "#cdcdcd"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
    )
}

prior_method_lookup <- c(
  "abc_score" = "ABC score",
  "eqtlgen" = "eQTLGen",
  "gnocchi" = "Gnocchi",
  "onek1k_r1" = "OneK1K (R1)",
  "onek1k_r2" = "OneK1k (R2+)",
  "polyfun" = "PolyFun",
  "unif" = "Uniform"
)

dataset_lookup <- c(
  "eqtl" = "eQTL",
  "pqtl" = "pQTL",
  "pqtl_eqtl" = "Both"
)
