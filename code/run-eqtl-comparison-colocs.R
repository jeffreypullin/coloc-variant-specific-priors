library(readr)
library(janitor)
library(dplyr)
library(patchwork)
library(tidyr)
library(purrr)
library(stringr)

hg38_data <- readRDS("/home/jp2045/coloc-estimated-eqtl-priors/data/tss-data/hg38-tss-data.rds")
load("/rds/project/cew54/rds-cew54-basis/People/CHRIS/eqtl_colocs_for_jeffrey.RData")
devtools::load_all("~/coloc")

set.seed(42)

lookup_tss <- function(trait) {
  gene_id <- unlist(strsplit(trait, ".", fixed = TRUE))[[2]]
  hg38_data |>
    filter(ensembl_gene_id == gene_id) |>
    pull(tss)
}

eqtl_method_comparison_data <- COLOC |>
  # The BLUEPRINT Bayes factors are empty.
  filter(!(substr(trait1, 1, 9) == "BLUEPRINT" | substr(trait2, 1, 9) == "BLUEPRINT")) |>
  as_tibble() |>
  mutate(
    gene1 = grepl("ENSG", trait1, fixed = TRUE),
    gene2 = grepl("ENSG", trait2, fixed = TRUE)
  ) |>
  filter(gene1 & gene2) |>
  slice_sample(n = 500) |>
  rowwise() |>
  mutate(tss1 = list(lookup_tss(trait1)), tss2 = list(lookup_tss(trait2))) |>
  filter(!(is.null(tss1) & !is.null(tss1))) |>
  select(trait1, trait2, block, tss1, tss2, hit1.margz, hit2.margz) |>
  mutate(bf1 = list(BF[[block]][[trait1]]), bf2 = list(BF[[block]][[trait2]])) |>
  mutate(
    coloc_unif = list(coloc.bf_bf(bf1, bf2)),
    coloc_non_unif_tss1 = list(coloc.bf_bf(bf1, bf2, tss1 = tss1)),
    coloc_non_unif_tss2 = list(coloc.bf_bf(bf1, bf2, tss2 = tss2)),
    coloc_non_unif_both = list(coloc.bf_bf(bf1, bf2, tss1 = tss1, tss2 = tss2)),
  ) |>
  mutate(
    h4_unif = median(coloc_unif$summary$PP.H4.abf),
    h4_non_unif_tss1 = median(coloc_non_unif_tss1$summary$PP.H4.abf),
    h4_non_unif_tss2 = median(coloc_non_unif_tss2$summary$PP.H4.abf),
    h4_non_unif_both = median(coloc_non_unif_both$summary$PP.H4.abf),
  ) |>
  filter(!is.na(h4_unif)) |>
  mutate(
    h4_diff_tss1 = h4_non_unif_tss1 - h4_unif,
    h4_diff_tss2 = h4_non_unif_tss2 - h4_unif,
    h4_diff_both = h4_non_unif_both - h4_unif
  )

saveRDS(eqtl_method_comparison_data,
        here::here("output", "data", "eqtl-method-comparison.rds"))


reduced_blocks <- unique(eqtl_method_comparison_data$block)

tmp <- eqtl_method_comparison_data |>
  ungroup() |>
  summarise(traits = list(c(trait1, trait2)), .by = "block")

reduced_blocks_trait_list <- tmp$traits
names(reduced_blocks_trait_list) <- tmp$block

reduced_bf_list <- BF[reduced_blocks]
for (i in 1:length(reduced_bf_list)) {
  block <- names(reduced_bf_list)[[i]]
  reduced_bf_list[[i]] <- reduced_bf_list[[i]][reduced_blocks_trait_list[[block]]]
}

saveRDS(reduced_bf_list, here::here("output", "data", "reduced-bf-list.rds"))
