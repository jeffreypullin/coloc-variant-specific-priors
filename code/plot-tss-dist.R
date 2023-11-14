library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(janitor)
library(biomaRt)
library(glue)
library(tools)
library(forcats)

eqtl_catalogue_metadata <- read_tsv(
  here::here("data", "eqtl-catalogue", "tabix_ftp_paths.tsv")
)

process_file_eqtl_catalogue <- function(path, mart) {

  raw_data <- read_tsv(path, col_select = c(gene_id, position, pvalue),
                       show_col_types = FALSE, progress = FALSE)

  gene_snp_data <- raw_data |>
    group_by(gene_id) |>
    arrange(pvalue) |>
    filter(row_number() == 1) |>
    ungroup()

  tss_data <- getBM(
    attributes = c("ensembl_gene_id", "transcription_start_site"),
    filters = "ensembl_gene_id",
    values = gene_snp_data$gene_id,
    mart = mart
  ) |>
    as_tibble() |>
    rename(tss = transcription_start_site) |>
    group_by(ensembl_gene_id) |>
    # NOTE: For genes with multiple TSSs we take the median of the listed TSSs.
    summarise(tss = median(tss))

  gene_snp_data |>
    left_join(tss_data, by = join_by(gene_id == ensembl_gene_id)) |>
    # NOTE: Remove genes that do not have a corresponding TSS.
    filter(!is.na(tss)) |>
    mutate(
      tss_distance = position - tss,
      abs_tss_distance = abs(tss_distance),
      log10_abs_tss_distance = log10(abs_tss_distance)
    ) |>
    filter(log10_abs_tss_distance > 0) |>
    pull(log10_abs_tss_distance)
}

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

eqtl_catalogue_data <- tibble(
  path = list.files(here::here("data", "eqtl-catalogue", "processed-sumstats"), full.names = TRUE)
) |>
  mutate(file = gsub("(.+).cc.tsv", "\\1", basename(path))) |>
  rowwise() |>
  mutate(log10_abs_tss_distance = list(process_file_eqtl_catalogue(path, mart))) |>
  ungroup() |>
  unnest(cols = log10_abs_tss_distance)

eqtl_catalogue_data |>
  ggplot(aes(log10_abs_tss_distance, file)) +
  geom_density_ridges() +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    y = "",
    x = "Absolute Log(10) distance from variant to TSS",
    title = "eQTL distance distribution - eQTL catalogue"
  ) +
  theme_bw() +
  theme(rect = element_rect(linewidth = 0))

eqtl_catalogue_data |>
  left_join(eqtl_catalogue_metadata, join_by(file == dataset_id)) |>
  mutate(y = 10^log10_abs_tss_distance) |>
  ggplot(aes(file, y, fill = study_label)) +
  geom_boxplot(outlier.alpha = 0.1) +
  coord_flip(ylim = c(0, 150000)) +
  labs(
    x = "Dataset",
    y = "Distance from the TSS"
  )
  theme_bw()

eqtl_catalogue_data |>
  left_join(eqtl_catalogue_metadata, join_by(file == dataset_id)) |>
  mutate(y = 10^log10_abs_tss_distance) |>
  ggplot(aes(study_label, y, fill = study_label)) +
  geom_boxplot(outlier.alpha = 0.1) +
  coord_flip(ylim = c(0, 150000)) +
  labs(
    x = "Dataset",
    y = "Distance from the TSS"
  ) +
  theme_bw()

eqtl_catalogue_data |>
  left_join(eqtl_catalogue_metadata, join_by(file == dataset_id)) |>
  mutate(y = 10^log10_abs_tss_distance) |>
  summarise(y = median(y), .by = c(file, sample_size)) |>
  ggplot(aes(sample_size, y)) +
  geom_point(size = 2) +
  labs(
    x = "Study sample size",
    y = "Median distance to TSS"
  ) +
  coord_cartesian(ylim = c(10000, 25000), xlim = c(0, 700)) +
  theme_bw()

# GTEx data

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

gtex_data <- tibble(path = list.files(here::here("data", "gtex-v8"), full.names = TRUE)) |>
  rowwise() |>
  mutate(tissue = gsub("_", " ", unlist(strsplit(basename(path), "[.]"))[[1]])) |>
  mutate(abs_tss_distance = list(process_file_gtex(path))) |>
  ungroup() |>
  unnest(cols = abs_tss_distance)

gtex_data |>
  mutate(tissue = fct_rev(factor(tissue))) |>
  ggplot(aes(tissue, abs_tss_distance)) +
  geom_boxplot(fill = "grey", outlier.alpha = 0.1) +
  labs(
    x = "Tissue",
    y = "Distance to TSS"
  ) +
  coord_flip() +
  theme_bw()

gtex_data |>
  ggplot(aes(abs_tss_distance)) +
  geom_histogram(binwidth = 1000, colour = "black", fill = "grey") +
  coord_cartesian(xlim = c(0, 500000)) +
  labs(
    y = "Count",
    x = "Absolute distance from variant to TSS",
    title = "eQTL distance distribution in GTEx v8"
  ) +
  theme_bw()

# eQTLGen data

raw_eqtlgen_data <- read_tsv(here::here("data", "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"),
                             col_select = c(gene_symbol, snp_pos, gene_pos, pvalue),
                             name_repair = make_clean_names)

eqtlgen_data <- raw_eqtlgen_data |>
  filter(pvalue < 5 * 10^-8) |>
  group_by(gene_symbol) |>
  arrange(pvalue) |>
  filter(row_number() == 1) |>
  mutate(abs_tss_distance = abs(snp_pos - gene_pos)) |>
  mutate(log10_abs_tss_distance = log10(abs_tss_distance)) |>
  filter(log10_abs_tss_distance > 0) |>
  ungroup()

eqtlgen_data |>
  ggplot(aes(log10_abs_tss_distance)) +
  geom_histogram(binwidth = 0.1, colour = "black", fill = "grey") +
  labs(
    y = "Count",
    x = "Absolute Log(10) distance from variant to TSS",
    title = "eQTL distance distribution in eQTLGen"
  ) +
  theme_bw()

bind_rows(
  gtex_data |>
    dplyr::select(log10_abs_tss_distance) |>
    mutate(study = "gtex_v8"),
  eqtlgen_data |>
    dplyr::select(log10_abs_tss_distance) |>
    mutate(study = "eqtlgen"),
  eqtl_catalogue_data |>
    dplyr::select(log10_abs_tss_distance) |>
    mutate(study = "eqtl_catalogue")
) |>
  ggplot(aes(study, log10_abs_tss_distance)) +
  geom_violin(fill = "grey") +
  geom_boxplot() +
  labs(
    x = "Study",
    y = "Absolute Log(10) distance from variant to TSS",
    title = "eQTL distance distribution across studies"
  ) +
  coord_flip() +
  theme_bw()
