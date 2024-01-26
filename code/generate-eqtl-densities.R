
source(here::here("renv/activate.R"))

library(dplyr)
library(biomaRt)
library(janitor)

density_to_tibble <- function(dens) {
  tibble(x = dens$x, y = dens$y)
}

hg19_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                        "https://feb2014.archive.ensembl.org")

hg19_tss_data <- getBM(
  attributes = c("ensembl_gene_id", "transcript_start", "transcript_end", "strand"),
  filters = "ensembl_gene_id",
  values = "",
  mart = hg19_mart
) |>
  as_tibble() |>
  mutate(transcription_start_site = if_else(strand == 1, transcript_start, transcript_end)) |>
  select(transcription_start_site, ensembl_gene_id, strand) |>
  rename(tss = transcription_start_site)

onek1k_data <- read_tsv("data/onek1k.tsv",
                        col_select = c(gene_id, pos, p_value, cell_type, round),
                        name_repair = make_clean_names,
                        show_col_types = FALSE, progress = FALSE)

onek1k_cd4nc_round_1 <- onek1k_data |>
  filter(p_value < 5 * 10^-8) |>
  # Filter to the most abundant cell type.
  filter(cell_type == "CD4 Naive/Central memory T cell") |>
  filter(round == 1) |>
  left_join(hg19_tss_data, by = join_by(gene_id == ensembl_gene_id),
    relationship = "many-to-many"
  ) |>
  # FIXME: Is this right?
  mutate(tss_distance = if_else(strand == 1, pos - tss, tss - pos)) |>
  filter(abs(tss_distance) <= 1e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 5) |>
  density_to_tibble()

onek1k_cd4nc_round_2 <- onek1k_data |>
  filter(p_value < 5 * 10^-8) |>
  # Filter to the most abundant cell type.
  filter(cell_type == "CD4 Naive/Central memory T cell") |>
  filter(round == 2) |>
  left_join(hg19_tss_data, by = join_by(gene_id == ensembl_gene_id)) |>
  # fixme: is this right?
  mutate(tss_distance = if_else(strand == 1, pos - tss, tss - pos)) |>
  filter(abs(tss_distance) <= 1e5) |>
  pull(tss_distance) |>
  density(bw = "sj", cut = 0, adjust = 5) |>
  density_to_tibble()

onek1k_cd4nc_round_3 <- onek1k_data |>
  filter(p_value < 5 * 10^-8) |>
  # Filter to the most abundant cell type.
  filter(cell_type == "CD4 Naive/Central memory T cell") |>
  filter(round >= 3) |>
  left_join(hg19_tss_data, by = join_by(gene_id == ensembl_gene_id),
    relationship = "many-to-many"
  ) |>
  # FIXME: Is this right?
  mutate(tss_distance = if_else(strand == 1, pos - tss, tss - pos)) |>
  filter(abs(tss_distance) <= 1e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 5) |>
  density_to_tibble()

eqtlgen_data <- read_tsv("data/eqtlgen.txt",
                         col_select = c(gene_symbol, snp_pos, gene_pos, pvalue),
                         name_repair = make_clean_names,
                         show_col_types = FALSE, progress = FALSE)

eqtlgen <- raw_eqtlgen_data |>
  filter(pvalue < 5 * 10^-8) |>
  group_by(gene_symbol) |>
  arrange(pvalue) |>
  filter(row_number() == 1) |>
  mutate(tss_distance = snp_pos - gene_pos) |>
  ungroup() |>
  filter(abs(tss_distance) <= 3e5) |>
  pull(tss_distance) |>
  density(bw = "SJ", cut = 0, adjust = 5) |>
  density_to_tibble()

saveRDS(eqtlgen, "output/densities/eqtlgen.rds")
saveRDS(onek1k_cd4nc_round_1, "output/densities/onek1k_cd4nc_round_1.rds")
saveRDS(onek1k_cd4nc_round_2, "output/densities/onek1k_cd4nc_round_2.rds")
saveRDS(onek1k_cd4nc_round_3, "output/densities/onek1k_cd4nc_round_3.rds")
