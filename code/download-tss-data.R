
source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

hg19_tss_data_path <- snakemake@output[["hg19_tss_data_path"]]
hg38_tss_data_path <- snakemake@output[["hg38_tss_data_path"]]

# HG38 data

hg38_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

hg38_tss_data <- getBM(
  attributes = c("ensembl_gene_id", "transcription_start_site"),
  filters = "ensembl_gene_id",
  values = "",
  mart = hg38_mart
) |>
  as_tibble() |>
  rename(tss = transcription_start_site) |>
  group_by(ensembl_gene_id) |>
  # NOTE: For genes with multiple TSSs we take the median of the listed TSSs.
  summarise(tss = median(tss))

saveRDS(hg38_tss_data, hg38_tss_data_path)

# HG19 data

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
  select(transcription_start_site, ensembl_gene_id) |>
  rename(tss = transcription_start_site) |>
  group_by(ensembl_gene_id) |>
  # NOTE: For genes with multiple TSSs we take the median of the listed TSSs.
  summarise(tss = median(tss))

saveRDS(hg19_tss_data, hg19_tss_data_path)
