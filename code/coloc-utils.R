
coloc_to_tibble <- function(coloc_out, suffix) {

  coloc_summary <- coloc_out$summary
  # SuSiE output.
  if (inherits(coloc_summary, "data.frame")) {
    out <- as.data.frame(coloc_summary)
  } else {
    out <- t(as.data.frame(coloc_summary))
  }
  colnames(out) <- paste0(colnames(out), "_", suffix)
  out
}

eqtl_catalouge_colnames <- c(
  "molecular_trait_id", "chromosome", "position", "ref",
  "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
  "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id",
  "median_tpm", "rsid"
)

prepare_coloc_dataset <- function(data) {
  data |>
    select(-rsid) |>
    distinct() |>
    mutate(id = paste(chromosome, position, sep = ":")) |>
    group_by(id) |>
    mutate(row_count = n()) |>
    ungroup() |>
    filter(row_count == 1) |>
    filter(!is.nan(se)) |>
    filter(!is.na(se)) |>
    filter(maf > 0 & maf < 1) |>
    mutate(
      maf = as.numeric(maf),
      beta = as.numeric(beta),
      se = as.numeric(se)
  )
}

filter_qtl_dataset <- function(data, trait_id, chrom, start_pos, end_pos) {
  data |>
    filter(molecular_trait_id == !!trait_id) |>
    filter(chromosome == !!chrom) |>
    filter(
      position >= !!start_pos &
      position <= !!end_pos
    )
}
