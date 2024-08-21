# How to deal with duplicate variants? 

# renv::install("liftOver")

logsum <- function(x) {
    my.max <- max(x)                              ##take out the maximum value in log form
    my.res <- my.max + log(sum(exp(x - my.max ))) 
    return(my.res)
}

hg38tohg19 <- function(pos, chrom) {

  path <- system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
  chain <- rtracklayer::import.chain(path)
  current <- tibble::tibble(chr = paste0("chr", chrom), start = pos, end = pos) |>
    GenomicRanges::makeGRangesFromDataFrame()
  GenomeInfoDb::seqlevelsStyle(current) <- "UCSC"
  current_19 <- rtracklayer::liftOver(current, chain)

  pos <- as.vector(GenomicRanges::start(current_19))
  unmapped_ind <- vapply(pos, function(x) length(x) == 0, logical(1))
  n_unmapped <- sum(unmapped_ind)
  pos[unmapped_ind] <- -seq(1, n_unmapped)
  pos <- as.numeric(pos)
  dup_ind <- duplicated(pos)
  pos[dup_ind] <- -seq(n_unmapped + 1, n_unmapped + sum(dup_ind) + 1)
  pos
}

compute_polyfun_prior_weights <- function(pos, chrom, data, build = "hg38") {

  if (build == "hg38") {
    pos <- hg38tohg19(pos, chrom)
  }

  stopifnot(chrom %in% data$CHR)

  min_snpvar <- max(data$snpvar_bin) / 100
  snpvar_sum1 <- sum(data$snpvar_bin)
  data <- data |>
    mutate(snpvar_bin = if_else(snpvar_bin < min_snpvar, min_snpvar, snpvar_bin))
  snpvar_sum2 <- sum(data$snpvar_bin)
  data$snpvar_bin <- data$snpvar_bin * (snpvar_sum1 / snpvar_sum2)

  impute_value <- min(data$snpvar_bin)
  pos_data <- tibble::tibble(chr = chrom, pos) |>
    left_join(
      data |>
        filter(CHR == chrom),
      by = join_by(pos == BP)
    ) |>
    mutate(snpvar_bin = if_else(is.na(snpvar_bin), impute_value, snpvar_bin)) |>
    summarise(snpvar_bin = mean(snpvar_bin), .by = c(chr, pos))

  weights <- pos_data$snpvar_bin / sum(pos_data$snpvar_bin)
  weights
}

compute_abc_prior_weights <- function(pos, chrom, gene_name, abc_data, build = "hg38") {

  # FIXME: Could just save the abc_data with hg38 coordinates.
  if (build == "hg38") {
    pos <- hg38tohg19(pos, chrom)
  }

  abc_data <- abc_data |>
    filter(chr == paste0("chr", chrom)) |>
    filter(target_gene == gene_name)

  pos_data <- tibble::tibble(pos) |>
    dplyr::left_join(
      abc_data,
       by = dplyr::join_by(
        between(pos, start, end, bounds = "[]")
      ) 
    )

  pos_data <- pos_data |>
    # Take the median ABC score over celltypes/contexts.
    summarise(abc_score = median(abc_score), .by = pos) |>
    mutate(abc_score = if_else(is.na(abc_score), 0.0075, abc_score))

  weights <- pos_data$abc_score / sum(pos_data$abc_score)
  weights
}

compute_gnocchi_prior_weights <- function(pos, chrom, gnocchi_data) {

  if (!(all(colnames(gnocchi_data) %in% c("chromosome", "start_pos", "end_pos", "score")))) {
    stop("Wrong column names.")
  }

  gnocchi_data <- gnocchi_data |>
    dplyr::filter(chromosome == paste0("chr", chrom)) |>
    select(-chromosome)

  pos_data <- tibble::tibble(chromosome = paste0("chr", chrom), pos)
  overlap_match_pos <- pos_data |>
    dplyr::left_join(
      gnocchi_data,
      by = dplyr::join_by(
        # If bounds is not "()" dplyr will add rows.
        between(pos, start_pos, end_pos, bounds = "()")
      )
    )

  closest_start_pos <- pos_data |> 
    dplyr::left_join(
      tidyr::pivot_longer(gnocchi_data, -score, 
        values_to = "start_pos") |>
        filter(name == "start_pos"),
      by = dplyr::join_by(
        closest(pos <= start_pos)
      )
    )
  
  closest_end_pos <- pos_data |> 
    dplyr::left_join(
      tidyr::pivot_longer(gnocchi_data, -score, 
        values_to = "end_pos") |>
        filter(name == "end_pos"),
      by = dplyr::join_by(
        closest(pos >= end_pos)
      )
    )

  pos_data <- bind_cols(
    overlap_match_pos |>
      select(-c(start_pos, end_pos)) |>
      rename(score_overlap = score),
    closest_start_pos |>
      select(-c(chromosome, pos, name)) |>
      rename(score_start_pos = score), 
    closest_end_pos |>
      select(-c(chromosome, pos, name)) |>
      rename(score_end_pos = score), 
    ) |>
    mutate(
      end_pos_dist = abs(pos - end_pos), 
      start_pos_dist = abs(pos - start_pos)
    ) |>
    mutate(score_closest = if_else(end_pos_dist < start_pos_dist, score_end_pos, score_start_pos)) |>
    mutate(
      score_closest = case_when(
        is.na(score_start_pos) ~ score_end_pos,
        is.na(score_end_pos) ~ score_start_pos,
        .default = score_closest
    )) |>
    mutate(score = if_else(is.na(score_overlap), score_closest, score_overlap), 
          .keep = "unused")

  if (any(is.na(pos_data$score))) {
    stop("Some SNPs do not lie in gnocchi annorated regions")
  }

  lsum <- logsum(pos_data$score)
  weights <- exp(pos_data$score - lsum)

  weights 
}

compute_eqtl_tss_dist_prior_weights <- function(pos, tss, density_data) {
  rel <- pos - tss
  x <- density_data$x
  y <- density_data$y
  closest <- numeric(length(pos))
  for (i in seq_along(pos)) {
    z <- which(abs(x - rel[[i]]) == min(abs(x - rel[[i]])))
    if (length(z) > 1) {
      closest[[i]] <- z[[1]]
    } else {
      closest[[i]] <- z
    }
  }
  # Normalise the values to sum to 1.
  out <- y[closest] / sum(y[closest])
  out
}

compute_rand_prior_weights <- function(pos) {
  as.vector(extraDistr::rdirichlet(1, rep(1, length(pos))))
}
