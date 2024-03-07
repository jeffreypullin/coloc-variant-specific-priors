
# Packages.

source(here::here("renv/activate.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(simGWAS)
  devtools::load_all("~/coloc")
})

# Functions.

simulate_dataset <- function(cv) {

  stopifnot(length(cv) == 1)
  g <- max(rnorm(100, sd = 0.05))

  n_0 <- 2000
  n_1 <- 2000
  n_rep <- 2000

  accept <- FALSE
  while (!accept) {

    print("Simulating...")
    z <- simulated_z_score(
      N0 = n_0,
      N1 = n_1,
      snps = snps,
      W = cv,
      gamma.W = g,
      freq = freq,
      nrep = n_rep
    )
    varbeta <- simulated_vbeta(
      N0 = n_0,
      N1 = n_1,
      snps = snps,
      W = cv,
      gamma.W = g,
      freq = freq,
      nrep = n_rep
    )

    beta <- z * sqrt(varbeta)
    accept_ind <- which(abs(z[, match(cv, snps)]) > 4.417)
    accept <- length(accept_ind) > 0
  }
  ind <- accept_ind[[1]]

  list(
    snp = snps,
    position = pos,
    beta = beta[ind, ],
    varbeta = varbeta[ind, ],
    z = z[ind, ],
    MAF = maf,
    LD = ld,
    type = "cc",
    N = n_0 + n_1,
    s = 0.5
  )
}

simulate_dataset_pair <- function(hyp, snps, cv1) {

  stopifnot(cv1 %in% snps)

  if (hyp == "h3") {
    cv1_ind <- which(snps == cv1)
    cv2 <- snps[sample(setdiff(seq_along(snps), cv1_ind), 1)]
  } else if (hyp == "h4") {
    cv2 <- cv1
  }

  dataset_1 <- simulate_dataset(cv1)
  dataset_2 <- simulate_dataset(cv2)

  list(dataset_1 = dataset_1, dataset_2 = dataset_2)
}

# Script.

set.seed(13022024)

legend_file <- snakemake@input[["leg_file"]]
haps_file <- snakemake@input[["haps_file"]]
sim_data_file <- snakemake@output[["sim_file"]]

legend <- read_delim(legend_file, show_col_types = FALSE)
haps <- matrix(scan(haps_file, what = 0), ncol = nrow(legend))
colnames(haps) <- legend$ID

if (nrow(legend) > 500) {
  sample_ind <- sort(sample(seq_len(nrow(legend)), 500))
  haps <- haps[, sample_ind]
  legend <- legend[sample_ind, ]
}

ld <- cor(haps)
ind <- which(apply(ld, 2, function(x) sum(x >= 0.7)) >= 20)
if (length(ind) > 0) {
  haps <- haps[, ind]
  legend <- legend[ind, ]
}

freq <- as.data.frame((haps + 1))
freq$Probability <- 1 / nrow(freq)

ld <- cor(haps)
snps <- colnames(haps)
maf <- pmin(colMeans(haps), 1 - colMeans(haps))
pos <- vapply(strsplit(snps, "-"), \(x) as.numeric(x[[2]]), numeric(1))

n_sims <- 100
prior_weights <- compute_eqtl_tss_dist_prior_weights(
  pos,
  round(median(pos), 0),
  read_rds("output/densities/eqtlgen.rds")
)
hyps <- sample(c("h3", "h4"), size = n_sims, replace = TRUE, prob = c(0.5, 0.5))

pp_h4 <- numeric(n_sims)
sim_out <- list()
for (i in seq_len(n_sims)) {

  cv1_ind <- sample(seq_along(snps), 1, prob = prior_weights)

  cv1 <- snps[cv1_ind]
  dataset_pair <- simulate_dataset_pair(hyps[[i]], snps, cv1)

  non_unif_coloc_out <- coloc.abf(
    dataset_pair$dataset_1,
    dataset_pair$dataset_2,
    prior_weights1 = prior_weights
  )$summary

  unif_coloc_out <- coloc.abf(
    dataset_pair$dataset_1,
    dataset_pair$dataset_2,
    prior_weights1 = NULL
  )$summary

  pp_h4 <- c(
    unif_coloc_out["PP.H4.abf"],
    non_unif_coloc_out["PP.H4.abf"]
  )
  prior <- c("unif", "non_unif")

  sim_out[[i]] <- tibble(hyp = hyps[[i]], pp_h4 = pp_h4, prior)
}

sim_data <- bind_rows(!!!sim_out)

saveRDS(sim_data, sim_data_file)
