# End-to-end pipeline for v0.6.
# Designed for RStudio / Quarto. Nothing auto-runs when this file is sourced.
#
# Primary v0.6 change:
#   finite i.i.d. Zipf, sticky Markov Zipf and canonical SSR all use V = 5000
#   and share the same theoretical stationary marginal p_j = 1/(j H_V).
#   SSR is therefore not calibrated to the fitted exponent. Calibration is
#   reserved for the latent-mixture construction and for the random-segmentation
#   / Simon compatibility examples.

.load_zipf_benchmark_code <- function() {
  needed <- c(
    "simulate_maxent_zipf", "simulate_markov_zipf", "simulate_ssr",
    "zipf_summary", "zipf_mle_summary", "calibrate_benchmark",
    "calibrate_latent_stable", "run_strict_replications",
    "pairwise_auc_bootstrap"
  )
  if (!all(vapply(needed, exists, logical(1), mode = "function"))) {
    source("R/generators.R")
    source("R/diagnostics.R")
    source("R/calibration.R")
  }
  invisible(TRUE)
}

.strict_alpha_status <- function(results) {
  rows <- lapply(split(seq_len(nrow(results)), results$generator), function(idx) {
    d <- results[idx, , drop = FALSE]
    a_ols <- d$alpha[is.finite(d$alpha)]
    a_mle <- d$alpha_mle[is.finite(d$alpha_mle)]
    data.frame(
      generator = d$generator[1L],
      mean_alpha_ols = mean(a_ols),
      sd_alpha_ols = if (length(a_ols) >= 2L) stats::sd(a_ols) else NA_real_,
      abs_bias_ols = abs(mean(a_ols) - 1),
      mean_alpha_mle = mean(a_mle),
      sd_alpha_mle = if (length(a_mle) >= 2L) stats::sd(a_mle) else NA_real_,
      abs_bias_mle = abs(mean(a_mle) - 1),
      n = nrow(d),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

zipf_probabilities <- function(V = 5000L, exponent = 1) {
  r <- seq_len(as.integer(V))
  w <- r^(-as.numeric(exponent))
  w / sum(w)
}

markov_tau_block <- function(rho, block_size = 250L) {
  rho <- as.numeric(rho)
  L <- as.integer(block_size)
  stopifnot(rho >= 0, rho < 1, L >= 2L)
  h <- seq_len(L - 1L)
  1 + 2 * sum((1 - h / L) * rho^h)
}

markov_effective_n <- function(n, rho, block_size = 250L) {
  as.numeric(n) / markov_tau_block(rho, block_size)
}

expected_vocabulary_iid <- function(n, V = 5000L, exponent = 1) {
  p <- zipf_probabilities(V, exponent)
  sum(1 - (1 - p)^as.integer(n))
}

expected_vocabulary_markov <- function(
    n, rho, V = 5000L, exponent = 1, block_size = 250L) {
  n <- as.integer(n)
  L <- as.integer(block_size)
  stopifnot(n %% L == 0L)
  B <- n %/% L
  p <- zipf_probabilities(V, exponent)
  absent_one_block <- (1 - p) * (1 - (1 - rho) * p)^(L - 1L)
  sum(1 - absent_one_block^B)
}

.run_rmin_sensitivity <- function(
    n_rep,
    first_seed,
    n_tokens,
    parameters,
    r_min_values = c(1L, 5L, 10L, 20L, 50L, 100L),
    r_max = 500L,
    min_count = 5L,
    save_path = NULL) {
  rows <- list()
  k <- 1L
  for (b in seq_len(n_rep)) {
    seed <- first_seed + b - 1L
    sims <- simulate_strict_benchmark_once(
      seed = seed, n_tokens = n_tokens, parameters = parameters
    )
    for (j in seq_along(sims)) {
      d <- sims[[j]]
      for (rmin in r_min_values) {
        z <- zipf_summary(d, r_min = rmin, r_max = r_max, min_count = min_count)
        zm <- zipf_mle_summary(d, r_min = rmin, r_max = r_max, min_count = min_count)
        rows[[k]] <- data.frame(
          generator = unique(d$generator)[1L],
          replicate = b,
          seed = seed,
          r_min = as.integer(rmin),
          r_max = as.integer(r_max),
          alpha = z$alpha,
          alpha_mle = zm$alpha_mle,
          n_fit = z$n_fit,
          effective_r_max = z$effective_r_max,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, save_path, row.names = FALSE)
  }
  out
}

.run_markov_rho_sensitivity <- function(
    n_rep,
    first_seed,
    n_tokens,
    vocabulary_size = 5000L,
    block_size = 250L,
    rho_values = c(0, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80),
    r_min = 10L,
    r_max = 500L,
    min_count = 5L,
    save_path = NULL) {
  rows <- list()
  k <- 1L
  for (rho in rho_values) {
    neff <- markov_effective_n(n_tokens, rho, block_size)
    ev <- expected_vocabulary_markov(
      n_tokens, rho, vocabulary_size, 1, block_size
    )
    for (b in seq_len(n_rep)) {
      seed <- first_seed + 10000L * k + b - 1L
      dm <- simulate_markov_zipf(
        n_tokens = n_tokens,
        vocabulary_size = vocabulary_size,
        exponent = 1,
        persistence = rho,
        block_size = block_size,
        seed = seed
      )
      zm <- zipf_summary(dm, r_min = r_min, r_max = r_max, min_count = min_count)
      zmm <- zipf_mle_summary(dm, r_min = r_min, r_max = r_max, min_count = min_count)

      di <- simulate_maxent_zipf(
        n_tokens = max(1000L, as.integer(round(neff))),
        vocabulary_size = vocabulary_size,
        exponent = 1,
        block_size = block_size,
        seed = seed + 500000L
      )
      zi <- zipf_summary(di, r_min = r_min, r_max = r_max, min_count = min_count)
      zim <- zipf_mle_summary(di, r_min = r_min, r_max = r_max, min_count = min_count)

      rows[[length(rows) + 1L]] <- data.frame(
        rho = rho,
        n_eff = neff,
        replicate = b,
        markov_alpha = zm$alpha,
        iid_neff_alpha = zi$alpha,
        markov_alpha_mle = zmm$alpha_mle,
        iid_neff_alpha_mle = zim$alpha_mle,
        markov_vocabulary = zm$vocabulary,
        expected_markov_vocabulary = ev,
        stringsAsFactors = FALSE
      )
    }
    k <- k + 1L
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, save_path, row.names = FALSE)
  }
  out
}

.run_ssr_boundary_sensitivity <- function(
    n_rep,
    first_seed,
    n_tokens,
    vocabulary_size = 5000L,
    block_size = 250L,
    n_surrogates = 20L,
    top_k = 50L,
    save_path = NULL) {
  rows <- vector("list", n_rep)
  for (b in seq_len(n_rep)) {
    seed <- first_seed + b - 1L
    d_nat <- simulate_ssr(
      n_tokens = n_tokens,
      vocabulary_size = vocabulary_size,
      initial_distribution = "stationary",
      seed = seed
    )
    d_fix <- d_nat
    d_fix$sequence_id <- ceiling(seq_len(nrow(d_fix)) / block_size)

    nat_nmi <- transition_information_excess(
      d_nat, top_k = top_k, n_surrogates = n_surrogates,
      seed = 700000L + seed
    )
    fix_nmi <- transition_information_excess(
      d_fix, top_k = top_k, n_surrogates = n_surrogates,
      seed = 800000L + seed
    )
    nat_d <- transition_direction(d_nat)
    fix_d <- transition_direction(d_fix)

    rows[[b]] <- data.frame(
      replicate = b,
      seed = seed,
      transition_nmi_excess_natural = nat_nmi$transition_nmi_excess,
      transition_nmi_excess_fixed = fix_nmi$transition_nmi_excess,
      D_natural = nat_d$toward_more_frequent,
      D_fixed = fix_d$toward_more_frequent,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, save_path, row.names = FALSE)
  }
  out
}

.paper_outputs <- function(results_dir = "results") {
  file.path(results_dir, c(
    "benchmark_strict_replications.rds",
    "benchmark_mimicking_replications.rds",
    "diagnostic_pairwise_auc_bootstrap_strict.csv",
    "nmi_topk_sensitivity_strict.csv",
    "alpha_sensitivity_strict_replications.csv",
    "alpha_sensitivity_mimicking_replications.csv",
    "alpha_rmin_sensitivity_strict.csv",
    "markov_rho_sensitivity.csv",
    "ssr_boundary_sensitivity.csv",
    "calibration_exponent_summary.csv",
    "strict_alpha_status.csv"
  ))
}

paper_outputs_exist <- function(results_dir = "results") {
  all(file.exists(.paper_outputs(results_dir)))
}

run_all <- function(
    n_rep = 100L,
    n_cal_seeds = 10L,
    n_latent_cal_seeds = 50L,
    n_surrogates = 20L,
    heaps_window = 7L,
    n_boot = 2000L,
    nmi_top_k_values = c(25L, 50L, 100L),
    primary_nmi_top_k = 50L,
    n_tokens = 100000L,
    zipf_r_max = 500L,
    vocabulary_size = 5000L,
    block_size = 250L,
    markov_persistence = 0.35,
    rho_values = c(0, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80),
    r_min_values = c(1L, 5L, 10L, 20L, 50L, 100L),
    results_dir = "results") {
  .load_zipf_benchmark_code()

  n_rep <- as.integer(n_rep)
  n_cal_seeds <- as.integer(n_cal_seeds)
  n_latent_cal_seeds <- as.integer(n_latent_cal_seeds)
  n_surrogates <- as.integer(n_surrogates)
  heaps_window <- as.integer(heaps_window)
  n_boot <- as.integer(n_boot)
  nmi_top_k_values <- sort(unique(as.integer(nmi_top_k_values)))
  primary_nmi_top_k <- as.integer(primary_nmi_top_k)
  vocabulary_size <- as.integer(vocabulary_size)
  block_size <- as.integer(block_size)
  r_min_values <- sort(unique(as.integer(r_min_values)))
  rho_values <- sort(unique(as.numeric(rho_values)))

  stopifnot(n_rep >= 2L, n_cal_seeds >= 2L, n_latent_cal_seeds >= 5L)
  stopifnot(n_surrogates >= 1L, n_boot >= 100L)
  stopifnot(primary_nmi_top_k %in% nmi_top_k_values)
  stopifnot(n_tokens %% block_size == 0L)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  core_cal_seeds <- 1001L + seq_len(n_cal_seeds) - 1L
  latent_cal_seeds <- 1001L + seq_len(n_latent_cal_seeds) - 1L
  latent_screening_seeds <- 1001L + seq_len(min(10L, n_latent_cal_seeds)) - 1L

  cat("\n=== 1. Calibration of tunable non-exact constructions ===\n")
  cat("Random segmentation and Simon use ", length(core_cal_seeds),
      " calibration seeds; latent mixture uses ", length(latent_cal_seeds),
      " for final selection. i.i.d. Zipf, Markov and SSR are not calibrated.\n",
      sep = "")

  cal_exp_core <- calibrate_benchmark(
    mode = "exponent",
    seeds = core_cal_seeds,
    n_tokens = n_tokens,
    alpha_tolerance = 0.01,
    refine_exponent = TRUE,
    zipf_r_max = zipf_r_max,
    models = c("random", "simon")
  )
  cal_exp_core$seed_count_by_model <- stats::setNames(
    rep(length(core_cal_seeds), length(cal_exp_core$tables)),
    names(cal_exp_core$tables)
  )

  cal_exp_latent <- calibrate_latent_stable(
    mode = "exponent",
    screening_seeds = latent_screening_seeds,
    final_seeds = latent_cal_seeds,
    n_tokens = n_tokens,
    alpha_tolerance = 0.005,
    zipf_r_max = zipf_r_max
  )

  cal_exp <- combine_calibration_objects(
    cal_exp_core, cal_exp_latent,
    save_path = file.path(results_dir, "calibration_exponent.rds")
  )
  cal_exp$version <- "0.6"
  saveRDS(cal_exp, file.path(results_dir, "calibration_exponent.rds"))
  print(write_calibration_outputs(cal_exp, results_dir))

  strict_parameters <- strict_benchmark_parameters(
    exponent_parameters = cal_exp$selected,
    n_tokens = n_tokens,
    maxent_vocabulary_size = vocabulary_size,
    markov_vocabulary_size = vocabulary_size,
    markov_persistence = markov_persistence
  )
  strict_parameters$maxent$block_size <- block_size
  strict_parameters$markov$block_size <- block_size
  strict_parameters$ssr$vocabulary_size <- vocabulary_size
  strict_parameters$ssr$initial_distribution <- "stationary"

  mimicking_parameters <- mimicking_benchmark_parameters(
    exponent_parameters = cal_exp$selected,
    n_tokens = n_tokens
  )

  cat("\n=== 2. Primary benchmark: exact marginal controls + latent mixture ===\n")
  strict <- run_strict_replications(
    n_rep = n_rep,
    first_seed = 2001L,
    n_tokens = n_tokens,
    parameters = strict_parameters,
    n_surrogates = n_surrogates,
    heaps_window = heaps_window,
    zipf_r_max = zipf_r_max,
    save_path = file.path(results_dir, "benchmark_strict_replications.csv"),
    curve_save_path = file.path(results_dir, "heaps_local_strict_replications.csv"),
    alpha_sensitivity_save_path = file.path(results_dir, "alpha_sensitivity_strict_replications.csv"),
    nmi_sensitivity_save_path = file.path(results_dir, "nmi_topk_sensitivity_strict.csv"),
    nmi_top_k_values = nmi_top_k_values,
    primary_nmi_top_k = primary_nmi_top_k
  )
  saveRDS(strict, file.path(results_dir, "benchmark_strict_replications.rds"))

  alpha_status <- .strict_alpha_status(strict)
  utils::write.csv(alpha_status, file.path(results_dir, "strict_alpha_status.csv"), row.names = FALSE)
  print(alpha_status)

  cat("\n=== 3. Bootstrap AUC / signed rank effects ===\n")
  auc_boot <- pairwise_auc_bootstrap(
    strict,
    n_boot = n_boot,
    conf_level = 0.95,
    seed = 20260808L
  )
  utils::write.csv(
    auc_boot,
    file.path(results_dir, "diagnostic_pairwise_auc_bootstrap_strict.csv"),
    row.names = FALSE
  )

  cat("\n=== 4. Compatible / mimicking benchmark ===\n")
  mimic <- run_mimicking_replications(
    n_rep = n_rep,
    first_seed = 2001L,
    n_tokens = n_tokens,
    parameters = mimicking_parameters,
    n_surrogates = n_surrogates,
    heaps_window = heaps_window,
    zipf_r_max = zipf_r_max,
    save_path = file.path(results_dir, "benchmark_mimicking_replications.csv"),
    curve_save_path = file.path(results_dir, "heaps_local_mimicking_replications.csv"),
    alpha_sensitivity_save_path = file.path(results_dir, "alpha_sensitivity_mimicking_replications.csv")
  )
  saveRDS(mimic, file.path(results_dir, "benchmark_mimicking_replications.rds"))

  cat("\n=== 5. Lower-rank cutoff sensitivity ===\n")
  rmin_sens <- .run_rmin_sensitivity(
    n_rep = n_rep,
    first_seed = 2001L,
    n_tokens = n_tokens,
    parameters = strict_parameters,
    r_min_values = r_min_values,
    r_max = zipf_r_max,
    save_path = file.path(results_dir, "alpha_rmin_sensitivity_strict.csv")
  )

  cat("\n=== 6. Markov persistence / effective-information check ===\n")
  rho_sens <- .run_markov_rho_sensitivity(
    n_rep = n_rep,
    first_seed = 12001L,
    n_tokens = n_tokens,
    vocabulary_size = vocabulary_size,
    block_size = block_size,
    rho_values = rho_values,
    r_min = 10L,
    r_max = zipf_r_max,
    save_path = file.path(results_dir, "markov_rho_sensitivity.csv")
  )

  cat("\n=== 7. SSR sequence-boundary sensitivity ===\n")
  ssr_boundary <- .run_ssr_boundary_sensitivity(
    n_rep = n_rep,
    first_seed = 22001L,
    n_tokens = n_tokens,
    vocabulary_size = vocabulary_size,
    block_size = block_size,
    n_surrogates = n_surrogates,
    top_k = primary_nmi_top_k,
    save_path = file.path(results_dir, "ssr_boundary_sensitivity.csv")
  )

  cat("\nDone. Render the QMD again without re-running to rebuild only tables/figures.\n")
  invisible(list(
    calibration = cal_exp,
    strict = strict,
    auc_bootstrap = auc_boot,
    mimicking = mimic,
    rmin_sensitivity = rmin_sens,
    rho_sensitivity = rho_sens,
    ssr_boundary_sensitivity = ssr_boundary
  ))
}

run_smoke_test <- function() {
  run_all(
    n_rep = 8L,
    n_cal_seeds = 3L,
    n_latent_cal_seeds = 8L,
    n_surrogates = 3L,
    heaps_window = 7L,
    n_boot = 250L,
    nmi_top_k_values = c(25L, 50L, 100L),
    rho_values = c(0, 0.35, 0.65),
    r_min_values = c(1L, 10L, 50L)
  )
}

run_paper_analysis <- function() {
  run_all(
    n_rep = 100L,
    n_cal_seeds = 10L,
    n_latent_cal_seeds = 50L,
    n_surrogates = 20L,
    heaps_window = 7L,
    n_boot = 2000L,
    nmi_top_k_values = c(25L, 50L, 100L),
    rho_values = c(0, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80),
    r_min_values = c(1L, 5L, 10L, 20L, 50L, 100L)
  )
}
