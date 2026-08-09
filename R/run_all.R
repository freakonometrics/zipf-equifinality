# End-to-end pipeline for v0.5.1.
# Designed for RStudio / Quarto.  Nothing auto-runs when this file is sourced.
#
# v0.5.1 keeps the v0.5 benchmark fixed except for one pre-specified change:
# latent-mixture calibration is stabilized by selecting its final parameter on
# 50 calibration seeds.  Held-out evaluation remains on seeds 2001+.

.load_zipf_benchmark_code <- function() {
  needed <- c(
    "simulate_maxent_zipf", "simulate_markov_zipf", "zipf_summary",
    "calibrate_benchmark", "calibrate_latent_stable",
    "run_strict_replications", "pairwise_auc_bootstrap"
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

.paper_outputs <- function(results_dir = "results") {
  file.path(results_dir, c(
    "benchmark_strict_replications.rds",
    "benchmark_mimicking_replications.rds",
    "benchmark_joint_strict_replications.rds",
    "diagnostic_pairwise_auc_bootstrap_strict.csv",
    "diagnostic_pairwise_auc_bootstrap_joint.csv",
    "nmi_topk_sensitivity_strict.csv",
    "alpha_sensitivity_strict_replications.csv",
    "alpha_sensitivity_mimicking_replications.csv",
    "calibration_exponent_summary.csv",
    "calibration_joint_summary.csv",
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
    target_vocabulary = 5000L,
    markov_persistence = 0.35,
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

  stopifnot(n_rep >= 2L, n_cal_seeds >= 2L, n_latent_cal_seeds >= 5L)
  stopifnot(n_surrogates >= 1L, n_boot >= 100L)
  stopifnot(primary_nmi_top_k %in% nmi_top_k_values)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  core_cal_seeds <- 1001L + seq_len(n_cal_seeds) - 1L
  latent_cal_seeds <- 1001L + seq_len(n_latent_cal_seeds) - 1L
  latent_screening_seeds <- 1001L + seq_len(min(10L, n_latent_cal_seeds)) - 1L

  cat("\n=== 1. Exponent-only calibration ===\n")
  cat("Core models use ", length(core_cal_seeds),
      " calibration seeds; latent mixture uses ", length(latent_cal_seeds),
      " for final selection.\n", sep = "")

  cal_exp_core <- calibrate_benchmark(
    mode = "exponent",
    seeds = core_cal_seeds,
    n_tokens = n_tokens,
    alpha_tolerance = 0.01,
    refine_exponent = TRUE,
    zipf_r_max = zipf_r_max,
    models = c("random", "simon", "ssr")
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
  print(write_calibration_outputs(cal_exp, results_dir))

  strict_parameters <- strict_benchmark_parameters(
    exponent_parameters = cal_exp$selected,
    n_tokens = n_tokens,
    maxent_vocabulary_size = 5000L,
    markov_vocabulary_size = 5000L,
    markov_persistence = markov_persistence
  )
  mimicking_parameters <- mimicking_benchmark_parameters(
    exponent_parameters = cal_exp$selected,
    n_tokens = n_tokens
  )

  cat("\n=== 2. Strict benchmark + non-SSR Markov control ===\n")
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

  cat("\n=== 5. Joint calibration on (alpha, V_n) ===\n")
  cal_joint_core <- calibrate_benchmark(
    mode = "joint",
    seeds = core_cal_seeds,
    n_tokens = n_tokens,
    target_alpha = 1,
    target_vocabulary = target_vocabulary,
    alpha_tolerance = 0.01,
    refine_joint = TRUE,
    zipf_r_max = zipf_r_max,
    models = c("maxent", "markov", "ssr")
  )
  cal_joint_core$seed_count_by_model <- stats::setNames(
    rep(length(core_cal_seeds), length(cal_joint_core$tables)),
    names(cal_joint_core$tables)
  )

  cal_joint_latent <- calibrate_latent_stable(
    mode = "joint",
    screening_seeds = latent_screening_seeds,
    final_seeds = latent_cal_seeds,
    n_tokens = n_tokens,
    target_alpha = 1,
    target_vocabulary = target_vocabulary,
    alpha_tolerance = 0.005,
    vocabulary_tolerance = 0.03,
    zipf_r_max = zipf_r_max
  )

  cal_joint <- combine_calibration_objects(
    cal_joint_core, cal_joint_latent,
    save_path = file.path(results_dir, "calibration_joint.rds")
  )
  print(write_calibration_outputs(cal_joint, results_dir))

  joint_parameters <- strict_joint_parameters(
    joint_parameters = cal_joint$selected,
    n_tokens = n_tokens
  )
  joint <- run_strict_replications(
    n_rep = n_rep,
    first_seed = 3001L,
    n_tokens = n_tokens,
    parameters = joint_parameters,
    n_surrogates = n_surrogates,
    heaps_window = heaps_window,
    zipf_r_max = zipf_r_max,
    save_path = file.path(results_dir, "benchmark_joint_strict_replications.csv"),
    primary_nmi_top_k = primary_nmi_top_k
  )
  joint$benchmark <- "joint matched alpha + V_n"
  utils::write.csv(joint, file.path(results_dir, "benchmark_joint_strict_replications.csv"), row.names = FALSE)
  saveRDS(joint, file.path(results_dir, "benchmark_joint_strict_replications.rds"))

  joint_auc_boot <- pairwise_auc_bootstrap(
    joint,
    diagnostics = c(
      "alpha", "alpha_mle", "vocabulary", "heaps_beta", "heaps_beta_change",
      "transition_toward_more_frequent", "transition_nmi_excess"
    ),
    n_boot = n_boot,
    conf_level = 0.95,
    seed = 20260809L
  )
  utils::write.csv(
    joint_auc_boot,
    file.path(results_dir, "diagnostic_pairwise_auc_bootstrap_joint.csv"),
    row.names = FALSE
  )

  cat("\nDone. Render the QMD again without re-running to rebuild only tables/figures.\n")
  invisible(list(
    exponent_calibration = cal_exp,
    strict = strict,
    auc_bootstrap = auc_boot,
    mimicking = mimic,
    joint_calibration = cal_joint,
    joint = joint,
    joint_auc_bootstrap = joint_auc_boot
  ))
}

run_smoke_test <- function() {
  run_all(
    n_rep = 10L,
    n_cal_seeds = 3L,
    n_latent_cal_seeds = 8L,
    n_surrogates = 3L,
    heaps_window = 7L,
    n_boot = 250L,
    nmi_top_k_values = c(25L, 50L, 100L)
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
    nmi_top_k_values = c(25L, 50L, 100L)
  )
}
