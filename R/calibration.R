# Calibration and replication protocol, v0.5.
#
# v0.5 adds a non-SSR Markov control with an exact finite-Zipf stationary
# marginal, and a strict joint-matching experiment that targets both the OLS
# Zipf exponent and observed final vocabulary size.
#
# Design A ("exponent") aggressively targets the marginal Zipf exponent only.
# It uses a broad mechanism-specific grid followed by a local refinement around
# the best coarse candidate. Vocabulary size and all structural diagnostics are
# held out.
#
# Design B ("joint") targets alpha and final vocabulary size jointly.
#
# The goal is not to force every mechanism to pass every design. Failure to
# enter the target tolerance under an enriched observation map is itself an
# identifying restriction.

benchmark_fallbacks <- function(mode = c("exponent", "joint"), n_tokens = 100000L) {
  mode <- match.arg(mode)
  if (mode == "exponent") {
    list(
      maxent = list(
        n_tokens = n_tokens, vocabulary_size = 5000L, exponent = 1, block_size = 250L
      ),
      markov = list(
        n_tokens = n_tokens, vocabulary_size = 5000L, exponent = 1,
        persistence = 0.35, block_size = 250L
      ),
      random = list(
        n_tokens = n_tokens, alphabet_size = 9L, p_delim = 0.56,
        letter_decay = 0.03, block_size = 250L
      ),
      simon = list(n_tokens = n_tokens, innovation_prob = 0.015),
      latent = list(
        n_sequences = 400L, sequence_length = as.integer(n_tokens / 400L),
        vocabulary_size = 5000L, scale_min = 1, scale_max = 8000
      ),
      ssr = list(n_tokens = n_tokens, vocabulary_size = 20000L)
    )
  } else {
    list(
      maxent = list(
        n_tokens = n_tokens, vocabulary_size = 5200L, exponent = 1, block_size = 250L
      ),
      markov = list(
        n_tokens = n_tokens, vocabulary_size = 5550L, exponent = 1,
        persistence = 0.35, block_size = 250L
      ),
      random = list(
        n_tokens = n_tokens, alphabet_size = 9L, p_delim = 0.62,
        letter_decay = 0.03, block_size = 250L
      ),
      simon = list(n_tokens = n_tokens, innovation_prob = 0.05),
      latent = list(
        n_sequences = 400L, sequence_length = as.integer(n_tokens / 400L),
        vocabulary_size = 5000L, scale_min = 1, scale_max = 8000
      ),
      ssr = list(n_tokens = n_tokens, vocabulary_size = 5000L)
    )
  }
}

.log_integer_grid <- function(from, to, length.out) {
  unique(as.integer(round(exp(seq(log(from), log(to), length.out = length.out)))))
}

candidate_grids <- function(mode = c("exponent", "joint")) {
  mode <- match.arg(mode)
  if (mode == "exponent") {
    list(
      random = expand.grid(
        alphabet_size = c(8L, 9L, 10L),
        p_delim = seq(0.46, 0.66, by = 0.04),
        letter_decay = c(0, 0.03, 0.06),
        KEEP.OUT.ATTRS = FALSE
      ),
      simon = data.frame(
        innovation_prob = unique(c(
          seq(0.002, 0.02, by = 0.003),
          seq(0.025, 0.10, by = 0.01)
        ))
      ),
      latent = expand.grid(
        vocabulary_size = c(5000L, 10000L),
        scale_min = c(0.5, 1, 2),
        scale_max = c(3000, 5000, 8000, 12000, 20000, 40000),
        KEEP.OUT.ATTRS = FALSE
      ),
      ssr = data.frame(
        vocabulary_size = .log_integer_grid(3000L, 200000L, 18L)
      )
    )
  } else {
    list(
      maxent = data.frame(
        vocabulary_size = seq(4800L, 6200L, by = 100L)
      ),
      markov = expand.grid(
        vocabulary_size = seq(5000L, 7500L, by = 100L),
        persistence = 0.35,
        KEEP.OUT.ATTRS = FALSE
      ),
      random = expand.grid(
        alphabet_size = c(8L, 9L, 10L),
        p_delim = c(0.56, 0.60, 0.62, 0.64, 0.66),
        letter_decay = c(0, 0.03, 0.06),
        KEEP.OUT.ATTRS = FALSE
      ),
      simon = data.frame(
        innovation_prob = seq(0.03, 0.07, by = 0.005)
      ),
      latent = expand.grid(
        vocabulary_size = c(4500L, 5000L, 5500L),
        scale_min = c(0.5, 1, 2),
        scale_max = c(5000, 8000, 12000, 20000),
        KEEP.OUT.ATTRS = FALSE
      ),
      ssr = data.frame(
        vocabulary_size = seq(4000L, 6500L, by = 250L)
      )
    )
  }
}

refinement_grid <- function(model, best) {
  if (model == "random") {
    p0 <- as.numeric(best$p_delim)
    data.frame(
      alphabet_size = as.integer(best$alphabet_size),
      p_delim = pmin(0.90, pmax(0.10, seq(p0 - 0.035, p0 + 0.035, by = 0.005))),
      letter_decay = as.numeric(best$letter_decay)
    )
  } else if (model == "simon") {
    p0 <- as.numeric(best$innovation_prob)
    lo <- max(0.0005, p0 - 0.02)
    hi <- min(0.20, p0 + 0.02)
    data.frame(innovation_prob = unique(seq(lo, hi, length.out = 21L)))
  } else if (model == "latent") {
    s0 <- as.numeric(best$scale_max)
    mult <- exp(seq(log(0.60), log(1.70), length.out = 15L))
    data.frame(
      vocabulary_size = as.integer(best$vocabulary_size),
      scale_min = as.numeric(best$scale_min),
      scale_max = unique(signif(s0 * mult, 8))
    )
  } else if (model == "ssr") {
    v0 <- as.numeric(best$vocabulary_size)
    mult <- exp(seq(log(0.45), log(2.25), length.out = 17L))
    data.frame(vocabulary_size = unique(pmax(100L, as.integer(round(v0 * mult)))))
  } else {
    stop("Unknown model for refinement: ", model)
  }
}

joint_refinement_grid <- function(model, best) {
  if (model == "maxent") {
    v0 <- as.integer(best$vocabulary_size)
    data.frame(vocabulary_size = seq(max(100L, v0 - 250L), v0 + 250L, by = 25L))
  } else if (model == "markov") {
    v0 <- as.integer(best$vocabulary_size)
    data.frame(
      vocabulary_size = seq(max(100L, v0 - 400L), v0 + 400L, by = 25L),
      persistence = as.numeric(best$persistence)
    )
  } else if (model == "latent") {
    v0 <- as.integer(best$vocabulary_size)
    s0 <- as.numeric(best$scale_max)
    expand.grid(
      vocabulary_size = unique(pmax(500L, v0 + c(-250L, -100L, 0L, 100L, 250L))),
      scale_min = as.numeric(best$scale_min),
      scale_max = unique(signif(s0 * exp(seq(log(0.75), log(1.35), length.out = 9L)), 8)),
      KEEP.OUT.ATTRS = FALSE
    )
  } else if (model == "ssr") {
    v0 <- as.integer(best$vocabulary_size)
    data.frame(vocabulary_size = seq(max(100L, v0 - 750L), v0 + 750L, by = 50L))
  } else {
    stop("No joint refinement grid for model: ", model)
  }
}

calibration_loss <- function(
    alpha,
    vocabulary,
    mode = c("exponent", "joint"),
    target_alpha = 1,
    target_vocabulary = 5000,
    alpha_scale = 0.025,
    log_vocabulary_scale = 0.03) {
  mode <- match.arg(mode)
  alpha_term <- ((alpha - target_alpha) / alpha_scale)^2
  if (mode == "exponent") return(alpha_term)
  vocab_term <- ((log(vocabulary) - log(target_vocabulary)) /
                   log_vocabulary_scale)^2
  alpha_term + vocab_term
}

.simulate_candidate <- function(model, pars, seed, n_tokens) {
  if (model == "maxent") {
    simulate_maxent_zipf(
      n_tokens = n_tokens,
      vocabulary_size = as.integer(pars$vocabulary_size),
      exponent = 1, block_size = 250L, seed = seed
    )
  } else if (model == "markov") {
    simulate_markov_zipf(
      n_tokens = n_tokens,
      vocabulary_size = as.integer(pars$vocabulary_size),
      exponent = 1,
      persistence = if (!is.null(pars$persistence)) as.numeric(pars$persistence) else 0.35,
      block_size = 250L, seed = seed
    )
  } else if (model == "random") {
    simulate_random_segmentation(
      n_tokens = n_tokens,
      alphabet_size = as.integer(pars$alphabet_size),
      p_delim = as.numeric(pars$p_delim),
      letter_decay = as.numeric(pars$letter_decay),
      block_size = 250L,
      seed = seed
    )
  } else if (model == "simon") {
    simulate_simon(
      n_tokens = n_tokens,
      innovation_prob = as.numeric(pars$innovation_prob),
      seed = seed
    )
  } else if (model == "latent") {
    simulate_latent_mixture(
      n_sequences = 400L,
      sequence_length = as.integer(n_tokens / 400L),
      vocabulary_size = if (!is.null(pars$vocabulary_size)) as.integer(pars$vocabulary_size) else 5000L,
      scale_min = if (!is.null(pars$scale_min)) as.numeric(pars$scale_min) else 1,
      scale_max = as.numeric(pars$scale_max),
      seed = seed
    )
  } else if (model == "ssr") {
    simulate_ssr(
      n_tokens = n_tokens,
      vocabulary_size = as.integer(pars$vocabulary_size),
      seed = seed
    )
  } else {
    stop("Unknown model: ", model)
  }
}

.evaluate_grid <- function(
    model,
    grid,
    seeds,
    n_tokens,
    mode,
    target_alpha,
    target_vocabulary,
    zipf_r_max = 500L,
    stage = "coarse") {
  rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    pars <- as.list(grid[i, , drop = FALSE])
    vals <- lapply(seeds, function(s) {
      z <- zipf_summary(.simulate_candidate(model, pars, s, n_tokens), r_max = zipf_r_max)
      z[, c("alpha", "vocabulary"), drop = FALSE]
    })
    vals <- do.call(rbind, vals)

    mean_alpha <- mean(vals$alpha, na.rm = TRUE)
    sd_alpha <- stats::sd(vals$alpha, na.rm = TRUE)
    mean_vocabulary <- mean(vals$vocabulary, na.rm = TRUE)
    sd_vocabulary <- stats::sd(vals$vocabulary, na.rm = TRUE)
    alpha_rmse <- sqrt(mean((vals$alpha - target_alpha)^2, na.rm = TRUE))

    rows[[i]] <- cbind(
      grid[i, , drop = FALSE],
      data.frame(
        stage = stage,
        mean_alpha = mean_alpha,
        sd_alpha = sd_alpha,
        abs_alpha_bias = abs(mean_alpha - target_alpha),
        alpha_rmse = alpha_rmse,
        mean_vocabulary = mean_vocabulary,
        sd_vocabulary = sd_vocabulary,
        loss = calibration_loss(
          mean_alpha, mean_vocabulary,
          mode = mode,
          target_alpha = target_alpha,
          target_vocabulary = target_vocabulary
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  out <- do.call(rbind, rows)
  if (mode == "exponent") {
    out <- out[order(out$abs_alpha_bias, out$alpha_rmse, out$sd_alpha), , drop = FALSE]
  } else {
    out <- out[order(out$loss), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.parameter_names <- function(model) {
  switch(
    model,
    maxent = "vocabulary_size",
    markov = c("vocabulary_size", "persistence"),
    random = c("alphabet_size", "p_delim", "letter_decay"),
    simon = "innovation_prob",
    latent = c("vocabulary_size", "scale_min", "scale_max"),
    ssr = "vocabulary_size",
    stop("Unknown model: ", model)
  )
}

.selected_parameter_list <- function(model, row, n_tokens) {
  if (model == "maxent") {
    list(
      n_tokens = n_tokens, vocabulary_size = as.integer(row$vocabulary_size),
      exponent = 1, block_size = 250L
    )
  } else if (model == "markov") {
    list(
      n_tokens = n_tokens, vocabulary_size = as.integer(row$vocabulary_size),
      exponent = 1,
      persistence = if (!is.null(row$persistence)) as.numeric(row$persistence) else 0.35,
      block_size = 250L
    )
  } else if (model == "random") {
    list(
      n_tokens = n_tokens,
      alphabet_size = as.integer(row$alphabet_size),
      p_delim = as.numeric(row$p_delim),
      letter_decay = as.numeric(row$letter_decay),
      block_size = 250L
    )
  } else if (model == "simon") {
    list(
      n_tokens = n_tokens,
      innovation_prob = as.numeric(row$innovation_prob)
    )
  } else if (model == "latent") {
    list(
      n_sequences = 400L,
      sequence_length = as.integer(n_tokens / 400L),
      vocabulary_size = as.integer(row$vocabulary_size),
      scale_min = as.numeric(row$scale_min),
      scale_max = as.numeric(row$scale_max)
    )
  } else if (model == "ssr") {
    list(
      n_tokens = n_tokens,
      vocabulary_size = as.integer(row$vocabulary_size)
    )
  } else {
    stop("Unknown model: ", model)
  }
}

calibrate_benchmark <- function(
    mode = c("exponent", "joint"),
    seeds = 1001:1010,
    n_tokens = 100000L,
    target_alpha = 1,
    target_vocabulary = 5000L,
    alpha_tolerance = 0.01,
    vocabulary_tolerance = 0.03,
    refine_exponent = TRUE,
    refine_joint = TRUE,
    zipf_r_max = 500L,
    models = NULL,
    save_path = NULL) {
  mode <- match.arg(mode)
  grids <- candidate_grids(mode)
  if (!is.null(models)) {
    missing_models <- setdiff(models, names(grids))
    if (length(missing_models)) {
      stop("Unknown calibration model(s): ", paste(missing_models, collapse = ", "))
    }
    grids <- grids[models]
  }
  tables <- list()
  selected <- list()

  for (model in names(grids)) {
    cat("  calibrating", model, "[", mode, "]\n")
    coarse <- .evaluate_grid(
      model = model,
      grid = grids[[model]],
      seeds = seeds,
      n_tokens = n_tokens,
      mode = mode,
      target_alpha = target_alpha,
      target_vocabulary = target_vocabulary,
      zipf_r_max = zipf_r_max,
      stage = "coarse"
    )

    tab <- coarse
    if (mode == "exponent" && isTRUE(refine_exponent)) {
      pnames <- .parameter_names(model)
      best_coarse <- as.list(coarse[1L, pnames, drop = FALSE])
      fine_grid <- refinement_grid(model, best_coarse)
      fine <- .evaluate_grid(
        model = model,
        grid = fine_grid,
        seeds = seeds,
        n_tokens = n_tokens,
        mode = mode,
        target_alpha = target_alpha,
        target_vocabulary = target_vocabulary,
        zipf_r_max = zipf_r_max,
        stage = "refined"
      )
      tab <- rbind(coarse, fine)
      if (mode == "exponent") {
        tab <- tab[order(tab$abs_alpha_bias, tab$alpha_rmse, tab$sd_alpha), , drop = FALSE]
      }
      rownames(tab) <- NULL
    }
    if (mode == "joint" && isTRUE(refine_joint) && model %in% c("maxent", "markov", "latent", "ssr")) {
      pnames <- .parameter_names(model)
      best_coarse <- as.list(coarse[1L, pnames, drop = FALSE])
      fine_grid <- joint_refinement_grid(model, best_coarse)
      fine <- .evaluate_grid(
        model = model, grid = fine_grid, seeds = seeds, n_tokens = n_tokens,
        mode = mode, target_alpha = target_alpha,
        target_vocabulary = target_vocabulary, zipf_r_max = zipf_r_max,
        stage = "refined"
      )
      tab <- rbind(coarse, fine)
      tab <- tab[order(tab$loss), , drop = FALSE]
      rownames(tab) <- NULL
    }

    tables[[model]] <- tab
    selected[[model]] <- .selected_parameter_list(model, tab[1L, , drop = FALSE], n_tokens)
  }

  status <- do.call(rbind, lapply(names(tables), function(model) {
    x <- tables[[model]][1L, , drop = FALSE]
    data.frame(
      model = model,
      mean_alpha = x$mean_alpha,
      sd_alpha = x$sd_alpha,
      abs_alpha_bias = x$abs_alpha_bias,
      mean_vocabulary = x$mean_vocabulary,
      relative_vocabulary_error = abs(x$mean_vocabulary - target_vocabulary) / target_vocabulary,
      within_alpha_tolerance = abs(x$mean_alpha - target_alpha) <= alpha_tolerance,
      within_vocabulary_tolerance =
        abs(x$mean_vocabulary - target_vocabulary) / target_vocabulary <= vocabulary_tolerance,
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    version = "0.5",
    mode = mode,
    target_alpha = target_alpha,
    target_vocabulary = target_vocabulary,
    alpha_tolerance = alpha_tolerance,
    vocabulary_tolerance = vocabulary_tolerance,
    zipf_r_max = zipf_r_max,
    seeds = seeds,
    n_tokens = n_tokens,
    selected = selected,
    tables = tables,
    status = status
  )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, save_path)
  }
  out
}

calibration_status <- function(calibration) {
  if (!is.null(calibration$status)) return(calibration$status)
  do.call(rbind, lapply(names(calibration$tables), function(model) {
    x <- calibration$tables[[model]][1L, , drop = FALSE]
    data.frame(
      model = model,
      mean_alpha = x$mean_alpha,
      sd_alpha = x$sd_alpha,
      abs_alpha_bias = abs(x$mean_alpha - calibration$target_alpha),
      mean_vocabulary = x$mean_vocabulary,
      within_alpha_tolerance = abs(x$mean_alpha - calibration$target_alpha) <= 0.01,
      stringsAsFactors = FALSE
    )
  }))
}

load_calibration_or_fallback <- function(
    mode = c("exponent", "joint"),
    n_tokens = 100000L,
    results_dir = "results") {
  mode <- match.arg(mode)
  path <- file.path(results_dir, paste0("calibration_", mode, ".rds"))
  if (file.exists(path)) {
    cal <- readRDS(path)
    return(cal$selected)
  }
  warning(
    "No saved ", mode, " calibration found at ", path,
    ". Using transparent fallback parameters. Run R/run_calibration.R first for the paper analysis."
  )
  benchmark_fallbacks(mode, n_tokens = n_tokens)
}

simulate_benchmark_once <- function(
    seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL,
    mode = c("exponent", "joint")) {
  mode <- match.arg(mode)
  if (is.null(parameters)) parameters <- benchmark_fallbacks(mode, n_tokens)

  list(
    random = do.call(simulate_random_segmentation,
                     c(parameters$random, list(seed = seed + 11L))),
    simon = do.call(simulate_simon,
                    c(parameters$simon, list(seed = seed + 23L))),
    latent = do.call(simulate_latent_mixture,
                     c(parameters$latent, list(seed = seed + 37L))),
    ssr = do.call(simulate_ssr,
                  c(parameters$ssr, list(seed = seed + 53L)))
  )
}

run_replications <- function(
    n_rep = 100L,
    first_seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL,
    mode = c("exponent", "joint"),
    n_surrogates = 5L,
    heaps_window = 7L,
    zipf_r_max = 500L,
    save_path = NULL,
    curve_save_path = NULL,
    alpha_sensitivity_save_path = NULL,
    alpha_r_max_values = c(300L, 500L, 1000L, 2000L, 3000L)) {
  mode <- match.arg(mode)
  if (is.null(parameters)) {
    parameters <- load_calibration_or_fallback(mode, n_tokens = n_tokens)
  }

  all_results <- vector("list", n_rep)
  all_curves <- if (!is.null(curve_save_path)) vector("list", n_rep) else NULL
  all_sensitivity <- if (!is.null(alpha_sensitivity_save_path)) vector("list", n_rep) else NULL

  for (b in seq_len(n_rep)) {
    seed <- first_seed + b - 1L
    sims <- simulate_benchmark_once(
      seed = seed,
      n_tokens = n_tokens,
      parameters = parameters,
      mode = mode
    )

    d <- do.call(rbind, lapply(seq_along(sims), function(j) {
      benchmark_diagnostics(
        sims[[j]],
        n_surrogates = n_surrogates,
        surrogate_seed = 100000L + 100L * seed + j,
        heaps_window = heaps_window,
        zipf_r_max = zipf_r_max
      )
    }))
    d$replicate <- b
    d$seed <- seed
    d$calibration_mode <- mode
    all_results[[b]] <- d

    if (!is.null(curve_save_path)) {
      curves_b <- do.call(rbind, lapply(sims, function(x) {
        h <- heaps_local_curve(x, window = heaps_window)
        h$generator <- unique(x$generator)[1L]
        h
      }))
      curves_b$replicate <- b
      curves_b$seed <- seed
      curves_b$calibration_mode <- mode
      all_curves[[b]] <- curves_b
    }

    if (!is.null(alpha_sensitivity_save_path)) {
      sens_b <- do.call(rbind, lapply(sims, function(x) {
        z <- zipf_sensitivity(x, r_max_values = alpha_r_max_values)
        z
      }))
      sens_b$replicate <- b
      sens_b$seed <- seed
      sens_b$calibration_mode <- mode
      all_sensitivity[[b]] <- sens_b
    }
  }

  out <- do.call(rbind, all_results)
  rownames(out) <- NULL

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, save_path, row.names = FALSE)
  }

  if (!is.null(curve_save_path)) {
    curves <- do.call(rbind, all_curves)
    rownames(curves) <- NULL
    dir.create(dirname(curve_save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(curves, curve_save_path, row.names = FALSE)
  }

  if (!is.null(alpha_sensitivity_save_path)) {
    sens <- do.call(rbind, all_sensitivity)
    rownames(sens) <- NULL
    dir.create(dirname(alpha_sensitivity_save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(sens, alpha_sensitivity_save_path, row.names = FALSE)
  }

  out
}

write_calibration_outputs <- function(calibration, results_dir = "results") {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  mode <- calibration$mode
  status <- calibration_status(calibration)

  summary_rows <- do.call(rbind, lapply(names(calibration$tables), function(model) {
    x <- calibration$tables[[model]][1L, , drop = FALSE]
    data.frame(
      model = model,
      stage = x$stage,
      mean_alpha = x$mean_alpha,
      sd_alpha = x$sd_alpha,
      abs_alpha_bias = x$abs_alpha_bias,
      alpha_rmse = x$alpha_rmse,
      mean_vocabulary = x$mean_vocabulary,
      sd_vocabulary = x$sd_vocabulary,
      loss = x$loss,
      stringsAsFactors = FALSE
    )
  }))

  status_cols <- intersect(
    c("model", "relative_vocabulary_error", "within_alpha_tolerance",
      "within_vocabulary_tolerance"),
    names(status)
  )
  summary_rows <- merge(summary_rows, status[, status_cols, drop = FALSE],
                        by = "model", all.x = TRUE, sort = FALSE)

  utils::write.csv(
    summary_rows,
    file.path(results_dir, paste0("calibration_", mode, "_summary.csv")),
    row.names = FALSE
  )

  for (model in names(calibration$tables)) {
    utils::write.csv(
      calibration$tables[[model]],
      file.path(results_dir, paste0("calibration_", mode, "_", model, "_grid.csv")),
      row.names = FALSE
    )
  }

  invisible(summary_rows)
}

# -------------------------------------------------------------------------
# v0.4 benchmark split
# -------------------------------------------------------------------------
# The strict benchmark contains mechanisms that are genuinely close under the
# prespecified Zipf statistic: an exact finite MaxEnt/Zipf distribution, the
# latent-mixture generator, and canonical sample-space reduction. Random
# segmentation and Simon reinforcement are retained as Zipf-mimicking / 
# Zipf-compatible mechanisms and studied separately as an observation-design
# sensitivity benchmark.

strict_benchmark_parameters <- function(
    exponent_parameters = NULL,
    n_tokens = 100000L,
    maxent_vocabulary_size = 5000L,
    markov_vocabulary_size = 5000L,
    markov_persistence = 0.35) {
  if (is.null(exponent_parameters)) {
    exponent_parameters <- load_calibration_or_fallback(
      "exponent", n_tokens = n_tokens
    )
  }
  list(
    maxent = list(
      n_tokens = n_tokens,
      vocabulary_size = as.integer(maxent_vocabulary_size),
      exponent = 1,
      block_size = 250L
    ),
    markov = list(
      n_tokens = n_tokens,
      vocabulary_size = as.integer(markov_vocabulary_size),
      exponent = 1,
      persistence = as.numeric(markov_persistence),
      block_size = 250L
    ),
    latent = exponent_parameters$latent,
    ssr = exponent_parameters$ssr
  )
}

strict_joint_parameters <- function(
    joint_parameters = NULL,
    n_tokens = 100000L) {
  if (is.null(joint_parameters)) {
    joint_parameters <- load_calibration_or_fallback("joint", n_tokens = n_tokens)
  }
  needed <- c("maxent", "markov", "latent", "ssr")
  missing <- setdiff(needed, names(joint_parameters))
  if (length(missing)) {
    stop("Joint calibration is missing: ", paste(missing, collapse = ", "))
  }
  joint_parameters[needed]
}

mimicking_benchmark_parameters <- function(
    exponent_parameters = NULL,
    n_tokens = 100000L) {
  if (is.null(exponent_parameters)) {
    exponent_parameters <- load_calibration_or_fallback(
      "exponent", n_tokens = n_tokens
    )
  }
  list(
    random = exponent_parameters$random,
    simon = exponent_parameters$simon
  )
}

simulate_strict_benchmark_once <- function(
    seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL) {
  if (is.null(parameters)) {
    parameters <- strict_benchmark_parameters(n_tokens = n_tokens)
  }
  list(
    maxent = do.call(
      simulate_maxent_zipf,
      c(parameters$maxent, list(seed = seed + 5L))
    ),
    markov = do.call(
      simulate_markov_zipf,
      c(parameters$markov, list(seed = seed + 19L))
    ),
    latent = do.call(
      simulate_latent_mixture,
      c(parameters$latent, list(seed = seed + 37L))
    ),
    ssr = do.call(
      simulate_ssr,
      c(parameters$ssr, list(seed = seed + 53L))
    )
  )
}

simulate_mimicking_benchmark_once <- function(
    seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL) {
  if (is.null(parameters)) {
    parameters <- mimicking_benchmark_parameters(n_tokens = n_tokens)
  }
  list(
    random = do.call(
      simulate_random_segmentation,
      c(parameters$random, list(seed = seed + 11L))
    ),
    simon = do.call(
      simulate_simon,
      c(parameters$simon, list(seed = seed + 23L))
    )
  )
}

.run_split_replications <- function(
    simulator,
    benchmark_label,
    n_rep = 100L,
    first_seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL,
    n_surrogates = 5L,
    heaps_window = 7L,
    zipf_r_max = 500L,
    save_path = NULL,
    curve_save_path = NULL,
    alpha_sensitivity_save_path = NULL,
    nmi_sensitivity_save_path = NULL,
    nmi_top_k_values = c(25L, 50L, 100L),
    primary_nmi_top_k = 50L,
    alpha_r_max_values = c(300L, 500L, 1000L, 2000L, 3000L)) {
  all_results <- vector("list", n_rep)
  all_curves <- if (!is.null(curve_save_path)) vector("list", n_rep) else NULL
  all_sensitivity <- if (!is.null(alpha_sensitivity_save_path)) vector("list", n_rep) else NULL
  all_nmi <- if (!is.null(nmi_sensitivity_save_path)) vector("list", n_rep) else NULL

  for (b in seq_len(n_rep)) {
    seed <- first_seed + b - 1L
    sims <- simulator(
      seed = seed,
      n_tokens = n_tokens,
      parameters = parameters
    )

    nmi_rows_b <- vector("list", length(sims))
    d_rows <- vector("list", length(sims))
    for (j in seq_along(sims)) {
      ti_primary <- NULL
      if (!is.null(nmi_sensitivity_save_path)) {
        ti_multi <- transition_information_excess_multi(
          sims[[j]],
          top_k_values = nmi_top_k_values,
          n_surrogates = n_surrogates,
          seed = 100000L + 100L * seed + j
        )
        ti_multi$generator <- unique(sims[[j]]$generator)[1L]
        nmi_rows_b[[j]] <- ti_multi
        if (!primary_nmi_top_k %in% ti_multi$top_k) {
          stop("primary_nmi_top_k must be included in nmi_top_k_values")
        }
        ti_primary <- ti_multi[ti_multi$top_k == primary_nmi_top_k, , drop = FALSE]
      }
      d_rows[[j]] <- benchmark_diagnostics(
        sims[[j]],
        n_surrogates = n_surrogates,
        surrogate_seed = 100000L + 100L * seed + j,
        heaps_window = heaps_window,
        zipf_r_max = zipf_r_max,
        nmi_top_k = primary_nmi_top_k,
        transition_info = ti_primary
      )
    }
    d <- do.call(rbind, d_rows)
    d$replicate <- b
    d$seed <- seed
    d$benchmark <- benchmark_label
    all_results[[b]] <- d

    if (!is.null(curve_save_path)) {
      curves_b <- do.call(rbind, lapply(sims, function(x) {
        h <- heaps_local_curve(x, window = heaps_window)
        h$generator <- unique(x$generator)[1L]
        h
      }))
      curves_b$replicate <- b
      curves_b$seed <- seed
      curves_b$benchmark <- benchmark_label
      all_curves[[b]] <- curves_b
    }

    if (!is.null(alpha_sensitivity_save_path)) {
      sens_b <- do.call(rbind, lapply(sims, function(x) {
        zipf_sensitivity(x, r_max_values = alpha_r_max_values)
      }))
      sens_b$replicate <- b
      sens_b$seed <- seed
      sens_b$benchmark <- benchmark_label
      all_sensitivity[[b]] <- sens_b
    }

    if (!is.null(nmi_sensitivity_save_path)) {
      nmi_b <- do.call(rbind, nmi_rows_b)
      nmi_b$replicate <- b
      nmi_b$seed <- seed
      nmi_b$benchmark <- benchmark_label
      all_nmi[[b]] <- nmi_b
    }
  }

  out <- do.call(rbind, all_results)
  rownames(out) <- NULL

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, save_path, row.names = FALSE)
  }
  if (!is.null(curve_save_path)) {
    curves <- do.call(rbind, all_curves)
    rownames(curves) <- NULL
    dir.create(dirname(curve_save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(curves, curve_save_path, row.names = FALSE)
  }
  if (!is.null(alpha_sensitivity_save_path)) {
    sens <- do.call(rbind, all_sensitivity)
    rownames(sens) <- NULL
    dir.create(dirname(alpha_sensitivity_save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(sens, alpha_sensitivity_save_path, row.names = FALSE)
  }
  if (!is.null(nmi_sensitivity_save_path)) {
    nmi <- do.call(rbind, all_nmi)
    rownames(nmi) <- NULL
    dir.create(dirname(nmi_sensitivity_save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(nmi, nmi_sensitivity_save_path, row.names = FALSE)
  }
  out
}

run_strict_replications <- function(
    n_rep = 100L,
    first_seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL,
    n_surrogates = 5L,
    heaps_window = 7L,
    zipf_r_max = 500L,
    save_path = NULL,
    curve_save_path = NULL,
    alpha_sensitivity_save_path = NULL,
    nmi_sensitivity_save_path = NULL,
    nmi_top_k_values = c(25L, 50L, 100L),
    primary_nmi_top_k = 50L) {
  if (is.null(parameters)) {
    parameters <- strict_benchmark_parameters(n_tokens = n_tokens)
  }
  .run_split_replications(
    simulator = simulate_strict_benchmark_once,
    benchmark_label = "strict equifinality",
    n_rep = n_rep,
    first_seed = first_seed,
    n_tokens = n_tokens,
    parameters = parameters,
    n_surrogates = n_surrogates,
    heaps_window = heaps_window,
    zipf_r_max = zipf_r_max,
    save_path = save_path,
    curve_save_path = curve_save_path,
    alpha_sensitivity_save_path = alpha_sensitivity_save_path,
    nmi_sensitivity_save_path = nmi_sensitivity_save_path,
    nmi_top_k_values = nmi_top_k_values,
    primary_nmi_top_k = primary_nmi_top_k
  )
}

run_mimicking_replications <- function(
    n_rep = 100L,
    first_seed = 2001L,
    n_tokens = 100000L,
    parameters = NULL,
    n_surrogates = 5L,
    heaps_window = 7L,
    zipf_r_max = 500L,
    save_path = NULL,
    curve_save_path = NULL,
    alpha_sensitivity_save_path = NULL) {
  if (is.null(parameters)) {
    parameters <- mimicking_benchmark_parameters(n_tokens = n_tokens)
  }
  .run_split_replications(
    simulator = simulate_mimicking_benchmark_once,
    benchmark_label = "Zipf-compatible / mimicking",
    n_rep = n_rep,
    first_seed = first_seed,
    n_tokens = n_tokens,
    parameters = parameters,
    n_surrogates = n_surrogates,
    heaps_window = heaps_window,
    zipf_r_max = zipf_r_max,
    save_path = save_path,
    curve_save_path = curve_save_path,
    alpha_sensitivity_save_path = alpha_sensitivity_save_path
  )
}

# -------------------------------------------------------------------------
# v0.5.1: stable latent-mixture calibration
# -------------------------------------------------------------------------
# The latent-mixture generator has materially larger Monte Carlo variability
# than the other benchmark constructions.  To avoid selecting a parameter
# that looks close to alpha=1 on only a handful of calibration seeds, v0.5.1
# uses a two-stage design:
#   (i) screen the broad latent grid on a small subset of calibration seeds;
#  (ii) evaluate a restricted candidate/refinement set on a larger,
#       prespecified 50-seed calibration panel (1001:1050 by default).
# Held-out evaluation seeds (2001+) are never used for model selection.

calibrate_latent_stable <- function(
    mode = c("exponent", "joint"),
    screening_seeds = 1001:1010,
    final_seeds = 1001:1050,
    n_tokens = 100000L,
    target_alpha = 1,
    target_vocabulary = 5000L,
    alpha_tolerance = 0.005,
    vocabulary_tolerance = 0.03,
    zipf_r_max = 500L,
    n_screen = 3L,
    save_path = NULL) {
  mode <- match.arg(mode)
  n_screen <- max(1L, as.integer(n_screen))
  screening_seeds <- as.integer(screening_seeds)
  final_seeds <- as.integer(final_seeds)
  if (!length(final_seeds)) stop("final_seeds must be non-empty")

  coarse_grid <- candidate_grids(mode)$latent
  cat("  latent stable calibration [", mode, "]: screening ",
      nrow(coarse_grid), " candidates on ", length(screening_seeds),
      " seeds\n", sep = "")
  coarse <- .evaluate_grid(
    model = "latent",
    grid = coarse_grid,
    seeds = screening_seeds,
    n_tokens = n_tokens,
    mode = mode,
    target_alpha = target_alpha,
    target_vocabulary = target_vocabulary,
    zipf_r_max = zipf_r_max,
    stage = "screening"
  )

  n_keep <- min(n_screen, nrow(coarse))
  top <- coarse[seq_len(n_keep), , drop = FALSE]
  candidate_rows <- list()
  k <- 1L
  for (i in seq_len(nrow(top))) {
    pnames <- .parameter_names("latent")
    best <- as.list(top[i, pnames, drop = FALSE])
    candidate_rows[[k]] <- as.data.frame(best, stringsAsFactors = FALSE)
    k <- k + 1L
    fine <- if (mode == "exponent") {
      refinement_grid("latent", best)
    } else {
      joint_refinement_grid("latent", best)
    }
    candidate_rows[[k]] <- fine
    k <- k + 1L
  }
  final_grid <- do.call(rbind, candidate_rows)
  final_grid <- unique(final_grid[, .parameter_names("latent"), drop = FALSE])
  rownames(final_grid) <- NULL

  cat("  latent stable calibration [", mode, "]: final selection among ",
      nrow(final_grid), " candidates on ", length(final_seeds),
      " calibration seeds\n", sep = "")
  final <- .evaluate_grid(
    model = "latent",
    grid = final_grid,
    seeds = final_seeds,
    n_tokens = n_tokens,
    mode = mode,
    target_alpha = target_alpha,
    target_vocabulary = target_vocabulary,
    zipf_r_max = zipf_r_max,
    stage = paste0(length(final_seeds), "-seed-final")
  )

  # Put the final 50-seed ranking first so existing output/reporting helpers
  # always read the actually selected candidate from row 1.
  tab <- rbind(final, coarse)
  rownames(tab) <- NULL
  best_row <- final[1L, , drop = FALSE]
  selected <- list(latent = .selected_parameter_list("latent", best_row, n_tokens))

  status <- data.frame(
    model = "latent",
    mean_alpha = best_row$mean_alpha,
    sd_alpha = best_row$sd_alpha,
    abs_alpha_bias = best_row$abs_alpha_bias,
    mean_vocabulary = best_row$mean_vocabulary,
    relative_vocabulary_error = abs(best_row$mean_vocabulary - target_vocabulary) / target_vocabulary,
    within_alpha_tolerance = abs(best_row$mean_alpha - target_alpha) <= alpha_tolerance,
    within_vocabulary_tolerance =
      abs(best_row$mean_vocabulary - target_vocabulary) / target_vocabulary <= vocabulary_tolerance,
    stringsAsFactors = FALSE
  )

  out <- list(
    version = "0.5.1",
    mode = mode,
    target_alpha = target_alpha,
    target_vocabulary = target_vocabulary,
    alpha_tolerance = alpha_tolerance,
    vocabulary_tolerance = vocabulary_tolerance,
    zipf_r_max = zipf_r_max,
    seeds = final_seeds,
    screening_seeds = screening_seeds,
    n_tokens = n_tokens,
    selected = selected,
    tables = list(latent = tab),
    status = status,
    seed_count_by_model = c(latent = length(final_seeds))
  )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, save_path)
  }
  out
}

combine_calibration_objects <- function(core, latent, save_path = NULL) {
  if (!identical(core$mode, latent$mode)) stop("Calibration modes do not match")
  out <- core
  out$version <- "0.5.1"
  out$selected$latent <- latent$selected$latent
  out$tables$latent <- latent$tables$latent

  core_status <- core$status[core$status$model != "latent", , drop = FALSE]
  out$status <- rbind(core_status, latent$status)
  rownames(out$status) <- NULL

  core_counts <- if (!is.null(core$seed_count_by_model)) core$seed_count_by_model else {
    stats::setNames(rep(length(core$seeds), length(core$tables)), names(core$tables))
  }
  latent_counts <- if (!is.null(latent$seed_count_by_model)) latent$seed_count_by_model else c(latent = length(latent$seeds))
  out$seed_count_by_model <- c(core_counts[names(core$tables)], latent_counts["latent"])
  out$seed_design <- list(
    core = core$seeds,
    latent_screening = latent$screening_seeds,
    latent_final = latent$seeds
  )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, save_path)
  }
  out
}

# Override the v0.5 writer to record how many calibration seeds supported each
# selected parameter.  The selected row remains the first row of each table.
write_calibration_outputs <- function(calibration, results_dir = "results") {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  mode <- calibration$mode
  status <- calibration_status(calibration)

  summary_rows <- do.call(rbind, lapply(names(calibration$tables), function(model) {
    x <- calibration$tables[[model]][1L, , drop = FALSE]
    n_seed <- if (!is.null(calibration$seed_count_by_model) && model %in% names(calibration$seed_count_by_model)) {
      as.integer(calibration$seed_count_by_model[[model]])
    } else {
      length(calibration$seeds)
    }
    data.frame(
      model = model,
      stage = x$stage,
      n_calibration_seeds = n_seed,
      mean_alpha = x$mean_alpha,
      sd_alpha = x$sd_alpha,
      abs_alpha_bias = x$abs_alpha_bias,
      alpha_rmse = x$alpha_rmse,
      mean_vocabulary = x$mean_vocabulary,
      sd_vocabulary = x$sd_vocabulary,
      loss = x$loss,
      stringsAsFactors = FALSE
    )
  }))

  status_cols <- intersect(
    c("model", "relative_vocabulary_error", "within_alpha_tolerance",
      "within_vocabulary_tolerance"),
    names(status)
  )
  summary_rows <- merge(summary_rows, status[, status_cols, drop = FALSE],
                        by = "model", all.x = TRUE, sort = FALSE)

  utils::write.csv(
    summary_rows,
    file.path(results_dir, paste0("calibration_", mode, "_summary.csv")),
    row.names = FALSE
  )

  for (model in names(calibration$tables)) {
    utils::write.csv(
      calibration$tables[[model]],
      file.path(results_dir, paste0("calibration_", mode, "_", model, "_grid.csv")),
      row.names = FALSE
    )
  }
  invisible(summary_rows)
}
