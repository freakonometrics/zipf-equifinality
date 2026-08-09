# Diagnostics for separating Zipf-equivalent generators.
# v0.6 retains:
#   - a restricted-domain discrete-likelihood estimate of the Zipf exponent on the
#     same prespecified rank window used by the descriptive OLS slope;
#   - a non-SSR Markov control diagnostic path;
#   - bootstrap uncertainty for pairwise AUC/rank-effect separation;
#   - NMI sensitivity over multiple top-K alphabet collapses.
# v0.4.1 changes:
#   - Zipf summaries now report the effective upper rank actually entering
#     the fit after the minimum-count filter is applied;
#   - sensitivity summaries carry the effective fit window so requested
#     cutoffs are not confused with count-limited cutoffs.
#
# v0.3.1 changes:
#   - the primary descriptive Zipf slope is estimated on ranks 10--500;
#   - alpha sensitivity over alternative upper-rank cutoffs is explicit;
#   - early/late Heaps exponents are estimated by direct window regressions;
#   - local Heaps curves remain descriptive and are displayed only after n >= 1000.
#
# v0.3 changes:
#   - local Heaps exponents are estimated by sliding local regressions;
#   - Monte Carlo Heaps curves can be summarised by median and bands;
#   - pairwise mechanism separation uses an AUC-based metric in [0, 1];
#   - numerical diagnostic summaries are produced for the manuscript table.

rank_frequency <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(
    rank = seq_along(tab),
    count = as.integer(tab),
    frequency = as.numeric(tab) / sum(tab),
    type = names(tab),
    stringsAsFactors = FALSE
  )
}

.rank_fit_data <- function(data, r_min = 10L, r_max = 500L, min_count = 5L) {
  rf <- rank_frequency(data$type)
  keep <- rf$rank >= r_min & rf$rank <= r_max & rf$count >= min_count
  list(rf = rf, fit_data = rf[keep, , drop = FALSE])
}

zipf_mle_summary <- function(
    data,
    r_min = 10L,
    r_max = 500L,
    min_count = 5L,
    alpha_bounds = c(0.20, 3.00)) {
  parts <- .rank_fit_data(data, r_min = r_min, r_max = r_max, min_count = min_count)
  fit_data <- parts$fit_data
  effective_r_max <- if (nrow(fit_data)) max(fit_data$rank) else NA_integer_

  if (nrow(fit_data) < 10L || sum(fit_data$count) < 20L) {
    return(data.frame(
      alpha_mle = NA_real_, loglik_mle = NA_real_,
      n_fit_mle = nrow(fit_data), effective_r_max_mle = effective_r_max
    ))
  }

  r <- fit_data$rank
  counts <- fit_data$count
  N <- sum(counts)
  sum_log_r <- sum(counts * log(r))

  # Restricted-domain multinomial likelihood on exactly the same empirical rank
  # window as the OLS summary: q_r(alpha) propto r^{-alpha}.  This is a
  # robustness estimator for the rank-frequency exponent, not a claim that
  # the full distribution has passed a power-law goodness-of-fit test.
  negloglik <- function(alpha) {
    log_z <- log(sum(r^(-alpha)))
    alpha * sum_log_r + N * log_z
  }

  opt <- stats::optimize(
    negloglik, interval = sort(as.numeric(alpha_bounds)),
    maximum = FALSE, tol = 1e-09
  )

  data.frame(
    alpha_mle = opt$minimum,
    loglik_mle = -opt$objective,
    n_fit_mle = nrow(fit_data),
    effective_r_max_mle = effective_r_max
  )
}

zipf_summary <- function(
    data,
    r_min = 10L,
    r_max = 500L,
    min_count = 5L) {
  parts <- .rank_fit_data(data, r_min = r_min, r_max = r_max, min_count = min_count)
  rf <- parts$rf
  fit_data <- parts$fit_data

  effective_r_max <- if (nrow(fit_data)) max(fit_data$rank) else NA_integer_

  if (nrow(fit_data) < 10L) {
    return(data.frame(
      alpha = NA_real_, r2 = NA_real_, n_fit = nrow(fit_data),
      effective_r_max = effective_r_max,
      vocabulary = nrow(rf), tokens = nrow(data)
    ))
  }

  # Deliberately descriptive: T0 is a log-log slope, not a stand-alone
  # goodness-of-fit test for a power law.
  fit <- stats::lm(log(frequency) ~ log(rank), data = fit_data)
  data.frame(
    alpha = -unname(stats::coef(fit)[2L]),
    r2 = summary(fit)$r.squared,
    n_fit = nrow(fit_data),
    effective_r_max = effective_r_max,
    vocabulary = nrow(rf),
    tokens = nrow(data)
  )
}

zipf_sensitivity <- function(
    data,
    r_max_values = c(300L, 500L, 1000L, 2000L, 3000L),
    r_min = 10L,
    min_count = 5L) {
  r_max_values <- sort(unique(as.integer(r_max_values)))
  rows <- lapply(r_max_values, function(rm) {
    z <- zipf_summary(
      data,
      r_min = r_min,
      r_max = rm,
      min_count = min_count
    )
    z$r_max <- rm
    z
  })
  out <- do.call(rbind, rows)
  out$generator <- unique(data$generator)[1L]
  rownames(out) <- NULL
  out
}

summarise_alpha_sensitivity_replications <- function(
    sensitivity,
    probs = c(0.10, 0.50, 0.90)) {
  required <- c("generator", "r_max", "alpha")
  if (!all(required %in% names(sensitivity))) {
    stop("sensitivity must contain: ", paste(required, collapse = ", "))
  }
  if (length(probs) != 3L) stop("probs must contain lower, median, upper probabilities")

  key <- interaction(
    sensitivity$generator, sensitivity$r_max,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(sensitivity)), key)
  rows <- lapply(groups, function(idx) {
    x <- sensitivity$alpha[idx]
    x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8)

    eff <- if ("effective_r_max" %in% names(sensitivity)) {
      as.numeric(sensitivity$effective_r_max[idx])
    } else {
      rep(as.numeric(sensitivity$r_max[idx[1L]]), length(idx))
    }
    eff <- eff[is.finite(eff)]
    eff_q <- if (length(eff)) {
      stats::quantile(eff, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
    } else {
      rep(NA_real_, 3L)
    }

    nfit <- if ("n_fit" %in% names(sensitivity)) {
      as.numeric(sensitivity$n_fit[idx])
    } else {
      rep(NA_real_, length(idx))
    }
    nfit <- nfit[is.finite(nfit)]

    data.frame(
      generator = sensitivity$generator[idx[1L]],
      r_max = sensitivity$r_max[idx[1L]],
      alpha_lower = q[1L],
      alpha_median = q[2L],
      alpha_upper = q[3L],
      effective_r_max_lower = eff_q[1L],
      effective_r_max_median = eff_q[2L],
      effective_r_max_upper = eff_q[3L],
      n_fit_median = if (length(nfit)) stats::median(nfit) else NA_real_,
      n_rep = length(x),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$generator, out$r_max), , drop = FALSE]
}

heaps_curve <- function(data, n_points = 35L, min_n = 500L) {
  n <- nrow(data)
  min_n <- min(as.integer(min_n), n)
  checkpoints <- unique(as.integer(round(exp(seq(
    log(max(2L, min_n)), log(n), length.out = n_points
  )))))

  cumulative_vocab <- cumsum(!duplicated(data$type))
  data.frame(
    n = checkpoints,
    vocabulary = cumulative_vocab[checkpoints]
  )
}

.local_loglog_slope <- function(x, y, index, window = 7L) {
  n <- length(x)
  window <- max(3L, as.integer(window))
  if (window %% 2L == 0L) window <- window + 1L
  half <- (window - 1L) %/% 2L

  lo <- max(1L, index - half)
  hi <- min(n, index + half)

  # Keep approximately the requested window near the boundaries.
  if ((hi - lo + 1L) < window && n >= window) {
    if (lo == 1L) hi <- window
    if (hi == n) lo <- n - window + 1L
  }

  idx <- lo:hi
  ok <- is.finite(x[idx]) & is.finite(y[idx]) & x[idx] > 0 & y[idx] > 0
  idx <- idx[ok]
  if (length(idx) < 3L || length(unique(x[idx])) < 3L) return(c(NA_real_, NA_real_))

  fit <- stats::lm(log(y[idx]) ~ log(x[idx]))
  c(
    beta = unname(stats::coef(fit)[2L]),
    r2 = summary(fit)$r.squared
  )
}

heaps_local_curve <- function(
    data,
    n_points = 35L,
    min_n = 500L,
    window = 7L) {
  h <- heaps_curve(data, n_points = n_points, min_n = min_n)
  if (nrow(h) < 3L) {
    h$beta_local <- NA_real_
    h$beta_local_r2 <- NA_real_
    return(h)
  }

  local <- t(vapply(
    seq_len(nrow(h)),
    function(j) .local_loglog_slope(h$n, h$vocabulary, j, window = window),
    numeric(2)
  ))

  h$beta_local <- local[, "beta"]
  h$beta_local_r2 <- local[, "r2"]
  h
}

heaps_summary <- function(
    data,
    n_points = 35L,
    min_n = 500L,
    window = 7L) {
  hc <- heaps_local_curve(
    data,
    n_points = n_points,
    min_n = min_n,
    window = window
  )

  if (nrow(hc) < 6L || length(unique(hc$vocabulary)) < 2L) {
    return(data.frame(
      beta = NA_real_, beta_early = NA_real_, beta_late = NA_real_,
      beta_change = NA_real_, r2 = NA_real_
    ))
  }

  fit <- stats::lm(log(vocabulary) ~ log(n), data = hc)

  # Direct regressions over disjoint early and late windows are more stable
  # than averages of local derivatives and have a simple interpretation.
  q <- max(3L, floor(nrow(hc) / 3L))
  early_idx <- seq_len(q)
  late_idx <- (nrow(hc) - q + 1L):nrow(hc)

  fit_early <- stats::lm(
    log(vocabulary) ~ log(n),
    data = hc[early_idx, , drop = FALSE]
  )
  fit_late <- stats::lm(
    log(vocabulary) ~ log(n),
    data = hc[late_idx, , drop = FALSE]
  )

  beta_early <- unname(stats::coef(fit_early)[2L])
  beta_late <- unname(stats::coef(fit_late)[2L])

  data.frame(
    beta = unname(stats::coef(fit)[2L]),
    beta_early = beta_early,
    beta_late = beta_late,
    beta_change = beta_late - beta_early,
    r2 = summary(fit)$r.squared
  )
}

summarise_heaps_local_replications <- function(
    curves,
    probs = c(0.10, 0.50, 0.90)) {
  required <- c("generator", "n", "beta_local")
  if (!all(required %in% names(curves))) {
    stop("curves must contain: ", paste(required, collapse = ", "))
  }
  if (length(probs) != 3L) stop("probs must contain lower, median, upper probabilities")

  key <- interaction(curves$generator, curves$n, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(curves)), key)
  rows <- lapply(groups, function(idx) {
    x <- curves$beta_local[idx]
    x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
    data.frame(
      generator = curves$generator[idx[1L]],
      n = curves$n[idx[1L]],
      beta_lower = q[1L],
      beta_median = q[2L],
      beta_upper = q[3L],
      n_rep = length(x),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$generator, out$n), , drop = FALSE]
}

transition_direction <- function(data) {
  if (nrow(data) < 2L) {
    return(data.frame(
      toward_more_frequent = NA_real_, toward_less_frequent = NA_real_,
      tie_rate = NA_real_
    ))
  }

  rf <- rank_frequency(data$type)
  rank_map <- setNames(rf$rank, rf$type)
  r <- unname(rank_map[data$type])

  same_sequence <- data$sequence_id[-nrow(data)] == data$sequence_id[-1L]
  current <- r[-length(r)][same_sequence]
  next_r <- r[-1L][same_sequence]

  non_tie <- next_r != current
  if (sum(non_tie) == 0L) {
    return(data.frame(
      toward_more_frequent = NA_real_, toward_less_frequent = NA_real_,
      tie_rate = 1
    ))
  }

  data.frame(
    toward_more_frequent = mean(next_r[non_tie] < current[non_tie]),
    toward_less_frequent = mean(next_r[non_tie] > current[non_tie]),
    tie_rate = mean(next_r == current)
  )
}

.collapse_top_types <- function(data, top_k = 50L) {
  tab <- sort(table(data$type), decreasing = TRUE)
  top <- names(tab)[seq_len(min(top_k, length(tab)))]
  x <- ifelse(data$type %in% top, data$type, "__OTHER__")
  factor(x, levels = c(top, "__OTHER__"))
}

.normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n <= 0L) return(NA_real_)

  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(k) {
    i <- nz[k, 1L]
    j <- nz[k, 2L]
    pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
  }, numeric(1)))

  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  denom <- sqrt(hx * hy)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  mi / denom
}

transition_information <- function(data, top_k = 50L) {
  if (nrow(data) < 2L) return(NA_real_)
  x <- .collapse_top_types(data, top_k = top_k)
  same_sequence <- data$sequence_id[-nrow(data)] == data$sequence_id[-1L]
  if (!any(same_sequence)) return(NA_real_)
  current <- x[-length(x)][same_sequence]
  next_x <- x[-1L][same_sequence]
  .normalized_mutual_information(current, next_x)
}

.shuffle_within_sequences <- function(data) {
  out <- data
  groups <- split(seq_len(nrow(data)), data$sequence_id)
  for (idx in groups) {
    if (length(idx) > 1L) out$type[idx] <- sample(data$type[idx], replace = FALSE)
  }
  out
}

transition_information_excess <- function(
    data,
    top_k = 50L,
    n_surrogates = 5L,
    seed = 9001L) {
  obs <- transition_information(data, top_k = top_k)
  if (!is.finite(obs) || n_surrogates < 1L) {
    return(data.frame(
      transition_nmi = obs,
      surrogate_nmi = NA_real_,
      transition_nmi_excess = NA_real_,
      transition_nmi_z = NA_real_
    ))
  }

  set.seed(seed)
  null <- vapply(seq_len(n_surrogates), function(b) {
    transition_information(.shuffle_within_sequences(data), top_k = top_k)
  }, numeric(1))

  null_mean <- mean(null, na.rm = TRUE)
  null_sd <- stats::sd(null, na.rm = TRUE)
  z <- if (is.finite(null_sd) && null_sd > 0) (obs - null_mean) / null_sd else NA_real_

  data.frame(
    transition_nmi = obs,
    surrogate_nmi = null_mean,
    transition_nmi_excess = obs - null_mean,
    transition_nmi_z = z
  )
}

transition_information_excess_multi <- function(
    data,
    top_k_values = c(25L, 50L, 100L),
    n_surrogates = 5L,
    seed = 9001L) {
  top_k_values <- sort(unique(as.integer(top_k_values)))
  top_k_values <- top_k_values[top_k_values >= 2L]
  if (!length(top_k_values)) stop("top_k_values must contain values >= 2")

  obs <- vapply(
    top_k_values,
    function(k) transition_information(data, top_k = k),
    numeric(1)
  )

  if (n_surrogates < 1L) {
    return(data.frame(
      top_k = top_k_values, transition_nmi = obs, surrogate_nmi = NA_real_,
      transition_nmi_excess = NA_real_, transition_nmi_z = NA_real_
    ))
  }

  set.seed(seed)
  null <- matrix(NA_real_, nrow = n_surrogates, ncol = length(top_k_values))
  for (b in seq_len(n_surrogates)) {
    shuffled <- .shuffle_within_sequences(data)
    null[b, ] <- vapply(
      top_k_values,
      function(k) transition_information(shuffled, top_k = k),
      numeric(1)
    )
  }

  null_mean <- colMeans(null, na.rm = TRUE)
  null_sd <- apply(null, 2L, stats::sd, na.rm = TRUE)
  z <- ifelse(is.finite(null_sd) & null_sd > 0, (obs - null_mean) / null_sd, NA_real_)

  data.frame(
    top_k = top_k_values,
    transition_nmi = obs,
    surrogate_nmi = null_mean,
    transition_nmi_excess = obs - null_mean,
    transition_nmi_z = z
  )
}

latent_conditioning <- function(
    data,
    bins = 4L,
    r_min = 10L,
    r_max = 500L,
    min_count = 5L) {
  if (all(is.na(data$latent_scale))) return(NULL)

  seq_meta <- unique(data[, c("sequence_id", "latent_scale")])
  seq_meta$context_bin <- pmin(
    bins,
    ceiling(rank(seq_meta$latent_scale, ties.method = "first") /
              nrow(seq_meta) * bins)
  )
  bin_map <- setNames(seq_meta$context_bin, seq_meta$sequence_id)
  context_bin <- unname(bin_map[as.character(data$sequence_id)])

  out <- lapply(seq_len(bins), function(b) {
    d <- data[context_bin == b, , drop = FALSE]
    z <- zipf_summary(d, r_min = r_min, r_max = r_max,
                      min_count = min_count)
    z$context_bin <- b
    z
  })
  do.call(rbind, out)
}

length_frequency_summary <- function(data) {
  ok <- !is.na(data$unit_length)
  if (!any(ok)) return(NULL)

  d <- data[ok, c("type", "unit_length"), drop = FALSE]
  counts <- table(d$type)
  first_length <- tapply(d$unit_length, d$type, function(x) x[1L])
  common <- intersect(names(counts), names(first_length))
  df <- data.frame(
    type = common,
    count = as.numeric(counts[common]),
    unit_length = as.numeric(first_length[common])
  )
  df <- df[df$count >= 2L, , drop = FALSE]
  if (nrow(df) < 10L) return(NULL)

  fit <- stats::lm(log(count) ~ unit_length, data = df)
  data.frame(
    slope = unname(stats::coef(fit)[2L]),
    r2 = summary(fit)$r.squared,
    n_types = nrow(df)
  )
}

benchmark_diagnostics <- function(
    data,
    n_surrogates = 5L,
    surrogate_seed = 9001L,
    heaps_window = 7L,
    zipf_r_max = 500L,
    nmi_top_k = 50L,
    transition_info = NULL) {
  z <- zipf_summary(data, r_max = zipf_r_max)
  zmle <- zipf_mle_summary(data, r_max = zipf_r_max)
  h <- heaps_summary(data, window = heaps_window)
  tr <- transition_direction(data)
  ti <- if (is.null(transition_info)) {
    transition_information_excess(
      data, top_k = nmi_top_k,
      n_surrogates = n_surrogates, seed = surrogate_seed
    )
  } else {
    transition_info[1L, , drop = FALSE]
  }
  lc <- latent_conditioning(data, r_max = zipf_r_max)
  lf <- length_frequency_summary(data)

  data.frame(
    generator = unique(data$generator)[1L],
    alpha = z$alpha,
    alpha_mle = zmle$alpha_mle,
    zipf_r2 = z$r2,
    zipf_n_fit = z$n_fit,
    zipf_effective_r_max = z$effective_r_max,
    zipf_mle_n_fit = zmle$n_fit_mle,
    zipf_mle_effective_r_max = zmle$effective_r_max_mle,
    vocabulary = z$vocabulary,
    heaps_beta = h$beta,
    heaps_beta_early = h$beta_early,
    heaps_beta_late = h$beta_late,
    heaps_beta_change = h$beta_change,
    transition_toward_more_frequent = tr$toward_more_frequent,
    transition_nmi = ti$transition_nmi,
    transition_nmi_excess = ti$transition_nmi_excess,
    transition_nmi_z = ti$transition_nmi_z,
    conditional_alpha_sd = if (is.null(lc)) NA_real_ else stats::sd(lc$alpha, na.rm = TRUE),
    conditional_mean_r2 = if (is.null(lc)) NA_real_ else mean(lc$r2, na.rm = TRUE),
    length_frequency_slope = if (is.null(lf)) NA_real_ else lf$slope,
    length_frequency_r2 = if (is.null(lf)) NA_real_ else lf$r2,
    stringsAsFactors = FALSE
  )
}

.empirical_auc <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx < 1L || ny < 1L) return(NA_real_)

  pooled <- c(x, y)
  r <- rank(pooled, ties.method = "average")
  u_x <- sum(r[seq_len(nx)]) - nx * (nx + 1) / 2
  u_x / (nx * ny)
}

pairwise_auc_separation <- function(
    results,
    diagnostics = c(
      "alpha", "alpha_mle", "vocabulary", "heaps_beta", "heaps_beta_change",
      "transition_toward_more_frequent", "transition_nmi_excess"
    )) {
  gens <- unique(results$generator)
  pairs <- utils::combn(gens, 2L, simplify = FALSE)
  rows <- list()
  k <- 1L

  for (pair in pairs) {
    for (variable in diagnostics) {
      if (!variable %in% names(results)) next
      x <- results[results$generator == pair[1L], variable]
      y <- results[results$generator == pair[2L], variable]
      auc <- .empirical_auc(x, y)
      if (!is.finite(auc)) next
      rows[[k]] <- data.frame(
        pair = paste(pair, collapse = " vs "),
        diagnostic = variable,
        auc = auc,
        separation = 2 * abs(auc - 0.5),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

pairwise_auc_bootstrap <- function(
    results,
    diagnostics = c(
      "alpha", "alpha_mle", "vocabulary", "heaps_beta", "heaps_beta_change",
      "transition_toward_more_frequent", "transition_nmi_excess"
    ),
    n_boot = 1000L,
    conf_level = 0.95,
    seed = 20260808L) {
  stopifnot(n_boot >= 100L, conf_level > 0, conf_level < 1)
  gens <- unique(results$generator)
  pairs <- utils::combn(gens, 2L, simplify = FALSE)
  alpha <- (1 - conf_level) / 2
  rows <- list()
  k_out <- 1L
  set.seed(seed)

  for (pair in pairs) {
    for (variable in diagnostics) {
      if (!variable %in% names(results)) next
      x <- results[results$generator == pair[1L], variable]
      y <- results[results$generator == pair[2L], variable]
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      if (length(x) < 2L || length(y) < 2L) next

      auc <- .empirical_auc(x, y)
      rank_effect <- 2 * auc - 1
      boot <- vapply(seq_len(n_boot), function(b) {
        xb <- sample(x, length(x), replace = TRUE)
        yb <- sample(y, length(y), replace = TRUE)
        2 * .empirical_auc(xb, yb) - 1
      }, numeric(1))
      q <- stats::quantile(
        boot, probs = c(alpha, 1 - alpha),
        na.rm = TRUE, names = FALSE, type = 8
      )

      rows[[k_out]] <- data.frame(
        pair = paste(pair, collapse = " vs "),
        diagnostic = variable,
        auc = auc,
        rank_effect = rank_effect,
        separation = abs(rank_effect),
        rank_effect_lower = q[1L],
        rank_effect_upper = q[2L],
        ci_contains_zero = q[1L] <= 0 & q[2L] >= 0,
        n_x = length(x),
        n_y = length(y),
        n_boot = n_boot,
        conf_level = conf_level,
        stringsAsFactors = FALSE
      )
      k_out <- k_out + 1L
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Backward-compatible alias kept so old notebooks fail gracefully rather than
# silently changing interpretation. v0.3 code should call pairwise_auc_separation().
pairwise_standardized_separation <- function(...) {
  warning("v0.3 replaces SMD by bounded AUC separation; returning AUC separation.")
  pairwise_auc_separation(...)
}

.quantile_safe <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, p, names = FALSE, type = 8))
}

diagnostic_summary_table <- function(
    results,
    diagnostics = c(
      "alpha", "alpha_mle", "vocabulary", "heaps_beta", "heaps_beta_change",
      "transition_nmi_excess", "transition_toward_more_frequent"
    )) {
  gens <- unique(results$generator)
  rows <- list()
  k <- 1L
  for (g in gens) {
    d <- results[results$generator == g, , drop = FALSE]
    for (variable in diagnostics) {
      if (!variable %in% names(d)) next
      x <- d[[variable]]
      x <- x[is.finite(x)]
      if (!length(x)) next
      rows[[k]] <- data.frame(
        generator = g,
        diagnostic = variable,
        mean = mean(x),
        sd = if (length(x) >= 2L) stats::sd(x) else NA_real_,
        median = stats::median(x),
        q10 = .quantile_safe(x, 0.10),
        q90 = .quantile_safe(x, 0.90),
        n = length(x),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
