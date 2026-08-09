# Figure helpers for v0.5. Requires ggplot2; the multipanel paper figure
# additionally uses patchwork.

plot_rank_frequency <- function(simulations) {
  dfs <- lapply(simulations, function(d) {
    rf <- rank_frequency(d$type)
    rf$generator <- unique(d$generator)[1L]
    rf
  })
  all <- do.call(rbind, dfs)

  ggplot2::ggplot(all, ggplot2::aes(rank, frequency, linetype = generator)) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Rank", y = "Token frequency", linetype = NULL,
      title = "Distinct generative constructions, similar rank-frequency scaling"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_alpha_replications <- function(results, title = NULL) {
  if (is.null(title)) title <- "Equifinality under the marginal exponent"
  ggplot2::ggplot(results, ggplot2::aes(x = generator, y = alpha)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.35) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = expression(hat(alpha)), title = title) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_alpha_estimator_comparison <- function(results) {
  needed <- c("generator", "alpha", "alpha_mle")
  if (!all(needed %in% names(results))) {
    stop("results must contain: ", paste(needed, collapse = ", "))
  }
  d <- rbind(
    data.frame(generator = results$generator, estimator = "OLS log-log", alpha = results$alpha),
    data.frame(generator = results$generator, estimator = "conditional MLE", alpha = results$alpha_mle)
  )
  d <- d[is.finite(d$alpha), , drop = FALSE]

  ggplot2::ggplot(d, ggplot2::aes(x = generator, y = alpha, linetype = estimator)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = interaction(generator, estimator)),
                          position = ggplot2::position_dodge(width = 0.7),
                          outlier.shape = NA) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = expression(hat(alpha)), linetype = NULL,
      title = "Equifinality is checked with both OLS and conditional likelihood"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_alpha_sensitivity <- function(
    sensitivity,
    probs = c(0.10, 0.50, 0.90),
    primary_r_max = 500L,
    truncation_ratio = 0.98) {
  d <- summarise_alpha_sensitivity_replications(sensitivity, probs = probs)
  d$fit_window <- ifelse(
    is.finite(d$effective_r_max_median) &
      d$effective_r_max_median < truncation_ratio * d$r_max,
    "count-limited",
    "requested cutoff reached"
  )

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = r_max, y = alpha_median, linetype = generator, group = generator)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = alpha_lower, ymax = alpha_upper, fill = generator),
      alpha = 0.12, linetype = 0, show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(shape = fit_window),
      size = 2.0, stroke = 0.8
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::geom_vline(xintercept = primary_r_max, linetype = 3) +
    ggplot2::scale_x_log10(breaks = sort(unique(d$r_max))) +
    ggplot2::scale_shape_manual(
      values = c("requested cutoff reached" = 16, "count-limited" = 1),
      breaks = c("requested cutoff reached", "count-limited")
    ) +
    ggplot2::labs(
      x = expression(r[max]), y = expression(hat(alpha)),
      linetype = NULL, shape = "Effective fit window",
      title = "The fitted Zipf exponent depends on the rank window",
      subtitle = paste0(
        "Median with 10–90% held-out Monte Carlo bands; dotted line marks the primary window. ",
        "Open points indicate that count >= 5 truncates the requested upper rank."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_alpha_vocabulary_tradeoff <- function(results) {
  ggplot2::ggplot(
    results,
    ggplot2::aes(x = vocabulary, y = alpha, shape = generator)
  ) +
    ggplot2::geom_point(alpha = 0.55) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::geom_vline(xintercept = 5000, linetype = 2) +
    ggplot2::labs(
      x = "Final vocabulary size", y = expression(hat(alpha)), shape = NULL,
      title = "Matching Zipf and vocabulary size can impose incompatible constraints"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_heaps_curves <- function(simulations) {
  dfs <- lapply(simulations, function(d) {
    h <- heaps_curve(d)
    h$generator <- unique(d$generator)[1L]
    h
  })
  all <- do.call(rbind, dfs)

  ggplot2::ggplot(all, ggplot2::aes(n, vocabulary, linetype = generator)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Tokens observed", y = "Distinct types",
      linetype = NULL, title = "Vocabulary growth separates mechanisms"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_heaps_local_bands <- function(
    curves,
    probs = c(0.10, 0.50, 0.90),
    min_display_n = 1000L) {
  d <- summarise_heaps_local_replications(curves, probs = probs)
  d <- d[d$n >= min_display_n, , drop = FALSE]

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = n, y = beta_median, linetype = generator, group = generator)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = beta_lower, ymax = beta_upper, fill = generator),
      alpha = 0.12, linetype = 0, show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Tokens observed", y = expression(beta(n)), linetype = NULL,
      title = "Local Heaps exponents reveal distinct accumulation dynamics",
      subtitle = "Median local slope with 10–90% Monte Carlo bands"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

# Kept as a convenience for one-run exploratory checks.
plot_heaps_local <- function(simulations, window = 7L) {
  dfs <- lapply(simulations, function(d) {
    h <- heaps_local_curve(d, window = window)
    h$generator <- unique(d$generator)[1L]
    h
  })
  all <- do.call(rbind, dfs)

  ggplot2::ggplot(all, ggplot2::aes(n, beta_local, linetype = generator)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Tokens observed", y = expression(beta(n)), linetype = NULL,
      title = "Smoothed local Heaps exponents"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_latent_conditioning <- function(data, bins = 4L) {
  stopifnot(!all(is.na(data$latent_scale)))
  seq_meta <- unique(data[, c("sequence_id", "latent_scale")])
  seq_meta$context_bin <- pmin(
    bins,
    ceiling(rank(seq_meta$latent_scale, ties.method = "first") /
              nrow(seq_meta) * bins)
  )
  bin_map <- setNames(seq_meta$context_bin, seq_meta$sequence_id)
  d <- data
  d$context_bin <- unname(bin_map[as.character(d$sequence_id)])

  pooled <- rank_frequency(d$type)
  pooled$group <- "pooled"
  conditional <- do.call(rbind, lapply(seq_len(bins), function(b) {
    x <- rank_frequency(d$type[d$context_bin == b])
    x$group <- paste0("latent-scale Q", b)
    x
  }))
  all <- rbind(pooled, conditional)

  ggplot2::ggplot(all, ggplot2::aes(rank, frequency, linetype = group)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Rank", y = "Frequency", linetype = NULL,
      title = "Pooling can create scaling absent within latent regimes"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_transition_direction <- function(simulations) {
  out <- do.call(rbind, lapply(simulations, function(d) {
    tr <- transition_direction(d)
    data.frame(
      generator = unique(d$generator)[1L],
      toward_more_frequent = tr$toward_more_frequent
    )
  }))

  ggplot2::ggplot(out, ggplot2::aes(generator, toward_more_frequent)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = NULL,
      y = "P(next token is more frequent | rank changes)",
      title = "Sample-space reduction leaves a directional transition signature"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_transition_information <- function(results) {
  ggplot2::ggplot(
    results,
    ggplot2::aes(x = generator, y = transition_nmi_excess)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.35) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Lag-1 NMI minus within-sequence shuffle expectation",
      title = "Order-sensitive information detects sequential dependence"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_nmi_sensitivity <- function(
    nmi_results,
    probs = c(0.10, 0.50, 0.90)) {
  required <- c("generator", "top_k", "transition_nmi_excess")
  if (!all(required %in% names(nmi_results))) {
    stop("nmi_results must contain: ", paste(required, collapse = ", "))
  }

  key <- interaction(nmi_results$generator, nmi_results$top_k, drop = TRUE)
  groups <- split(seq_len(nrow(nmi_results)), key)
  rows <- lapply(groups, function(idx) {
    x <- nmi_results$transition_nmi_excess[idx]
    x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
    data.frame(
      generator = nmi_results$generator[idx[1L]],
      top_k = nmi_results$top_k[idx[1L]],
      lower = q[1L], median = q[2L], upper = q[3L],
      stringsAsFactors = FALSE
    )
  })
  d <- do.call(rbind, rows)

  ggplot2::ggplot(
    d, ggplot2::aes(x = top_k, y = median, linetype = generator, group = generator)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = generator),
      alpha = 0.12, linetype = 0, show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$top_k))) +
    ggplot2::labs(
      x = "Number of retained frequent types K",
      y = "Excess lag-1 NMI", linetype = NULL,
      title = "Sequence-information conclusions are stable to alphabet collapse",
      subtitle = "Median with 10–90% held-out Monte Carlo bands"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_bootstrap_rank_effect <- function(bootstrap_results, diagnostics = c("alpha", "alpha_mle")) {
  d <- bootstrap_results[bootstrap_results$diagnostic %in% diagnostics, , drop = FALSE]
  if (!nrow(d)) stop("No requested diagnostics in bootstrap_results")
  d$diagnostic_label <- .pretty_diagnostic(d$diagnostic)
  d$pair_label <- .pretty_pair(d$pair)

  ggplot2::ggplot(
    d, ggplot2::aes(x = rank_effect, y = pair_label)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = rank_effect_lower, xend = rank_effect_upper,
        y = pair_label, yend = pair_label
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ diagnostic_label, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      x = "Signed rank effect 2 AUC - 1 (bootstrap interval)", y = NULL,
      title = "Uncertainty in exponent-based discrimination"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

.pretty_diagnostic <- function(x) {
  labels <- c(
    alpha = "Zipf exponent (OLS)",
    alpha_mle = "Zipf exponent (MLE)",
    vocabulary = "Final vocabulary",
    heaps_beta = "Global Heaps beta",
    heaps_beta_change = "Heaps beta change",
    transition_toward_more_frequent = "Transition direction",
    transition_nmi_excess = "Excess lag-1 NMI"
  )
  out <- unname(labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

.pretty_pair <- function(x) {
  out <- x
  out <- gsub("sample-space reduction", "SSR", out, fixed = TRUE)
  out <- gsub("random segmentation", "random", out, fixed = TRUE)
  out <- gsub("Simon reinforcement", "Simon", out, fixed = TRUE)
  out <- gsub("latent mixture", "latent", out, fixed = TRUE)
  out <- gsub("MaxEnt Zipf", "finite Zipf", out, fixed = TRUE)
  out <- gsub("finite Zipf reference", "finite Zipf", out, fixed = TRUE)
  out <- gsub("Zipf Markov control", "Markov", out, fixed = TRUE)
  out
}

plot_separation_heatmap_auc <- function(separation) {
  if (!nrow(separation)) stop("No separation values available")
  d <- separation
  d$diagnostic_label <- .pretty_diagnostic(d$diagnostic)
  d$pair_label <- .pretty_pair(d$pair)

  ggplot2::ggplot(d, ggplot2::aes(diagnostic_label, pair_label, fill = separation)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", separation)),
      size = 3
    ) +
    ggplot2::scale_fill_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "Separation",
      title = "Auxiliary diagnostics separate mechanisms hidden by the Zipf exponent",
      subtitle = "2|AUC - 0.5|: 0 = indistinguishable, 1 = perfect separation"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

# Backward-compatible name for existing calls.
plot_separation_heatmap <- function(separation, ...) {
  plot_separation_heatmap_auc(separation)
}

plot_calibration_frontier <- function(calibration, model = "simon") {
  d <- calibration$tables[[model]]
  if (is.null(d)) stop("Model not found in calibration object: ", model)
  ggplot2::ggplot(d, ggplot2::aes(mean_vocabulary, mean_alpha, shape = stage)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = calibration$target_alpha, linetype = 2) +
    ggplot2::geom_vline(xintercept = calibration$target_vocabulary, linetype = 2) +
    ggplot2::labs(
      x = "Mean final vocabulary on calibration seeds",
      y = expression(E(hat(alpha))), shape = NULL,
      title = paste("Calibration frontier:", model)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_main_figure <- function(
    p_rank,
    p_alpha,
    p_heaps,
    p_latent,
    p_mi,
    p_sep) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install package 'patchwork' to build the multipanel paper figure.")
  }

  compact_theme <- ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7),
    legend.key.width = grid::unit(1.0, "lines"),
    plot.margin = ggplot2::margin(5, 8, 5, 5)
  )

  a <- p_rank +
    ggplot2::labs(title = "Distinct generative constructions, the same Zipf target") +
    compact_theme +
    ggplot2::theme(legend.position = "bottom")
  b <- p_alpha +
    ggplot2::labs(title = "Marginal equifinality under the OLS exponent") +
    compact_theme
  c <- p_heaps +
    ggplot2::labs(title = "Vocabulary accumulation diverges") +
    compact_theme +
    ggplot2::theme(legend.position = "none")
  d <- p_latent +
    ggplot2::labs(title = "Pooling creates marginal scaling") +
    compact_theme +
    ggplot2::theme(legend.position = "bottom")
  e <- p_mi +
    ggplot2::labs(title = "Sequence information reveals temporal structure") +
    compact_theme
  f <- p_sep +
    ggplot2::labs(
      title = "Auxiliary diagnostics recover generative information",
      subtitle = "2|AUC - 0.5|: 0 = none, 1 = perfect"
    ) +
    compact_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 7),
      legend.position = "right"
    )

  patchwork::wrap_plots(
    a, b, c, d, e, f,
    ncol = 2,
    widths = c(1, 1.08)
  ) + patchwork::plot_annotation(tag_levels = "A")
}

plot_observation_design_figure <- function(p_sensitivity, p_tradeoff) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install package 'patchwork' to build the observation-design figure.")
  }
  compact <- ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7),
    plot.margin = ggplot2::margin(5, 8, 5, 5)
  )
  a <- p_sensitivity +
    ggplot2::labs(title = "The apparent exponent depends on the rank window") +
    compact
  b <- p_tradeoff +
    ggplot2::labs(title = "A second marginal constraint shrinks the model class") +
    compact
  patchwork::wrap_plots(a, b, ncol = 2, widths = c(1.08, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
}

# -------------------------------------------------------------------------
# v0.5.1 paper-ready figure helpers
# -------------------------------------------------------------------------

plot_transition_direction_replications <- function(results) {
  needed <- c("generator", "transition_toward_more_frequent")
  if (!all(needed %in% names(results))) {
    stop("results must contain: ", paste(needed, collapse = ", "))
  }
  d <- results[is.finite(results$transition_toward_more_frequent), , drop = FALSE]
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = generator, y = transition_toward_more_frequent)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.8) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = NULL,
      y = "P(next rank is more frequent | rank changes)",
      title = "Directionality distinguishes SSR from generic Markov dependence"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_sequence_diagnostic_figure <- function(results) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install package 'patchwork' to build the sequence diagnostic figure.")
  }
  p1 <- plot_transition_information(results) +
    ggplot2::labs(title = "A. Excess lag-1 information detects order dependence")
  p2 <- plot_transition_direction_replications(results) +
    ggplot2::labs(title = "B. Directionality isolates SSR-like contraction")
  p1 + p2 + patchwork::plot_layout(widths = c(1, 1))
}

plot_bootstrap_selected_diagnostics <- function(
    bootstrap_results,
    diagnostics = c(
      "alpha", "alpha_mle",
      "transition_nmi_excess", "transition_toward_more_frequent"
    )) {
  d <- bootstrap_results[bootstrap_results$diagnostic %in% diagnostics, , drop = FALSE]
  if (!nrow(d)) stop("No requested diagnostics in bootstrap_results")
  d$diagnostic_label <- .pretty_diagnostic(d$diagnostic)
  d$pair_label <- .pretty_pair(d$pair)
  d$diagnostic_label <- factor(
    d$diagnostic_label,
    levels = .pretty_diagnostic(diagnostics)
  )

  ggplot2::ggplot(d, ggplot2::aes(x = rank_effect, y = pair_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = rank_effect_lower, xmax = rank_effect_upper),
      height = 0.18
    ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ diagnostic_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      x = "Signed rank effect 2 AUC - 1 (95% bootstrap interval)",
      y = NULL,
      title = "Bootstrap uncertainty separates weak marginal evidence from strong structural diagnostics"
    ) +
    ggplot2::theme_minimal(base_size = 10)
}

plot_joint_sequence_robustness <- function(results) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install package 'patchwork' to build the joint robustness figure.")
  }
  p1 <- plot_alpha_vocabulary_tradeoff(results) +
    ggplot2::labs(title = "A. Joint matching on exponent and final vocabulary")
  p2 <- plot_transition_information(results) +
    ggplot2::labs(title = "B. Sequence information persists after joint matching")
  p1 + p2 + patchwork::plot_layout(widths = c(1, 1))
}
