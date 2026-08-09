# Generators for the Zipf equifinality benchmark.
# All functions return the same core token-level schema:
# token_index, type, sequence_id, generator, plus optional metadata.



simulate_maxent_zipf <- function(
    n_tokens = 100000L,
    vocabulary_size = 5000L,
    exponent = 1,
    block_size = 250L,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(vocabulary_size >= 2L, exponent > 0, block_size >= 2L)

  ranks <- seq_len(vocabulary_size)
  weights <- ranks^(-exponent)
  prob <- weights / sum(weights)
  token <- sample.int(
    vocabulary_size, size = n_tokens, replace = TRUE, prob = prob
  )

  data.frame(
    token_index = seq_len(n_tokens),
    type = sprintf("m%06d", token),
    sequence_id = ceiling(seq_len(n_tokens) / block_size),
    generator = "finite Zipf reference",
    unit_length = NA_integer_,
    latent_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}


simulate_markov_zipf <- function(
    n_tokens = 100000L,
    vocabulary_size = 5000L,
    exponent = 1,
    persistence = 0.35,
    block_size = 250L,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(vocabulary_size >= 2L, exponent > 0)
  stopifnot(persistence >= 0, persistence < 1, block_size >= 2L)

  ranks <- seq_len(vocabulary_size)
  weights <- ranks^(-exponent)
  prob <- weights / sum(weights)

  # Transition kernel:
  #   P_ij = persistence * 1{i=j} + (1-persistence) * p_j.
  # The finite Zipf marginal p is stationary. Sequence boundaries are
  # restarted from p so each block also starts in stationarity.
  sequence_id <- ceiling(seq_len(n_tokens) / block_size)
  is_boundary <- c(TRUE, sequence_id[-1L] != sequence_id[-n_tokens])
  refresh <- is_boundary | (stats::runif(n_tokens) > persistence)
  refresh[1L] <- TRUE

  refresh_states <- sample.int(
    vocabulary_size, size = sum(refresh), replace = TRUE, prob = prob
  )
  token <- refresh_states[cumsum(refresh)]

  data.frame(
    token_index = seq_len(n_tokens),
    type = sprintf("k%06d", token),
    sequence_id = sequence_id,
    generator = "Zipf Markov control",
    unit_length = NA_integer_,
    latent_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}

.resample_positive_geometric <- function(n, prob) {
  x <- stats::rgeom(n, prob = prob)
  zero <- which(x == 0L)
  while (length(zero) > 0L) {
    x[zero] <- stats::rgeom(length(zero), prob = prob)
    zero <- which(x == 0L)
  }
  x
}

simulate_random_segmentation <- function(
    n_tokens = 100000L,
    alphabet_size = 8L,
    p_delim = 0.62,
    letter_decay = 0.05,
    block_size = 250L,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(alphabet_size >= 2L, alphabet_size <= 26L)
  stopifnot(p_delim > 0, p_delim < 1)

  alphabet <- letters[seq_len(alphabet_size)]
  weights <- exp(-letter_decay * (seq_len(alphabet_size) - 1L))
  letter_prob <- weights / sum(weights)

  # rgeom gives the number of non-delimiter symbols before the delimiter.
  # Empty units are discarded by conditioning on length >= 1.
  unit_length <- .resample_positive_geometric(n_tokens, p_delim)

  type <- vapply(
    unit_length,
    function(L) paste0(sample(alphabet, size = L, replace = TRUE,
                             prob = letter_prob), collapse = ""),
    character(1)
  )

  data.frame(
    token_index = seq_len(n_tokens),
    type = type,
    sequence_id = ceiling(seq_len(n_tokens) / block_size),
    generator = "random segmentation",
    unit_length = unit_length,
    latent_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}

simulate_simon <- function(
    n_tokens = 100000L,
    innovation_prob = 0.05,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(innovation_prob > 0, innovation_prob < 1)

  token <- integer(n_tokens)
  token[1L] <- 1L
  next_type <- 2L

  if (n_tokens >= 2L) {
    for (t in 2L:n_tokens) {
      if (stats::runif(1L) < innovation_prob) {
        token[t] <- next_type
        next_type <- next_type + 1L
      } else {
        # Sampling a previous token position uniformly is exactly
        # proportional-to-current-abundance reuse.
        token[t] <- token[sample.int(t - 1L, size = 1L)]
      }
    }
  }

  data.frame(
    token_index = seq_len(n_tokens),
    type = sprintf("s%06d", token),
    sequence_id = 1L,
    generator = "Simon reinforcement",
    unit_length = NA_integer_,
    latent_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}

.rtrunc_discrete_exponential <- function(n, scale, max_rank) {
  # P(R <= r | s) = [1-exp(-r/s)] / [1-exp(-V/s)], r=1,...,V.
  u <- stats::runif(n)
  denom <- 1 - exp(-max_rank / scale)
  r <- ceiling(-scale * log(1 - u * denom))
  pmin(pmax(as.integer(r), 1L), as.integer(max_rank))
}

simulate_latent_mixture <- function(
    n_sequences = 400L,
    sequence_length = 250L,
    vocabulary_size = 5000L,
    scale_min = 1,
    scale_max = 8000,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(scale_min > 0, scale_max > scale_min)

  scales <- exp(stats::runif(
    n_sequences,
    min = log(scale_min), max = log(scale_max)
  ))

  token_list <- vector("list", n_sequences)
  for (j in seq_len(n_sequences)) {
    r <- .rtrunc_discrete_exponential(
      n = sequence_length,
      scale = scales[j],
      max_rank = vocabulary_size
    )
    token_list[[j]] <- data.frame(
      type = sprintf("l%06d", r),
      sequence_id = j,
      latent_scale = scales[j],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, token_list)
  rownames(out) <- NULL
  out$token_index <- seq_len(nrow(out))
  out$generator <- "latent mixture"
  out$unit_length <- NA_integer_
  out[, c("token_index", "type", "sequence_id", "generator",
          "unit_length", "latent_scale")]
}

simulate_ssr <- function(
    n_tokens = 100000L,
    vocabulary_size = 5000L,
    initial_distribution = c("stationary", "uniform_restart"),
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(vocabulary_size >= 2L)
  initial_distribution <- match.arg(initial_distribution)

  # Canonical SSR transition rule:
  #   i > 1: X_{t+1} ~ Uniform{1,...,i-1}
  #   i = 1: X_{t+1} ~ Uniform{1,...,V}
  # Its stationary distribution is exactly pi_j = 1/(j H_V).
  if (initial_distribution == "stationary") {
    p <- 1 / seq_len(vocabulary_size)
    p <- p / sum(p)
    current <- sample.int(vocabulary_size, size = 1L, prob = p)
  } else {
    current <- sample.int(vocabulary_size, size = 1L)
  }

  state <- integer(n_tokens)
  sequence_id <- integer(n_tokens)
  seq_id <- 1L

  for (t in seq_len(n_tokens)) {
    state[t] <- current
    sequence_id[t] <- seq_id

    if (current <= 1L) {
      # Canonical restart from the full sample space. sequence_id marks
      # natural SSR cycles so restart transitions can be excluded explicitly
      # from within-sequence diagnostics when desired.
      seq_id <- seq_id + 1L
      current <- sample.int(vocabulary_size, size = 1L)
    } else {
      current <- sample.int(current - 1L, size = 1L)
    }
  }

  data.frame(
    token_index = seq_len(n_tokens),
    type = sprintf("r%06d", state),
    sequence_id = sequence_id,
    generator = "sample-space reduction",
    unit_length = NA_integer_,
    latent_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}
