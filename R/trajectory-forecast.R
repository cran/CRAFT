trajectory_forecast <- function(
  trajectory_draws,
  last_level = NULL,
  recon_mode = NULL,  # NULL | "differences" | "scaled_growth" | "log_growth"
  
  # smooth_quantile_distribution params
  probs = seq(0.001, 0.999, length.out = 1000),
  spar = NULL,
  density_spar = NULL,
  eps = 1e-6,
  enforce_monotone = TRUE,
  plot_result = FALSE,
  
  # robustness params
  min_obs = 10L,
  min_unique = 10L,
  jitter_if_needed = TRUE,
  jitter_scale = 1e-8,
  zero_density_outside = TRUE,
  
  seed = 1L,
  verbose = FALSE,
  return_diagnostics = TRUE
) {
  stopifnot(is.data.frame(trajectory_draws))
  if (ncol(trajectory_draws) < 1L) stop("`trajectory_draws` must have at least one column.")
  
  if (!exists("smooth_quantile_distribution", mode = "function")) {
    stop("smooth_quantile_distribution() must be defined before calling trajectory_forecast().")
  }
  
  if (!is.null(recon_mode) && is.null(last_level)) {
    stop("If recon_mode is provided, you must also provide last_level.")
  }
  
  if (!is.null(recon_mode) &&
      !recon_mode %in% c("differences", "scaled_growth", "log_growth")) {
    stop("recon_mode must be one of: 'differences', 'scaled_growth', 'log_growth' or NULL.")
  }
  
  min_obs <- max(10L, as.integer(min_obs))
  min_unique <- max(3L, as.integer(min_unique))
  jitter_scale <- max(0, jitter_scale)
  
  if (!is.null(last_level)) {
    last_level <- as.numeric(last_level)
    
    if (length(last_level) == 1L) {
      last_level <- rep(last_level, ncol(trajectory_draws))
    } else if (length(last_level) != ncol(trajectory_draws)) {
      stop("last_level must be either a scalar or a vector with length equal to ncol(trajectory_draws).")
    }
    
    if (any(!is.finite(last_level))) {
      stop("last_level must contain only finite values.")
    }
    
    if (!is.null(recon_mode) && recon_mode == "log_growth" && any(last_level <= 0)) {
      stop("last_level must be > 0 for 'log_growth'.")
    }
    
    if (!is.null(recon_mode) && recon_mode == "scaled_growth" && any(last_level == 0)) {
      stop("last_level cannot contain 0 for 'scaled_growth'.")
    }
  }
  
  # ------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------
  
  .safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(1e-8)
    
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) s <- stats::mad(x)
    if (!is.finite(s) || s <= 0) s <- 1e-8
    
    max(s, 1e-8)
  }
  
  .sample_skew <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3L) return(NA_real_)
    
    s <- .safe_sd(x)
    mean((x - mean(x))^3) / (s^3 + 1e-12)
  }
  
  .sample_excess_kurtosis <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 4L) return(NA_real_)
    
    s <- .safe_sd(x)
    mean((x - mean(x))^4) / (s^4 + 1e-12) - 3
  }
  
  .reconstruct_to_level <- function(x, last_level_j, recon_mode) {
    x <- as.numeric(x)
    
    if (is.null(recon_mode)) {
      return(x)
    }
    
    if (recon_mode == "differences") {
      return(last_level_j + x)
    }
    
    if (recon_mode == "scaled_growth") {
      return(last_level_j * (1 + x))
    }
    
    if (recon_mode == "log_growth") {
      return(last_level_j * exp(x))
    }
    
    stop("Unknown recon_mode.")
  }
  
  .empty_distribution <- function() {
    list(
      data = data.frame(),
      dfun = function(z) rep(0, length(z)),
      pfun = function(z) rep(NA_real_, length(z)),
      qfun = function(p) rep(NA_real_, length(p)),
      rfun = function(n) rep(NA_real_, n),
      q_fit = NULL,
      d_fit = NULL,
      original_x = numeric(0),
      fallback = "empty"
    )
  }
  
  .gaussian_fallback_distribution <- function(x, reason = "fallback_gaussian") {
    x <- x[is.finite(x)]
    
    mu <- if (length(x)) mean(x) else 0
    sig <- .safe_sd(x)
    
    grid <- seq(mu - 4 * sig, mu + 4 * sig, length.out = 512L)
    
    dfun <- function(z) {
      stats::dnorm(as.numeric(z), mean = mu, sd = sig)
    }
    
    pfun <- function(z) {
      stats::pnorm(as.numeric(z), mean = mu, sd = sig)
    }
    
    qfun <- function(p) {
      p <- pmin(pmax(as.numeric(p), eps), 1 - eps)
      stats::qnorm(p, mean = mu, sd = sig)
    }
    
    rfun <- function(n) {
      if (length(n) != 1 || n <= 0) stop("n must be a positive integer.")
      stats::rnorm(as.integer(n), mean = mu, sd = sig)
    }
    
    list(
      data = data.frame(
        x = grid,
        density = dfun(grid),
        cdf = pfun(grid),
        q_smooth = grid
      ),
      dfun = dfun,
      pfun = pfun,
      qfun = qfun,
      rfun = rfun,
      q_fit = NULL,
      d_fit = NULL,
      original_x = x,
      fallback = reason
    )
  }
  
  .prepare_for_smoothing <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    
    info <- list(
      original_n = length(x),
      augmented = FALSE,
      jittered = FALSE
    )
    
    if (!length(x)) {
      return(list(x = x, info = info))
    }
    
    # If too short, bootstrap with replacement to reach min_obs
    if (length(x) < min_obs) {
      x <- sample(x, size = min_obs, replace = TRUE)
      info$augmented <- TRUE
    }
    
    # If too discrete / constant, add tiny jitter
    if (jitter_if_needed && length(unique(x)) < min_unique) {
      sc <- max(.safe_sd(x), abs(stats::median(x)), 1)
      x <- x + stats::rnorm(length(x), mean = 0, sd = jitter_scale * sc)
      info$jittered <- TRUE
    }
    
    list(x = x, info = info)
  }
  
  .zero_wrap_density <- function(dist) {
    if (!zero_density_outside) return(dist)
    if (!is.data.frame(dist$data) || !nrow(dist$data) || !"x" %in% names(dist$data)) {
      return(dist)
    }
    
    lo <- min(dist$data$x, na.rm = TRUE)
    hi <- max(dist$data$x, na.rm = TRUE)
    old_dfun <- dist$dfun
    
    dist$dfun <- function(z) {
      z <- as.numeric(z)
      out <- old_dfun(z)
      out[!is.finite(z) | z < lo | z > hi] <- 0
      pmax(0, out)
    }
    
    dist
  }
  
  .fit_smooth_distribution <- function(x) {
    x_raw <- as.numeric(x)
    x_raw <- x_raw[is.finite(x_raw)]
    
    if (!length(x_raw)) {
      return(list(
        dist = .empty_distribution(),
        diagnostics = list(
          original_n = 0L,
          fitted_n = 0L,
          augmented = FALSE,
          jittered = FALSE,
          fallback = "empty",
          support = c(NA_real_, NA_real_),
          mean = NA_real_,
          median = NA_real_,
          sd = NA_real_,
          sample_skew = NA_real_,
          sample_excess_kurtosis = NA_real_
        )
      ))
    }
    
    prep <- .prepare_for_smoothing(x_raw)
    x_fit <- prep$x
    
    fit <- try(
      smooth_quantile_distribution(
        x = x_fit,
        probs = probs,
        spar = spar,
        density_spar = density_spar,
        eps = eps,
        enforce_monotone = enforce_monotone,
        plot_result = plot_result
      ),
      silent = TRUE
    )
    
    fallback_reason <- NA_character_
    
    # Second attempt with stronger jitter if the first smoothing fails
    if (inherits(fit, "try-error") && jitter_if_needed) {
      sc <- max(.safe_sd(x_fit), abs(stats::median(x_fit)), 1)
      x_retry <- x_fit + stats::rnorm(length(x_fit), mean = 0, sd = 1000 * jitter_scale * sc)
      
      fit <- try(
        smooth_quantile_distribution(
          x = x_retry,
          probs = probs,
          spar = spar,
          density_spar = density_spar,
          eps = eps,
          enforce_monotone = enforce_monotone,
          plot_result = plot_result
        ),
        silent = TRUE
      )
      
      prep$info$jittered <- TRUE
    }
    
    if (inherits(fit, "try-error")) {
      fallback_reason <- "gaussian_after_smoothing_failure"
      fit <- .gaussian_fallback_distribution(x_raw, reason = fallback_reason)
    } else {
      fit <- .zero_wrap_density(fit)
      fit$fallback <- NA_character_
    }
    
    sx <- if (is.data.frame(fit$data) && nrow(fit$data) && "x" %in% names(fit$data)) {
      range(fit$data$x, na.rm = TRUE)
    } else {
      range(x_raw, na.rm = TRUE)
    }
    
    list(
      dist = fit,
      diagnostics = list(
        original_n = length(x_raw),
        fitted_n = length(x_fit),
        augmented = prep$info$augmented,
        jittered = prep$info$jittered,
        fallback = fit$fallback,
        support = sx,
        mean = mean(x_raw),
        median = stats::median(x_raw),
        sd = .safe_sd(x_raw),
        sample_skew = .sample_skew(x_raw),
        sample_excess_kurtosis = .sample_excess_kurtosis(x_raw)
      )
    )
  }
  
  # ------------------------------------------------------------
  # Main loop
  # ------------------------------------------------------------
  
  set.seed(as.integer(seed))
  
  out <- vector("list", ncol(trajectory_draws))
  diagnostics <- vector("list", ncol(trajectory_draws))
  
  out_names <- colnames(trajectory_draws)
  if (is.null(out_names) || any(!nzchar(out_names))) {
    out_names <- paste0("t", seq_len(ncol(trajectory_draws)))
  }
  
  names(out) <- out_names
  names(diagnostics) <- out_names
  
  for (j in seq_len(ncol(trajectory_draws))) {
    if (verbose) {
      message("Fitting smooth empirical distribution for column ", out_names[j],
              " [", j, " / ", ncol(trajectory_draws), "]")
    }
    
    xj_original <- as.numeric(trajectory_draws[[j]])
    xj_original <- xj_original[is.finite(xj_original)]
    
    xj_target <- if (!is.null(recon_mode)) {
      .reconstruct_to_level(
        x = xj_original,
        last_level_j = last_level[j],
        recon_mode = recon_mode
      )
    } else {
      xj_original
    }
    
    xj_target <- xj_target[is.finite(xj_target)]
    
    fitted <- .fit_smooth_distribution(xj_target)
    
    out[[j]] <- fitted$dist
    
    diagnostics[[j]] <- c(
      list(
        column = out_names[j],
        recon_mode = recon_mode,
        last_level = if (!is.null(last_level)) last_level[j] else NULL,
        original_input_support = if (length(xj_original)) range(xj_original) else c(NA_real_, NA_real_),
        target_scale = if (is.null(recon_mode)) "input" else "level"
      ),
      fitted$diagnostics
    )
  }
  
  if (return_diagnostics) {
    attr(out, "diagnostics") <- diagnostics
  }
  
  attr(out, "settings") <- list(
    probs = probs,
    spar = spar,
    density_spar = density_spar,
    eps = eps,
    enforce_monotone = enforce_monotone,
    min_obs = min_obs,
    min_unique = min_unique,
    jitter_if_needed = jitter_if_needed,
    jitter_scale = jitter_scale,
    zero_density_outside = zero_density_outside,
    recon_mode = recon_mode
  )
  
  class(out) <- c("smooth_trajectory_forecast", "new_trajectory_forecast", class(out))
  
  out
}
