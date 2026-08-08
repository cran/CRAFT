#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(file.path(getwd(), "CRAFT"), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "inst", "benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_craft <- function(root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
  } else {
    r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
    invisible(lapply(r_files, source))
  }
}

load_craft(root)

availability <- data.frame(
  item = c("CRAFT", "changepoint", "randomForest", "vars", "rugarch", "rmgarch", "copula"),
  available = c(
    exists("craft_fit", mode = "function"),
    requireNamespace("changepoint", quietly = TRUE),
    requireNamespace("randomForest", quietly = TRUE),
    requireNamespace("vars", quietly = TRUE),
    requireNamespace("rugarch", quietly = TRUE),
    requireNamespace("rmgarch", quietly = TRUE),
    requireNamespace("copula", quietly = TRUE)
  ),
  role = c(
    "model under test",
    "CRAFT regime segmentation dependency",
    "optional tree ensemble baseline",
    "package VAR baseline, unavailable here if FALSE",
    "GARCH baseline, unavailable here if FALSE",
    "DCC-GARCH baseline, unavailable here if FALSE",
    "copula scenario baseline, unavailable here if FALSE"
  )
)

write.csv(availability, file.path(out_dir, "availability.csv"), row.names = FALSE)

simulate_series <- function(n = 240, m = 3, seed = 20260802) {
  set.seed(seed)
  mus <- rbind(
    c(0.0007, 0.0004, 0.0002),
    c(-0.0012, -0.0008, 0.0003),
    c(0.0010, 0.0012, -0.0002)
  )
  vols <- rbind(
    c(0.009, 0.011, 0.008),
    c(0.021, 0.025, 0.016),
    c(0.013, 0.010, 0.019)
  )
  cors <- list(
    matrix(c(1.00, 0.35, 0.15, 0.35, 1.00, 0.25, 0.15, 0.25, 1.00), 3, 3),
    matrix(c(1.00, 0.78, 0.55, 0.78, 1.00, 0.62, 0.55, 0.62, 1.00), 3, 3),
    matrix(c(1.00, -0.10, 0.20, -0.10, 1.00, -0.15, 0.20, -0.15, 1.00), 3, 3)
  )
  regime <- rep(rep(1:3, times = c(35, 25, 30)), length.out = n)
  regime <- regime[c(seq(1, n, by = 2), seq(2, n, by = 2))]
  regime <- regime[seq_len(n)]
  ret <- matrix(NA_real_, n, m)
  for (i in seq_len(n)) {
    g <- regime[i]
    L <- chol(diag(vols[g, ], m) %*% cors[[g]] %*% diag(vols[g, ], m))
    z <- as.numeric(stats::rnorm(m) %*% L)
    if (i > 1) z <- z + 0.18 * ret[i - 1, ]
    ret[i, ] <- mus[g, ] + z
  }
  ret <- pmax(ret, -0.22)
  levels <- apply(1 + ret, 2, cumprod)
  levels <- as.data.frame(levels)
  names(levels) <- paste0("asset_", seq_len(m))
  list(levels = levels, returns = ret, regime = regime)
}

profile_matrix <- function(series, window, horizon) {
  emb <- trajectory_embedding(series, window)
  future_cols <- unlist(lapply(names(series), function(a) paste0(a, "_cum_lead_", horizon)), use.names = FALSE)
  list(
    B = as.matrix(emb$past_trajectories),
    Y = as.matrix(emb$future_trajectories[, future_cols, drop = FALSE]),
    latest = as.numeric(emb$latest_trajectory),
    future_cols = future_cols
  )
}

standardize_with <- function(x, center, scale) {
  scale[!is.finite(scale) | scale <= 0] <- 1
  sweep(sweep(as.matrix(x), 2, center, "-"), 2, scale, "/")
}

distance_weights <- function(B_train, b_now, bandwidth = NULL, k = 35) {
  center <- colMeans(B_train, na.rm = TRUE)
  scale <- apply(B_train, 2, stats::sd, na.rm = TRUE)
  Z <- standardize_with(B_train, center, scale)
  z <- as.numeric(standardize_with(matrix(b_now, nrow = 1), center, scale))
  d <- sqrt(rowSums((sweep(Z, 2, z, "-"))^2))
  k <- min(length(d), max(5L, k))
  keep <- order(d, decreasing = FALSE)[seq_len(k)]
  dd <- d[keep]
  if (is.null(bandwidth)) {
    bandwidth <- stats::median(dd[is.finite(dd) & dd > 0], na.rm = TRUE)
  }
  if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
  w <- exp(-0.5 * (dd / bandwidth)^2)
  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) w[] <- 1
  list(index = keep, weight = w / sum(w), distance = dd)
}

draw_rows <- function(Y, idx, w, n_draws) {
  take <- sample(seq_along(idx), size = n_draws, replace = TRUE, prob = w)
  Y[idx[take], , drop = FALSE]
}

draw_unconditional <- function(Y_train, n_draws) {
  Y_train[sample(seq_len(nrow(Y_train)), n_draws, replace = TRUE), , drop = FALSE]
}

draw_direct_analog <- function(B_train, b_now, Y_train, n_draws) {
  dw <- distance_weights(B_train, b_now, k = 35)
  draw_rows(Y_train, dw$index, dw$weight, n_draws)
}

draw_svd_analog <- function(B_train, b_now, Y_train, n_draws, n_pc = 3) {
  center <- colMeans(B_train, na.rm = TRUE)
  scale <- apply(B_train, 2, stats::sd, na.rm = TRUE)
  Z <- standardize_with(B_train, center, scale)
  pca <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pca$x), nrow(pca$rotation))
  score <- pca$x[, seq_len(n_pc), drop = FALSE]
  z_now <- as.numeric(standardize_with(matrix(b_now, nrow = 1), center, scale) %*% pca$rotation[, seq_len(n_pc), drop = FALSE])
  dw <- distance_weights(score, z_now, k = 35)
  draw_rows(Y_train, dw$index, dw$weight, n_draws)
}

draw_var_bootstrap <- function(series, horizon, n_draws, p = 2) {
  x <- as.matrix(series)
  lr <- diff(log(pmax(x, .Machine$double.eps)))
  n <- nrow(lr)
  m <- ncol(lr)
  p <- min(p, n - 5L)
  if (p < 1L) stop("not enough observations for VAR")
  Y <- lr[(p + 1L):n, , drop = FALSE]
  X <- do.call(cbind, lapply(seq_len(p), function(j) lr[(p + 1L - j):(n - j), , drop = FALSE]))
  X <- cbind(1, X)
  coef <- qr.solve(X, Y)
  resid <- Y - X %*% coef
  last_lags <- lr[(n - p + 1L):n, , drop = FALSE]
  out <- matrix(NA_real_, n_draws, m)
  for (i in seq_len(n_draws)) {
    lags <- last_lags
    sim <- matrix(NA_real_, horizon, m)
    for (h in seq_len(horizon)) {
      lag_vec <- as.numeric(t(lags[p:1, , drop = FALSE]))
      eps <- resid[sample(seq_len(nrow(resid)), 1L), ]
      next_ret <- as.numeric(c(1, lag_vec) %*% coef + eps)
      sim[h, ] <- next_ret
      lags <- rbind(lags[-1L, , drop = FALSE], next_ret)
    }
    out[i, ] <- exp(colSums(sim)) - 1
  }
  colnames(out) <- colnames(series)
  out
}

draw_random_forest <- function(B_train, b_now, Y_train, n_draws, seed) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest package is not installed")
  }
  set.seed(seed)
  X <- as.data.frame(B_train)
  names(X) <- paste0("x", seq_len(ncol(X)))
  newx <- as.data.frame(matrix(b_now, nrow = 1))
  names(newx) <- names(X)
  pred <- matrix(NA_real_, n_draws, ncol(Y_train))
  for (j in seq_len(ncol(Y_train))) {
    dat <- X
    dat$y <- Y_train[, j]
    fit <- randomForest::randomForest(y ~ ., data = dat, ntree = n_draws)
    all_tree <- stats::predict(fit, newdata = newx, predict.all = TRUE)$individual
    pred[, j] <- as.numeric(all_tree[1, ])
  }
  colnames(pred) <- colnames(Y_train)
  pred
}

extract_craft_horizon <- function(fit, assets, horizon) {
  draws <- fit$trajectory_draws
  cols <- paste0(assets, "_cum_lead_", horizon)
  if (all(cols %in% names(draws))) return(as.matrix(draws[, cols, drop = FALSE]))
  cols <- grep(paste0("_cum_lead_", horizon, "$"), names(draws), value = TRUE)
  if (length(cols) >= length(assets)) return(as.matrix(draws[, cols[seq_along(assets)], drop = FALSE]))
  stop("Could not find requested horizon columns in CRAFT trajectory draws")
}

crps_one <- function(x, y, max_pairs = 1200) {
  x <- as.numeric(x)
  n <- length(x)
  a <- mean(abs(x - y))
  if (n > max_pairs) x2 <- sample(x, max_pairs) else x2 <- x
  b <- mean(abs(sample(x2, length(x2), replace = TRUE) - sample(x2, length(x2), replace = TRUE)))
  a - 0.5 * b
}

pinball_loss <- function(draws, y, probs = c(0.05, 0.10, 0.50, 0.90, 0.95)) {
  q <- apply(draws, 2, stats::quantile, probs = probs, na.rm = TRUE)
  vals <- numeric(0)
  for (a in seq_len(ncol(draws))) {
    for (i in seq_along(probs)) {
      u <- y[a] - q[i, a]
      vals <- c(vals, (probs[i] - as.numeric(u < 0)) * u)
    }
  }
  mean(vals)
}

energy_score <- function(draws, y, max_pairs = 1200) {
  D <- as.matrix(draws)
  n <- nrow(D)
  a <- mean(sqrt(rowSums((sweep(D, 2, y, "-"))^2)))
  if (n > max_pairs) D <- D[sample(seq_len(n), max_pairs), , drop = FALSE]
  i <- sample(seq_len(nrow(D)), nrow(D), replace = TRUE)
  j <- sample(seq_len(nrow(D)), nrow(D), replace = TRUE)
  b <- mean(sqrt(rowSums((D[i, , drop = FALSE] - D[j, , drop = FALSE])^2)))
  a - 0.5 * b
}

variogram_score <- function(draws, y, p = 0.5) {
  D <- as.matrix(draws)
  m <- ncol(D)
  out <- 0
  for (i in seq_len(m - 1L)) {
    for (j in (i + 1L):m) {
      obs <- abs(y[i] - y[j])^p
      pred <- mean(abs(D[, i] - D[, j])^p)
      out <- out + (obs - pred)^2
    }
  }
  out
}

score_draws <- function(draws, y, train_y, realized_pair_threshold) {
  D <- as.matrix(draws)
  y <- as.numeric(y)
  crps <- mean(vapply(seq_len(ncol(D)), function(j) crps_one(D[, j], y[j]), numeric(1)))
  med <- apply(D, 2, stats::median, na.rm = TRUE)
  lo <- apply(D, 2, stats::quantile, probs = 0.05, na.rm = TRUE)
  hi <- apply(D, 2, stats::quantile, probs = 0.95, na.rm = TRUE)
  port <- rowMeans(D)
  port_y <- mean(y)
  var05 <- as.numeric(stats::quantile(port, 0.05, na.rm = TRUE))
  es05 <- mean(port[port <= var05])
  pred_co <- mean(D[, 1] <= realized_pair_threshold[1] & D[, 2] <= realized_pair_threshold[2])
  obs_co <- as.numeric(y[1] <= realized_pair_threshold[1] & y[2] <= realized_pair_threshold[2])
  data.frame(
    crps = crps,
    pinball = pinball_loss(D, y),
    directional_accuracy = mean(sign(med) == sign(y)),
    interval90_coverage = mean(y >= lo & y <= hi),
    interval90_width = mean(hi - lo),
    energy = energy_score(D, y),
    variogram = variogram_score(D, y),
    portfolio_var05_breach = as.numeric(port_y < var05),
    portfolio_es05_gap = ifelse(port_y < var05, port_y - es05, 0),
    tail_coexceed_brier = (pred_co - obs_co)^2,
    stringsAsFactors = FALSE
  )
}

effective_n <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  p <- w / sum(w)
  1 / sum(p^2)
}

first_or_na <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  x[1L]
}

key_rows <- function(x) {
  if (is.null(x) || !nrow(x)) return(character(0))
  apply(as.data.frame(x, stringsAsFactors = FALSE), 1, paste, collapse = "|")
}

hac_test <- function(diff, lag = 0L) {
  diff <- diff[is.finite(diff)]
  n <- length(diff)
  if (n < 4L) {
    return(c(n = n, mean_diff = mean(diff), se = NA_real_, t = NA_real_, p_value = NA_real_))
  }
  x <- diff - mean(diff)
  gamma0 <- sum(x * x) / n
  v <- gamma0
  max_lag <- min(as.integer(lag), n - 2L)
  if (max_lag > 0L) {
    for (l in seq_len(max_lag)) {
      g <- sum(x[(l + 1L):n] * x[seq_len(n - l)]) / n
      v <- v + 2 * (1 - l / (max_lag + 1L)) * g
    }
  }
  se <- sqrt(max(v, 0) / n)
  tval <- mean(diff) / se
  p <- 2 * stats::pt(-abs(tval), df = n - 1L)
  c(n = n, mean_diff = mean(diff), se = se, t = tval, p_value = p)
}

set.seed(20260802)
cfg <- list(n = 240, m = 3, window = 8, horizon = 5, n_draws = 250, n_origins = 18)
sim <- simulate_series(n = cfg$n, m = cfg$m)
series <- sim$levels
assets <- names(series)
origins <- unique(round(seq(120, cfg$n - cfg$horizon - 1L, length.out = cfg$n_origins)))

metrics <- list()
draw_store <- list()
diagnostics <- list()

models <- c("CRAFT", "direct_analog", "svd_analog", "unconditional_bootstrap", "ols_var_bootstrap")
if (requireNamespace("randomForest", quietly = TRUE)) models <- c(models, "random_forest")

for (origin in origins) {
  train <- series[seq_len(origin), , drop = FALSE]
  y <- as.numeric(series[origin + cfg$horizon, ] / series[origin, ] - 1)
  names(y) <- assets
  prof <- profile_matrix(train, cfg$window, cfg$horizon)
  B_train <- prof$B
  Y_train <- prof$Y
  b_now <- prof$latest
  colnames(Y_train) <- assets
  tail_thr <- apply(Y_train[, 1:2, drop = FALSE], 2, stats::quantile, probs = 0.10, na.rm = TRUE)

  origin_draws <- list()

  fit <- tryCatch(
    craft_fit(
      train,
      window = cfg$window,
      horizon = cfg$horizon,
      n_draws = cfg$n_draws,
      n_factors = 1,
      min_segment = 8,
      max_regimes_per_factor = 3,
      n_testing = 0,
      valid_frac = 0.20,
      valid_split = "last",
      tail_penalty_lambda = 0,
      sampler_weight = "uniform",
      return_train_probs = FALSE,
      verbose = FALSE,
      seed = 1000 + origin
    ),
    error = function(e) e
  )
  if (!inherits(fit, "error")) {
    origin_draws$CRAFT <- extract_craft_horizon(fit, assets, cfg$horizon)
    past_keys <- key_rows(fit$past_labels)
    latest_key <- key_rows(fit$latest_regime_labels)
    key_counts <- table(past_keys)
    sampler_prob <- fit$sampler_trace$probability
    sampler_prob <- sampler_prob[is.finite(sampler_prob) & sampler_prob > 0]
    sampler_eff_n <- if (length(sampler_prob)) {
      prob <- sampler_prob / sum(sampler_prob)
      1 / sum(prob^2)
    } else {
      NA_real_
    }
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      origin = origin,
      model = "CRAFT",
      n_backward_regime_keys = length(key_counts),
      min_observations_per_key = ifelse(length(key_counts), min(key_counts), NA_integer_),
      median_observations_per_key = ifelse(length(key_counts), stats::median(key_counts), NA_real_),
      max_observations_per_key = ifelse(length(key_counts), max(key_counts), NA_integer_),
      exact_current_key_observations = ifelse(length(latest_key), sum(past_keys == latest_key), NA_integer_),
      sampled_unique_paths = length(unique(fit$sampler_trace$index)),
      sampler_probability_effective_n = sampler_eff_n,
      sampler_weight_effective_n = first_or_na(fit$diagnostics$sampler_weight_effective_n),
      posterior_entropy_mean = first_or_na(mean(fit$diagnostics$posterior_entropy, na.rm = TRUE)),
      latest_unseen_total = first_or_na(fit$diagnostics$latest_unseen_total),
      stringsAsFactors = FALSE
    )
  }

  set.seed(2000 + origin)
  origin_draws$direct_analog <- draw_direct_analog(B_train, b_now, Y_train, cfg$n_draws)
  set.seed(3000 + origin)
  origin_draws$svd_analog <- draw_svd_analog(B_train, b_now, Y_train, cfg$n_draws)
  set.seed(4000 + origin)
  origin_draws$unconditional_bootstrap <- draw_unconditional(Y_train, cfg$n_draws)
  set.seed(5000 + origin)
  origin_draws$ols_var_bootstrap <- draw_var_bootstrap(train, cfg$horizon, cfg$n_draws, p = 2)
  if ("random_forest" %in% models) {
    origin_draws$random_forest <- draw_random_forest(B_train, b_now, Y_train, cfg$n_draws, seed = 6000 + origin)
  }

  for (model in names(origin_draws)) {
    sc <- score_draws(origin_draws[[model]], y, Y_train, tail_thr)
    metrics[[length(metrics) + 1L]] <- cbind(origin = origin, model = model, sc)
    draw_store[[paste(origin, model, sep = "::")]] <- list(draws = origin_draws[[model]], y = y)
  }
}

metric_df <- do.call(rbind, metrics)
row.names(metric_df) <- NULL
write.csv(metric_df, file.path(out_dir, "forecast_metrics_by_origin.csv"), row.names = FALSE)

mean_by_model <- aggregate(
  metric_df[, setdiff(names(metric_df), c("origin", "model")), drop = FALSE],
  by = list(model = metric_df$model),
  FUN = mean,
  na.rm = TRUE
)

cor_error <- vapply(mean_by_model$model, function(model) {
  pred_cors <- list()
  ys <- list()
  for (origin in origins) {
    key <- paste(origin, model, sep = "::")
    if (!is.null(draw_store[[key]])) {
      pred_cors[[length(pred_cors) + 1L]] <- stats::cor(draw_store[[key]]$draws)
      ys[[length(ys) + 1L]] <- draw_store[[key]]$y
    }
  }
  if (length(pred_cors) < 2L) return(NA_real_)
  pred_cor <- Reduce("+", pred_cors) / length(pred_cors)
  ymat <- do.call(rbind, ys)
  obs_cor <- stats::cor(ymat)
  mean(abs(pred_cor[upper.tri(pred_cor)] - obs_cor[upper.tri(obs_cor)]))
}, numeric(1))

mean_by_model$cross_asset_correlation_error <- cor_error[mean_by_model$model]
write.csv(mean_by_model, file.path(out_dir, "model_summary.csv"), row.names = FALSE)

loss_metrics <- c("crps", "pinball", "energy", "variogram", "tail_coexceed_brier")
tests <- list()
base <- "direct_analog"
for (model in setdiff(unique(metric_df$model), base)) {
  for (metric in loss_metrics) {
    wide <- merge(
      metric_df[metric_df$model == base, c("origin", metric)],
      metric_df[metric_df$model == model, c("origin", metric)],
      by = "origin",
      suffixes = c("_base", "_model")
    )
    diff <- wide[[paste0(metric, "_base")]] - wide[[paste0(metric, "_model")]]
    h <- hac_test(diff, lag = cfg$horizon - 1L)
    tests[[length(tests) + 1L]] <- data.frame(
      baseline = base,
      model = model,
      metric = metric,
      n = unname(h["n"]),
      mean_loss_improvement = unname(h["mean_diff"]),
      hac_se = unname(h["se"]),
      t = unname(h["t"]),
      p_value = unname(h["p_value"]),
      interpretation = ifelse(unname(h["mean_diff"]) > 0, "model lower loss than direct analog", "model higher loss than direct analog"),
      stringsAsFactors = FALSE
    )
  }
}

test_df <- do.call(rbind, tests)
write.csv(test_df, file.path(out_dir, "comparison_tests.csv"), row.names = FALSE)

if (length(diagnostics)) {
  diag_df <- do.call(rbind, diagnostics)
  write.csv(diag_df, file.path(out_dir, "craft_regime_support_diagnostics.csv"), row.names = FALSE)
}

config_df <- data.frame(
  key = names(cfg),
  value = unlist(cfg),
  row.names = NULL
)
write.csv(config_df, file.path(out_dir, "benchmark_config.csv"), row.names = FALSE)

cat("\nModel summary\n")
print(mean_by_model[order(mean_by_model$crps), ], row.names = FALSE)
cat("\nComparison tests versus direct_analog; positive mean_loss_improvement favors the listed model\n")
print(test_df[order(test_df$metric, -test_df$mean_loss_improvement), ], row.names = FALSE)
cat("\nWrote results to:", out_dir, "\n")

