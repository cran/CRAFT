# ============================================================
# Helper: align named probability lists
# ============================================================

.align_named_probs <- function(x, target_names, object_name = "object") {
  if (!is.list(x)) {
    stop("`", object_name, "` must be a named list.")
  }

  if (is.null(names(x))) {
    stop("`", object_name, "` must have names.")
  }

  missing_names <- setdiff(target_names, names(x))
  extra_names <- setdiff(names(x), target_names)

  if (length(missing_names) > 0L || length(extra_names) > 0L) {
    stop(
      "`", object_name, "` names do not match target regime names.\n",
      "Missing: ", paste(missing_names, collapse = ", "), "\n",
      "Extra: ", paste(extra_names, collapse = ", ")
    )
  }

  x[target_names]
}


# ============================================================
# Helper: sanitize weights for regime_sampler()
# ============================================================

.sanitize_regime_weights <- function(w, target_names) {
  w <- as.numeric(w)

  if (length(w) != length(target_names)) {
    stop("Regime weight vector has wrong length.")
  }

  names(w) <- target_names

  w[!is.finite(w)] <- NA_real_

  if (all(is.na(w))) {
    w[] <- 1
  } else {
    w[is.na(w)] <- mean(w, na.rm = TRUE)
  }

  # Avoid fully zero or negative pathological weights
  w[w < 0] <- 0

  if (sum(w) <= 0) {
    w[] <- 1
  }

  w
}


# ============================================================
# Helper: posterior entropy diagnostics
# ============================================================

.prob_entropy01 <- function(prob_mat) {
  p <- as.matrix(prob_mat)
  p <- pmax(p, 1e-12)
  k <- ncol(p)

  if (k <= 1L) {
    return(rep(0, nrow(p)))
  }

  -rowSums(p * log(p)) / log(k)
}


.effective_n <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0

  if (sum(w) <= 0) {
    return(NA_real_)
  }

  p <- w / sum(w)
  1 / sum(p^2)
}


.align_trajectory_columns <- function(trajectory_draws, asset_names, window) {
  trajectory_draws <- as.data.frame(trajectory_draws)

  expected <- unlist(
    lapply(asset_names, function(asset) {
      paste0(asset, "_cum_lead_", seq_len(window))
    }),
    use.names = FALSE
  )

  used_names <- FALSE

  if (!is.null(colnames(trajectory_draws)) && all(expected %in% colnames(trajectory_draws))) {
    trajectory_draws <- trajectory_draws[, expected, drop = FALSE]
    used_names <- TRUE
  } else {
    trajectory_size <- length(asset_names) * window

    if (ncol(trajectory_draws) < trajectory_size) {
      stop(
        "`trajectory_draws` has fewer columns than expected.\n",
        "Expected at least: ", trajectory_size,
        "\nActual: ", ncol(trajectory_draws)
      )
    }

    trajectory_draws <- trajectory_draws[, seq_len(trajectory_size), drop = FALSE]
    colnames(trajectory_draws) <- expected
  }

  list(
    trajectory_draws = trajectory_draws,
    used_names = used_names,
    expected_names = expected
  )
}


# ============================================================
# Helper: soft tail penalty for future-trajectory sampling
# ============================================================

.v3_normalize_prob <- function(p, eps = 1e-12) {
  p <- as.numeric(p)

  if (!length(p)) {
    return(p)
  }

  p[!is.finite(p) | p < 0] <- 0

  if (sum(p) <= eps) {
    p[] <- 1 / length(p)
  } else {
    p <- p / sum(p)
  }

  p
}


.v3_softmax <- function(x) {
  x <- as.numeric(x)

  if (!length(x)) {
    return(x)
  }

  x[!is.finite(x)] <- -Inf

  if (all(!is.finite(x))) {
    return(rep(1 / length(x), length(x)))
  }

  z <- exp(x - max(x, na.rm = TRUE))

  if (sum(z) <= 0 || any(!is.finite(z))) {
    return(rep(1 / length(x), length(x)))
  }

  z / sum(z)
}


.numeric_summary <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(c(
      min = NA_real_,
      q25 = NA_real_,
      median = NA_real_,
      mean = NA_real_,
      q75 = NA_real_,
      max = NA_real_
    ))
  }

  c(
    min = min(x),
    q25 = unname(stats::quantile(x, 0.25, names = FALSE)),
    median = stats::median(x),
    mean = mean(x),
    q75 = unname(stats::quantile(x, 0.75, names = FALSE)),
    max = max(x)
  )
}


.tail_penalty_scores <- function(
  X,
  threshold = 3,
  power = 2,
  aggregation = c("mean", "max", "sum")
) {
  aggregation <- match.arg(aggregation)

  X <- as.data.frame(X, stringsAsFactors = FALSE)
  X[] <- lapply(X, as.numeric)

  if (!nrow(X) || !ncol(X)) {
    return(list(
      scores = rep(0, nrow(X)),
      columns = names(X),
      center = numeric(0),
      scale = numeric(0),
      summary = .numeric_summary(rep(0, nrow(X)))
    ))
  }

  mat <- as.matrix(X)

  center <- apply(mat, 2, stats::median, na.rm = TRUE)
  center[!is.finite(center)] <- 0

  scale <- apply(mat, 2, stats::mad, na.rm = TRUE)
  fallback_scale <- apply(mat, 2, stats::sd, na.rm = TRUE)
  bad_scale <- !is.finite(scale) | scale <= 0
  scale[bad_scale] <- fallback_scale[bad_scale]
  scale[!is.finite(scale) | scale <= 0] <- 1

  z <- abs(sweep(sweep(mat, 2, center, "-"), 2, scale, "/"))
  excess <- pmax(z - threshold, 0)^power
  excess[!is.finite(excess)] <- 0

  scores <- switch(
    aggregation,
    mean = rowMeans(excess, na.rm = TRUE),
    max = apply(excess, 1, max, na.rm = TRUE),
    sum = rowSums(excess, na.rm = TRUE)
  )

  scores[!is.finite(scores)] <- 0

  list(
    scores = as.numeric(scores),
    columns = names(X),
    center = center,
    scale = scale,
    summary = .numeric_summary(scores)
  )
}


.regime_sampler_tail_penalized <- function(
  trafo_df,
  cluster_df,
  prob_vectors,
  weights = NULL,
  tail_penalty = NULL,
  tail_penalty_lambda = 0,
  size = 1000L,
  seed = NULL,
  eps = 1e-12
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  trafo_df <- as.data.frame(trafo_df, stringsAsFactors = FALSE)
  cluster_df <- as.data.frame(cluster_df, stringsAsFactors = FALSE)

  if (nrow(trafo_df) != nrow(cluster_df)) {
    stop("`trafo_df` and `cluster_df` must have the same number of rows.")
  }

  if (!nrow(trafo_df)) {
    stop("Cannot sample from an empty `trafo_df`.")
  }

  target_names <- names(cluster_df)

  if (is.null(target_names) || any(!nzchar(target_names))) {
    target_names <- paste0("target_", seq_len(ncol(cluster_df)))
    names(cluster_df) <- target_names
  }

  if (is.null(weights)) {
    weights <- rep(1, length(target_names))
  }

  weights <- as.numeric(weights)

  if (length(weights) != length(target_names)) {
    weights <- rep(1, length(target_names))
  }

  weights[!is.finite(weights) | weights < 0] <- 0

  if (sum(weights) <= 0) {
    weights[] <- 1
  }

  weights <- weights / sum(weights)
  names(weights) <- target_names

  base_log_score <- rep(0, nrow(cluster_df))

  for (target in target_names) {
    pv_obj <- prob_vectors[[target]]

    if (is.null(pv_obj)) {
      next
    }

    if (is.data.frame(pv_obj) || is.matrix(pv_obj)) {
      pv <- pv_obj[1, , drop = TRUE]
    } else {
      pv <- pv_obj
    }

    pv <- as.numeric(pv)
    nm <- colnames(pv_obj)

    if (is.null(nm)) {
      nm <- names(pv_obj)
    }

    if (is.null(nm) && !is.null(colnames(as.matrix(pv_obj)))) {
      nm <- colnames(as.matrix(pv_obj))
    }

    if (is.null(nm) || length(nm) != length(pv)) {
      nm <- as.character(sort(unique(cluster_df[[target]])))
    }

    names(pv) <- nm
    pv <- .v3_normalize_prob(pv, eps = eps)

    states <- as.character(cluster_df[[target]])
    p_row <- pv[states]
    p_row[!is.finite(p_row)] <- eps
    base_log_score <- base_log_score + weights[target] * log(pmax(p_row, eps))
  }

  if (is.null(tail_penalty)) {
    tail_penalty <- rep(0, nrow(cluster_df))
  }

  tail_penalty <- as.numeric(tail_penalty)

  if (length(tail_penalty) != nrow(cluster_df)) {
    stop("`tail_penalty` must have one value per candidate row.")
  }

  tail_penalty[!is.finite(tail_penalty) | tail_penalty < 0] <- 0

  final_log_score <- base_log_score - tail_penalty_lambda * tail_penalty
  prob <- .v3_softmax(final_log_score)

  idx <- sample.int(
    nrow(trafo_df),
    size = as.integer(size),
    replace = TRUE,
    prob = prob
  )

  draws <- trafo_df[idx, , drop = FALSE]
  rownames(draws) <- NULL

  list(
    draws = draws,
    index = idx,
    sampling_prob = prob,
    base_log_score = base_log_score,
    final_log_score = final_log_score,
    tail_penalty = tail_penalty
  )
}


# ============================================================
# Helper: transition-regime prediction engine
# ============================================================

.transition_row_weights <- function(row_weights, n) {
  if (is.null(row_weights)) return(rep(1, n))
  row_weights <- as.numeric(row_weights)
  if (length(row_weights) != n) stop("`row_weights` must have length equal to nrow(X).")
  if (any(!is.finite(row_weights)) || any(row_weights < 0)) stop("`row_weights` must contain finite non-negative values.")
  if (sum(row_weights) <= 0) stop("`row_weights` must contain at least one positive value.")
  row_weights / mean(row_weights)
}

.transition_row_key <- function(X, sep = "\r") {
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  if (!ncol(X)) return(rep("__NO_FEATURES__", nrow(X)))
  parts <- lapply(X, function(z) {
    z <- as.character(z)
    z[is.na(z)] <- "__NA__"
    z
  })
  do.call(paste, c(parts, sep = sep))
}

.transition_prob_to_class <- function(prob_df) {
  prob_mat <- as.matrix(prob_df)
  colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
}

.transition_make_split <- function(n, valid_frac, valid_split, seed = NULL) {
  if (n < 2L) return(list(train_idx = seq_len(n), valid_idx = seq_len(n)))
  valid_frac <- as.numeric(valid_frac)[1L]
  if (!is.finite(valid_frac) || valid_frac <= 0) valid_frac <- 0.2
  n_valid <- max(1L, floor(n * valid_frac))
  n_valid <- min(n - 1L, n_valid)
  if (!is.null(seed)) set.seed(seed)
  if (valid_split == "last") {
    valid_idx <- seq.int(n - n_valid + 1L, n)
  } else if (valid_split == "random") {
    valid_idx <- sort(sample.int(n, n_valid))
  } else if (valid_split == "blocked") {
    start <- sample.int(n - n_valid + 1L, 1L)
    valid_idx <- seq.int(start, start + n_valid - 1L)
  } else {
    stop("Unknown `valid_split`.")
  }
  list(train_idx = setdiff(seq_len(n), valid_idx), valid_idx = valid_idx)
}

.transition_classification_metrics <- function(Y_valid, pred_df, target_levels, row_weights = NULL) {
  Y_valid <- as.data.frame(Y_valid, stringsAsFactors = FALSE)
  pred_df <- as.data.frame(pred_df, stringsAsFactors = FALSE)
  if (is.null(row_weights)) row_weights <- rep(1, nrow(Y_valid))
  row_weights <- as.numeric(row_weights)
  row_weights[!is.finite(row_weights) | row_weights < 0] <- 0
  out <- lapply(names(target_levels), function(target) {
    lev <- target_levels[[target]]
    truth <- factor(as.character(Y_valid[[target]]), levels = lev)
    pred <- factor(as.character(pred_df[[target]]), levels = lev)
    ok <- !is.na(truth) & !is.na(pred)
    correct <- pred[ok] == truth[ok]
    w <- row_weights[ok]
    acc <- if (length(correct)) mean(correct) else NA_real_
    weighted_acc <- if (length(correct) && sum(w) > 0) sum(w * as.numeric(correct)) / sum(w) else NA_real_
    present <- lev[lev %in% as.character(truth[ok])]
    recalls <- vapply(present, function(cl) {
      idx <- ok & truth == cl
      if (!any(idx)) return(NA_real_)
      mean(pred[idx] == truth[idx])
    }, numeric(1))
    weighted_recalls <- vapply(present, function(cl) {
      idx <- ok & truth == cl
      ww <- row_weights[idx]
      if (!any(idx) || sum(ww) <= 0) return(NA_real_)
      sum(ww * as.numeric(pred[idx] == truth[idx])) / sum(ww)
    }, numeric(1))
    data.frame(target = target, accuracy = acc, balanced_accuracy = mean(recalls, na.rm = TRUE), weighted_accuracy = weighted_acc, weighted_balanced_accuracy = mean(weighted_recalls, na.rm = TRUE), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.transition_valid_joint_accuracy <- function(Y_valid, pred_df) {
  Y_valid <- as.data.frame(Y_valid, stringsAsFactors = FALSE)
  pred_df <- as.data.frame(pred_df, stringsAsFactors = FALSE)
  if (!nrow(Y_valid)) return(NA_real_)
  row_ok <- rep(TRUE, nrow(Y_valid))
  for (target in names(Y_valid)) row_ok <- row_ok & as.character(pred_df[[target]]) == as.character(Y_valid[[target]])
  mean(row_ok, na.rm = TRUE)
}

.fit_transition_target <- function(X, y, levels_y, row_weights, alpha = 1, min_count = 1, key_sep = "\r") {
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  y <- factor(as.character(y), levels = levels_y)
  row_weights <- .transition_row_weights(row_weights, nrow(X))
  global_counts <- setNames(rep(0, length(levels_y)), levels_y)
  for (cl in levels_y) {
    idx <- !is.na(y) & y == cl
    global_counts[cl] <- sum(row_weights[idx])
  }
  global_prob <- (global_counts + alpha) / (sum(global_counts) + alpha * length(levels_y))
  keys <- .transition_row_key(X, sep = key_sep)
  key_levels <- sort(unique(keys[!is.na(y)]))
  counts <- vector("list", length(key_levels))
  names(counts) <- key_levels
  for (key in key_levels) {
    idx_key <- keys == key & !is.na(y)
    key_counts <- setNames(rep(0, length(levels_y)), levels_y)
    for (cl in levels_y) {
      idx <- idx_key & y == cl
      key_counts[cl] <- sum(row_weights[idx])
    }
    counts[[key]] <- key_counts
  }
  predict_prob <- function(newdata) {
    newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
    new_keys <- .transition_row_key(newdata, sep = key_sep)
    mat <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(levels_y))
    colnames(mat) <- levels_y
    for (i in seq_len(nrow(newdata))) {
      key_counts <- counts[[new_keys[i]]]
      key_total <- if (is.null(key_counts)) 0 else sum(key_counts)
      if (is.null(key_counts) || key_total < min_count) {
        prob <- global_prob
      } else {
        prob <- (key_counts + alpha * global_prob) / (key_total + alpha)
      }
      mat[i, ] <- prob[levels_y]
    }
    as.data.frame(mat, stringsAsFactors = FALSE)
  }
  list(engine = "transition", levels = levels_y, alpha = alpha, min_count = min_count, key_sep = key_sep, global_prob = global_prob, counts = counts, predict_prob = predict_prob)
}

.fit_transition_classifier <- function(X, Y, row_weights = NULL, valid_frac = 0.2, valid_split = "last", seed = NULL, transition_alpha = 1, transition_min_count = 1, transition_key_sep = "\r") {
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  Y <- as.data.frame(Y, stringsAsFactors = FALSE)
  row_weights <- .transition_row_weights(row_weights, nrow(X))
  target_levels <- lapply(Y, function(z) levels(factor(z)))
  names(target_levels) <- names(Y)
  split <- .transition_make_split(nrow(X), valid_frac, valid_split, seed)
  final_models <- vector("list", ncol(Y)); names(final_models) <- names(Y)
  valid_pred_list <- setNames(vector("list", ncol(Y)), names(Y))
  for (target in names(Y)) {
    eval_model <- .fit_transition_target(X[split$train_idx, , drop = FALSE], Y[[target]][split$train_idx], target_levels[[target]], row_weights[split$train_idx], transition_alpha, transition_min_count, transition_key_sep)
    valid_prob <- eval_model$predict_prob(X[split$valid_idx, , drop = FALSE])
    valid_pred_list[[target]] <- .transition_prob_to_class(valid_prob)
    final_models[[target]] <- .fit_transition_target(X, Y[[target]], target_levels[[target]], row_weights, transition_alpha, transition_min_count, transition_key_sep)
  }
  valid_preds <- as.data.frame(valid_pred_list, stringsAsFactors = FALSE)
  metrics <- .transition_classification_metrics(Y[split$valid_idx, , drop = FALSE], valid_preds, target_levels, row_weights[split$valid_idx])
  valid_joint_acc <- .transition_valid_joint_accuracy(Y[split$valid_idx, names(Y), drop = FALSE], valid_preds)
  predict <- function(newdata, type = c("prob", "class")) {
    type <- match.arg(type)
    newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
    probs <- lapply(names(final_models), function(target) final_models[[target]]$predict_prob(newdata))
    names(probs) <- names(final_models)
    if (type == "prob") return(probs)
    as.data.frame(setNames(lapply(probs, .transition_prob_to_class), names(probs)), stringsAsFactors = FALSE)
  }
  out <- list(engine = "transition", models = final_models, target_levels = target_levels, predict = predict, diagnostics = list(metrics = metrics, valid_joint_acc = valid_joint_acc, valid_split = valid_split, train_idx = split$train_idx, valid_idx = split$valid_idx, temperature = NA_real_, transition = list(alpha = transition_alpha, min_count = transition_min_count, key_sep = transition_key_sep)))
  class(out) <- c("regime_classifier_transition", "regime_classifier")
  out
}

# ============================================================
# Helper: internal equal-distance backtesting
# ============================================================

.v5_backtest_horizons <- function(window, horizon) {
  if (is.null(horizon)) seq_len(window) else as.integer(horizon)
}

.v5_select_backtest_origins <- function(n_obs, n_testing, window, horizon, min_segment, backtest_min_train = NULL) {
  n_testing <- as.integer(n_testing)[1L]
  if (!is.finite(n_testing) || n_testing < 1L) return(integer(0))
  max_h <- max(.v5_backtest_horizons(window, horizon))
  last_origin <- n_obs - max_h
  min_train <- if (is.null(backtest_min_train)) max(window + 2L, window + 2L * as.integer(min_segment), 3L * window) else as.integer(backtest_min_train)[1L]
  min_origin <- max(window + 1L, min_train)
  if (!is.finite(min_origin) || min_origin > last_origin) return(integer(0))
  if (n_testing == 1L) return(as.integer(last_origin))
  origins <- unique(as.integer(round(seq(min_origin, last_origin, length.out = n_testing))))
  origins <- origins[origins >= min_origin & origins <= last_origin]
  sort(origins)
}

.v5_bind_rows <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(data.frame())
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.v5_backtest_one_origin <- function(task) {
  tryCatch({
    series <- as.data.frame(task$series, stringsAsFactors = FALSE)
    volumes <- task$volumes
    origin <- as.integer(task$origin)
    point_id <- as.integer(task$point_id)
    horizons <- task$horizons
    fit_args <- task$fit_args
    train_series <- series[seq_len(origin), , drop = FALSE]
    train_volumes <- NULL
    if (!is.null(volumes)) {
      volumes <- as.data.frame(volumes, stringsAsFactors = FALSE)
      train_volumes <- volumes[seq_len(origin), , drop = FALSE]
    }
    fit_args$series <- train_series
    fit_args$volumes <- train_volumes
    fit_args$seed <- as.integer(fit_args$seed) + point_id
    fit_args$n_testing <- 0L
    fit_args$backtest_parallel <- FALSE
    fit_args$backtest_workers <- 1L
    fit_args$return_train_probs <- FALSE
    fit_args$verbose <- FALSE
    fit <- do.call(craft_fit, fit_args)
    asset_names <- colnames(series)
    last_level <- as.numeric(series[origin, asset_names, drop = TRUE])
    names(last_level) <- asset_names
    rows <- list()
    row_id <- 0L
    for (asset in asset_names) {
      base_level <- last_level[[asset]]
      for (h in horizons) {
        col <- paste0(asset, "_cum_lead_", h)
        row_id <- row_id + 1L
        actual_level <- as.numeric(series[origin + h, asset])
        actual_return <- if (is.finite(base_level) && base_level != 0) actual_level / base_level - 1 else NA_real_
        draw_return <- if (col %in% names(fit$trajectory_draws)) as.numeric(fit$trajectory_draws[[col]]) else numeric(0)
        draw_return <- draw_return[is.finite(draw_return)]
        if (length(draw_return)) {
          q <- stats::quantile(draw_return, probs = c(0.05, 0.50, 0.95), names = FALSE, na.rm = TRUE)
          pred_return_mean <- mean(draw_return, na.rm = TRUE)
          level_draw <- base_level * (1 + draw_return)
          q_level <- stats::quantile(level_draw, probs = c(0.05, 0.50, 0.95), names = FALSE, na.rm = TRUE)
          pred_level_mean <- mean(level_draw, na.rm = TRUE)
          p_actual <- if (is.finite(actual_return)) mean(draw_return <= actual_return) else NA_real_
        } else {
          q <- rep(NA_real_, 3L)
          pred_return_mean <- NA_real_
          q_level <- rep(NA_real_, 3L)
          pred_level_mean <- NA_real_
          p_actual <- NA_real_
        }
        rows[[row_id]] <- data.frame(
          point_id = point_id,
          origin_index = origin,
          horizon = h,
          asset = asset,
          actual_return = actual_return,
          pred_return_mean = pred_return_mean,
          pred_return_q05 = q[1L],
          pred_return_median = q[2L],
          pred_return_q95 = q[3L],
          actual_level = actual_level,
          pred_level_mean = pred_level_mean,
          pred_level_q05 = q_level[1L],
          pred_level_median = q_level[2L],
          pred_level_q95 = q_level[3L],
          return_abs_error = abs(q[2L] - actual_return),
          return_sq_error = (q[2L] - actual_return)^2,
          level_abs_error = abs(q_level[2L] - actual_level),
          level_sq_error = (q_level[2L] - actual_level)^2,
          covered_90 = is.finite(actual_return) && actual_return >= q[1L] && actual_return <= q[3L],
          direction_correct = is.finite(actual_return) && is.finite(q[2L]) && sign(actual_return) == sign(q[2L]),
          p_actual = p_actual,
          fit_valid_acc = fit$valid_acc,
          fit_valid_bal_acc = fit$valid_bal_acc,
          fit_valid_joint_acc = fit$valid_joint_acc,
          stringsAsFactors = FALSE
        )
      }
    }
    list(ok = TRUE, point_id = point_id, origin_index = origin, predictions = .v5_bind_rows(rows), error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, point_id = task$point_id, origin_index = task$origin, predictions = data.frame(), error = conditionMessage(e))
  })
}

.v5_backtest_summary <- function(predictions) {
  if (!nrow(predictions)) {
    empty <- data.frame()
    return(list(overall = empty, by_asset = empty, by_horizon = empty, by_asset_horizon = empty))
  }
  summarize <- function(d) {
    data.frame(
      n = nrow(d),
      return_mae = mean(d$return_abs_error, na.rm = TRUE),
      return_rmse = sqrt(mean(d$return_sq_error, na.rm = TRUE)),
      level_mae = mean(d$level_abs_error, na.rm = TRUE),
      level_rmse = sqrt(mean(d$level_sq_error, na.rm = TRUE)),
      coverage_90 = mean(as.numeric(d$covered_90), na.rm = TRUE),
      direction_accuracy = mean(as.numeric(d$direction_correct), na.rm = TRUE),
      mean_p_actual = mean(d$p_actual, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  overall <- summarize(predictions)
  by_asset <- do.call(rbind, lapply(split(predictions, predictions$asset), function(d) cbind(asset = d$asset[1L], summarize(d), stringsAsFactors = FALSE)))
  rownames(by_asset) <- NULL
  by_horizon <- do.call(rbind, lapply(split(predictions, predictions$horizon), function(d) cbind(horizon = d$horizon[1L], summarize(d), stringsAsFactors = FALSE)))
  rownames(by_horizon) <- NULL
  key <- paste(predictions$asset, predictions$horizon, sep = "\r")
  by_asset_horizon <- do.call(rbind, lapply(split(predictions, key), function(d) cbind(asset = d$asset[1L], horizon = d$horizon[1L], summarize(d), stringsAsFactors = FALSE)))
  rownames(by_asset_horizon) <- NULL
  list(overall = overall, by_asset = by_asset, by_horizon = by_horizon, by_asset_horizon = by_asset_horizon)
}

.v5_parallel_export_names <- function(export_env) {
  candidates <- c(
    "craft_fit", "trajectory_embedding", "trajectory_reframer",
    "svd_changepoint_regimes", "trajectory_forecast", "smooth_quantile_distribution",
    "time_volume_weights", "denoiser",
    ".align_named_probs", ".sanitize_regime_weights", ".prob_entropy01", ".effective_n",
    ".align_trajectory_columns", ".v3_normalize_prob", ".v3_softmax", ".numeric_summary",
    ".tail_penalty_scores", ".regime_sampler_tail_penalized",
    ".transition_row_weights", ".transition_row_key", ".transition_prob_to_class",
    ".transition_make_split", ".transition_classification_metrics",
    ".transition_valid_joint_accuracy", ".fit_transition_target", ".fit_transition_classifier",
    ".v5_backtest_horizons", ".v5_select_backtest_origins", ".v5_bind_rows",
    ".v5_backtest_one_origin", ".v5_backtest_summary"
  )
  candidates[vapply(candidates, exists, logical(1), envir = export_env, inherits = TRUE)]
}

.v5_run_internal_backtest <- function(series, volumes, n_testing, fit_args, window, horizon, min_segment, seed, backtest_parallel = TRUE, backtest_workers = NULL, backtest_min_train = NULL) {
  horizons <- .v5_backtest_horizons(window, horizon)
  origins <- .v5_select_backtest_origins(n_obs = nrow(series), n_testing = n_testing, window = window, horizon = horizon, min_segment = min_segment, backtest_min_train = backtest_min_train)
  settings <- list(requested_n_testing = n_testing, actual_n_testing = length(origins), origins = origins, horizons = horizons, parallel = isTRUE(backtest_parallel), workers = NA_integer_, min_train = if (is.null(backtest_min_train)) NA_integer_ else as.integer(backtest_min_train)[1L])
  if (!length(origins)) {
    return(list(settings = settings, predictions = data.frame(), failures = data.frame(), point_results = list(), summary = .v5_backtest_summary(data.frame())))
  }
  tasks <- lapply(seq_along(origins), function(i) list(point_id = i, origin = origins[i], series = series, volumes = volumes, horizons = horizons, fit_args = fit_args))
  use_parallel <- isTRUE(backtest_parallel) && length(tasks) > 1L
  if (use_parallel) {
    workers <- if (is.null(backtest_workers)) max(1L, min(length(tasks), parallel::detectCores(logical = TRUE) - 1L)) else as.integer(backtest_workers)[1L]
    if (!is.finite(workers) || workers < 1L) workers <- 1L
    workers <- min(workers, length(tasks))
  } else {
    workers <- 1L
  }
  settings$workers <- workers
  if (use_parallel && workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    export_env <- parent.env(environment())
    export_names <- .v5_parallel_export_names(export_env)
    parallel::clusterExport(cl, varlist = export_names, envir = export_env)
    results <- parallel::parLapply(cl, tasks, .v5_backtest_one_origin)
  } else {
    results <- lapply(tasks, .v5_backtest_one_origin)
  }
  predictions <- .v5_bind_rows(lapply(results, function(x) if (isTRUE(x$ok)) x$predictions else NULL))
  failures <- .v5_bind_rows(lapply(results, function(x) {
    if (isTRUE(x$ok)) return(NULL)
    data.frame(point_id = x$point_id, origin_index = x$origin_index, error = x$error, stringsAsFactors = FALSE)
  }))
  list(settings = settings, predictions = predictions, failures = failures, point_results = results, summary = .v5_backtest_summary(predictions))
}
# ============================================================
# Main function: craft_fit()
# ============================================================

craft_fit <- function(
  series, volumes = NULL,
  window = 20,
  n_draws = 1000,
  n_factors = 3, 
  var_threshold = 0.95,
  horizon = NULL,
  seed = 42,

  # Transition engine parameters
  transition_alpha = 1,
  transition_min_count = 1,
  transition_key_sep = "\r",

  # Validation
  valid_frac = 0.20,
  valid_split = c("last", "random", "blocked"),

  # Internal backtesting
  n_testing = 0L,
  backtest_parallel = TRUE,
  backtest_workers = NULL,
  backtest_min_train = NULL,

  # Weight used by regime_sampler()
  sampler_weight = c("uniform",
    "weighted_bal_acc",
    "bal_acc",
    "weighted_acc",
    "acc"
  ),
  tail_penalty_lambda = 0,
  tail_penalty_threshold = 3,
  tail_penalty_power = 2,
  tail_penalty_aggregation = c("mean", "max", "sum"),
  tail_penalty_scope = c("forecast", "all"),
  volume_alpha = 1,  max_regimes_per_factor = 6,
  factor_selection = "fixed",
  min_segment = 30L,
  past_regime_args = list(),
  future_regime_args = list(),

  # Diagnostics
  return_train_probs = TRUE,
  verbose = TRUE
) {

  sampler_weight <- match.arg(sampler_weight)
  valid_split <- match.arg(valid_split)
  tail_penalty_aggregation <- match.arg(tail_penalty_aggregation)
  tail_penalty_scope <- match.arg(tail_penalty_scope)

  transition_alpha <- as.numeric(transition_alpha)[1L]
  transition_min_count <- as.numeric(transition_min_count)[1L]

  if (!is.finite(transition_alpha) || transition_alpha < 0) {
    stop("`transition_alpha` must be a finite non-negative number.")
  }

  if (!is.finite(transition_min_count) || transition_min_count < 0) {
    stop("`transition_min_count` must be a finite non-negative number.")
  }

  if (!is.character(transition_key_sep) || length(transition_key_sep) != 1L || is.na(transition_key_sep)) {
    stop("`transition_key_sep` must be a single non-NA character string.")
  }

  n_testing <- as.integer(n_testing)[1L]

  if (!is.finite(n_testing) || n_testing < 0L) {
    stop("`n_testing` must be a finite non-negative integer.")
  }

  if (!is.null(backtest_workers)) {
    backtest_workers <- as.integer(backtest_workers)[1L]

    if (!is.finite(backtest_workers) || backtest_workers < 1L) {
      stop("`backtest_workers` must be NULL or a positive integer.")
    }
  }

  if (!is.logical(backtest_parallel) || length(backtest_parallel) != 1L || is.na(backtest_parallel)) {
    stop("`backtest_parallel` must be TRUE or FALSE.")
  }

  if (!is.null(backtest_min_train)) {
    backtest_min_train <- as.integer(backtest_min_train)[1L]

    if (!is.finite(backtest_min_train) || backtest_min_train < 1L) {
      stop("`backtest_min_train` must be NULL or a positive integer.")
    }
  }

  # ============================================================
  # 0) checks
  # ============================================================
  if (!is.matrix(series) && !is.data.frame(series)) {
    stop("`series` must be a matrix or data.frame.")
  }

  series <- as.data.frame(series, stringsAsFactors = FALSE)
  series[] <- lapply(series, as.numeric)

  if (nrow(series) <= window) {
    stop("`nrow(series)` must be greater than `window`.")
  }

  if (window < 1) {
    stop("`window` must be >= 1.")
  }

  if (n_draws < 1) {
    stop("`n_draws` must be >= 1.")
  }

  tail_penalty_lambda <- as.numeric(tail_penalty_lambda)[1L]
  tail_penalty_threshold <- as.numeric(tail_penalty_threshold)[1L]
  tail_penalty_power <- as.numeric(tail_penalty_power)[1L]

  if (!is.finite(tail_penalty_lambda) || tail_penalty_lambda < 0) {
    stop("`tail_penalty_lambda` must be a finite non-negative number.")
  }

  if (!is.finite(tail_penalty_threshold) || tail_penalty_threshold < 0) {
    stop("`tail_penalty_threshold` must be a finite non-negative number.")
  }

  if (!is.finite(tail_penalty_power) || tail_penalty_power <= 0) {
    stop("`tail_penalty_power` must be a finite positive number.")
  }

  all_horizons <- is.null(horizon)
  if (!all_horizons) {
    horizon <- as.integer(horizon)

    if (length(horizon) != 1L || !is.finite(horizon) || horizon < 1L) {
      stop("`horizon` must be NULL or a positive integer.")
    }

    if (horizon > window) {
      stop("`horizon` cannot exceed `window`.")
    }
  }
  
  n_feats <- ncol(series)
  
  if(all_horizons){horizon_index <- 1:(n_feats * window)}
  if(!all_horizons){horizon_index <- seq(horizon, n_feats * window, window)}

  if (is.null(colnames(series))) {
    colnames(series) <- paste0("series_", seq_len(n_feats))
  }

  if (!is.null(volumes)) {
    if (!is.matrix(volumes) && !is.data.frame(volumes)) {
      stop("`volumes` must be NULL, a matrix, or a data.frame.")
    }

    volumes <- as.data.frame(volumes, stringsAsFactors = FALSE)
    volumes[] <- lapply(volumes, as.numeric)

    if (nrow(volumes) != nrow(series)) {
      stop("`volumes` must have the same number of rows as `series`.")
    }

    if (ncol(volumes) != ncol(series)) {
      stop("`volumes` must have the same number of columns as `series`.")
    }

    colnames(volumes) <- colnames(series)
  }

  # ============================================================
  # 1) trajectory_embedding
  # ============================================================

  trajectory_model <- trajectory_embedding(
    series,
    window
  )
  
  # ============================================================
  # 2) multiregime encoding
  # ============================================================

  past_regime_call <- modifyList(
    list(
      X = trajectory_model$past_trajectories,
      reduction_criteria = factor_selection,
      n_dim = n_factors,
	  var_threshold = var_threshold,
      max_regimes_per_dim = max_regimes_per_factor,
      minseglen = min_segment
    ),
    past_regime_args
  )

  future_regime_call <- modifyList(
    list(
      X = trajectory_model$future_trajectories[, horizon_index],
      reduction_criteria = factor_selection,
      n_dim = n_factors,
	  var_threshold = var_threshold,
      max_regimes_per_dim = max_regimes_per_factor,
      minseglen = min_segment
    ),
    future_regime_args
  )

  past_regime_model <- do.call(svd_changepoint_regimes, past_regime_call)

  future_regime_model <- do.call(svd_changepoint_regimes, future_regime_call)

  past_labels <- past_regime_model$regimes

  future_labels <- future_regime_model$regimes

  latest_regime_labels_raw <- past_regime_model$transform_newdata(trajectory_model$latest_trajectory)
  latest_regime_labels <- latest_regime_labels_raw$regimes

  past_labels <- as.data.frame(past_labels, stringsAsFactors = FALSE)
  future_labels <- as.data.frame(future_labels, stringsAsFactors = FALSE)
  latest_regime_labels <- as.data.frame(latest_regime_labels, stringsAsFactors = FALSE)

  # ============================================================
  # 2b) enforce regime-column alignment
  # ============================================================

  past_names <- names(past_labels)
  latest_names <- names(latest_regime_labels)

  if (!setequal(past_names, latest_names)) {
    stop(
      "`past_labels` and `latest_regime_labels` do not contain the same regime columns.\n",
      "Only in past_labels: ",
      paste(setdiff(past_names, latest_names), collapse = ", "),
      "\nOnly in latest_regime_labels: ",
      paste(setdiff(latest_names, past_names), collapse = ", ")
    )
  }

  latest_regime_labels <- latest_regime_labels[, past_names, drop = FALSE]

  # Stable factor levels for X
  past_labels[] <- lapply(past_labels, function(z) {
    factor(z, levels = sort(unique(z[!is.na(z)])))
  })

  latest_unseen_counts <- vapply(seq_along(latest_regime_labels), function(j) {
    z <- as.character(latest_regime_labels[[j]])
    sum(!is.na(z) & !z %in% levels(past_labels[[j]]))
  }, integer(1))

  names(latest_unseen_counts) <- past_names

  if (any(latest_unseen_counts > 0L) && isTRUE(verbose)) {
    warning(
      "`latest_regime_labels` contains regime values not seen during training; ",
      "they will be encoded as all-zero inputs. Counts: ",
      paste(
        paste0(names(latest_unseen_counts), "=", latest_unseen_counts),
        collapse = ", "
      )
    )
  }

  # Stable factor levels for Y
  future_labels[] <- lapply(future_labels, function(z) {
    factor(z, levels = sort(unique(z[!is.na(z)])))
  })

  # New X must use training X levels
  latest_regime_labels[] <- lapply(seq_along(latest_regime_labels), function(j) {
    factor(latest_regime_labels[[j]], levels = levels(past_labels[[j]]))
  })

  latest_regime_labels <- as.data.frame(latest_regime_labels, stringsAsFactors = FALSE)
  names(latest_regime_labels) <- past_names

  # ============================================================
  # 3) Torch bridge
  # ============================================================
  row_weights <- NULL
  if (!is.null(volumes)) {
    row_weights <- time_volume_weights(volumes, alpha = volume_alpha)
    row_weights <- row_weights[trajectory_model$time_index]
  }

  classifier <- .fit_transition_classifier(
    X = past_labels,
    Y = future_labels,
    row_weights = row_weights,
    valid_frac = valid_frac,
    valid_split = valid_split,
    seed = seed,
    transition_alpha = transition_alpha,
    transition_min_count = transition_min_count,
    transition_key_sep = transition_key_sep
  )

  latest_regime_probabilities <- classifier$predict(
    latest_regime_labels,
    type = "prob"
  )

  if (isTRUE(return_train_probs)) {
    train_probs <- classifier$predict(
      past_labels,
      type = "prob"
    )
  } else {
    train_probs <- NULL
  }

  target_names <- names(future_labels)

  # ============================================================
  # 3b) enforce posterior-list name alignment
  # ============================================================

  latest_regime_probabilities <- .align_named_probs(
    x = latest_regime_probabilities,
    target_names = target_names,
    object_name = "latest_regime_probabilities"
  )

  if (!is.null(train_probs)) {
    train_probs <- .align_named_probs(
      x = train_probs,
      target_names = target_names,
      object_name = "train_probs"
    )
  }

  # ============================================================
  # 3c) derive sampling weights from transition diagnostics
  # ============================================================

  classifier_metrics <- classifier$diagnostics$metrics

  if (sampler_weight == "bal_acc") {
    sampler_weights <- classifier_metrics$balanced_accuracy
  }

  if (sampler_weight == "weighted_bal_acc") {
    sampler_weights <- classifier_metrics$weighted_balanced_accuracy
  }

  if (sampler_weight == "acc") {
    sampler_weights <- classifier_metrics$accuracy
  }

  if (sampler_weight == "weighted_acc") {
    sampler_weights <- classifier_metrics$weighted_accuracy
  }

  if (sampler_weight == "uniform") {
    sampler_weights <- rep(1, length(target_names))
  }

  names(sampler_weights) <- classifier_metrics$target

  sampler_weights <- sampler_weights[target_names]

  sampler_weights <- .sanitize_regime_weights(
    w = sampler_weights,
    target_names = target_names
  )

  # ============================================================
  # 3d) derive row-level soft tail penalty for candidate futures
  # ============================================================

  tail_penalty_cols <- if (tail_penalty_scope == "all") {
    seq_len(ncol(trajectory_model$future_trajectories))
  } else {
    horizon_index
  }

  tail_penalty_info <- .tail_penalty_scores(
    X = trajectory_model$future_trajectories[, tail_penalty_cols, drop = FALSE],
    threshold = tail_penalty_threshold,
    power = tail_penalty_power,
    aggregation = tail_penalty_aggregation
  )

  # ============================================================
  # 4) sample future trajectories  # ============================================================

  sampler_draw <- .regime_sampler_tail_penalized(
    trafo_df = trajectory_model$future_trajectories,
    cluster_df = future_labels,
    prob_vectors = latest_regime_probabilities,
    weights = sampler_weights,
    tail_penalty = tail_penalty_info$scores,
    tail_penalty_lambda = tail_penalty_lambda,
    size = n_draws,
    seed = seed
  )

  trajectory_draws <- sampler_draw$draws
  trajectory_draws <- as.data.frame(trajectory_draws)
  trajectory_alignment <- .align_trajectory_columns(
    trajectory_draws = trajectory_draws,
    asset_names = colnames(series),
    window = window
  )

  trajectory_draws <- trajectory_alignment$trajectory_draws
  
  corr_by_horizon <- lapply(seq_len(window), function(h) {
    cols <- paste0(colnames(series), "_cum_lead_", h)
    stats::cor(trajectory_draws[, cols, drop = FALSE])
  })
  names(corr_by_horizon) <- paste0("t", 1:window)
  corr_by_horizon <- lapply(corr_by_horizon, function(mat) {rownames(mat) <- colnames(mat) <- colnames(series); return(mat)}) 
  
  # ============================================================
  # 5) forecast reconstruction
  # ============================================================

  level_values <- as.numeric(tail(series, 1))

  if (isTRUE(all_horizons)) {

    split_draws <- lapply(colnames(series), function(asset) {
      cols <- paste0(asset, "_cum_lead_", seq_len(window))
      trajectory_draws[, cols, drop = FALSE]
    })

    names(split_draws) <- colnames(series)

    level_dists <- lapply(seq_len(n_feats), function(j) {
      trajectory_forecast(
        split_draws[[j]],
        last_level = level_values[j],
        recon_mode = "scaled_growth"
      )
    })

    return_dists <- lapply(seq_len(n_feats), function(j) {
      trajectory_forecast(
        split_draws[[j]],
        last_level = NULL,
        recon_mode = NULL
      )
    })

    names(level_dists) <- colnames(series)
    names(return_dists) <- colnames(series)
  }

  if (!isTRUE(all_horizons)) {

    split_draws <- lapply(colnames(series), function(asset) {
      col <- paste0(asset, "_cum_lead_", horizon)
      trajectory_draws[, col, drop = FALSE]
    })

    names(split_draws) <- colnames(series)

    level_dists <- lapply(seq_len(n_feats), function(j) {
      trajectory_forecast(
        split_draws[[j]],
        last_level = level_values[j],
        recon_mode = "scaled_growth"
      )
    })

    return_dists <- lapply(seq_len(n_feats), function(j) {
      trajectory_forecast(
        split_draws[[j]],
        last_level = NULL,
        recon_mode = NULL
      )
    })

    # Keep a single named horizon layer, as in your original design
    return_dists <- lapply(seq_len(n_feats), function(j) {
      x <- return_dists[[j]]
      names(x) <- paste0("t", horizon)
      x
    })

    level_dists <- lapply(seq_len(n_feats), function(j) {
      x <- level_dists[[j]]
      names(x) <- paste0("t", horizon)
      x
    })

    names(return_dists) <- colnames(series)
    names(level_dists) <- colnames(series)
  }

  # ============================================================
  # 6) diagnostics summary
  # ============================================================

  valid_acc <- mean(
    classifier_metrics$accuracy,
    na.rm = TRUE
  )

  valid_bal_acc <- mean(
    classifier_metrics$balanced_accuracy,
    na.rm = TRUE
  )

  valid_joint_acc <- classifier$diagnostics$valid_joint_acc

  posterior_entropy <- vapply(
    latest_regime_probabilities,
    function(p) mean(.prob_entropy01(p), na.rm = TRUE),
    numeric(1)
  )

  class_balance <- lapply(future_labels, function(z) {
    tab <- table(z, useNA = "ifany")
    data.frame(
      class = names(tab),
      n = as.integer(tab),
      prop = as.numeric(tab) / sum(tab),
      stringsAsFactors = FALSE
    )
  })

  diagnostics <- list(
    latest_unseen_counts = latest_unseen_counts,
    latest_unseen_total = sum(latest_unseen_counts),
    posterior_entropy = posterior_entropy,
    sampler_weight_effective_n = .effective_n(sampler_weights),
    trajectory_columns_used_names = trajectory_alignment$used_names,
    expected_trajectory_columns = trajectory_alignment$expected_names,
    class_balance = class_balance,
    transition = classifier$diagnostics$transition,
    valid_split = valid_split,    temperature = classifier$diagnostics$temperature,
    factor_selection = factor_selection,
    min_segment = min_segment,
    tail_penalty = list(
      enabled = tail_penalty_lambda > 0,
      lambda = tail_penalty_lambda,
      threshold = tail_penalty_threshold,
      power = tail_penalty_power,
      aggregation = tail_penalty_aggregation,
      scope = tail_penalty_scope,
      columns = tail_penalty_info$columns,
      candidate_summary = tail_penalty_info$summary,
      selected_summary = .numeric_summary(
        tail_penalty_info$scores[sampler_draw$index]
      ),
      sampling_probability_summary = .numeric_summary(
        sampler_draw$sampling_prob
      )
    )
  )

  backtest <- NULL

  if (n_testing > 0L) {
    backtest_fit_args <- list(
      window = window,
      n_draws = n_draws,
      n_factors = n_factors,
      var_threshold = var_threshold,
      horizon = horizon,
      seed = seed,
      transition_alpha = transition_alpha,
      transition_min_count = transition_min_count,
      transition_key_sep = transition_key_sep,
      valid_frac = valid_frac,
      valid_split = valid_split,
      n_testing = 0L,
      backtest_parallel = FALSE,
      backtest_workers = 1L,
      backtest_min_train = backtest_min_train,
      sampler_weight = sampler_weight,
      tail_penalty_lambda = tail_penalty_lambda,
      tail_penalty_threshold = tail_penalty_threshold,
      tail_penalty_power = tail_penalty_power,
      tail_penalty_aggregation = tail_penalty_aggregation,
      tail_penalty_scope = tail_penalty_scope,
      volume_alpha = volume_alpha,
      max_regimes_per_factor = max_regimes_per_factor,
      factor_selection = factor_selection,
      min_segment = min_segment,
      past_regime_args = past_regime_args,
      future_regime_args = future_regime_args,
      return_train_probs = FALSE,
      verbose = FALSE
    )

    backtest <- .v5_run_internal_backtest(
      series = series,
      volumes = volumes,
      n_testing = n_testing,
      fit_args = backtest_fit_args,
      window = window,
    horizon = horizon,
      min_segment = min_segment,
      seed = seed,
      backtest_parallel = backtest_parallel,
      backtest_workers = backtest_workers,
      backtest_min_train = backtest_min_train
    )

    diagnostics$backtest <- backtest$summary
    diagnostics$backtest_settings <- backtest$settings
  }

  # ============================================================
  # 7) output
  # ============================================================

  out <- list(
    call = match.call(),

    window = window,
    horizon = horizon,
    last_level = stats::setNames(level_values, colnames(series)),

    trajectories = trajectory_model,
    past_regimes = past_regime_model,
    future_regimes = future_regime_model,

    past_labels = past_labels,
    future_labels = future_labels,
    latest_regime_labels = latest_regime_labels,

    targets = target_names,
    transition = classifier$diagnostics$transition,

    latest_regime_probabilities = latest_regime_probabilities,    train_probs = train_probs,

    sampler_weights = sampler_weights,
    sampler_weight = sampler_weight,

    trajectory_draws = trajectory_draws,
    sampler_trace = list(
      index = sampler_draw$index,
      probability = sampler_draw$sampling_prob,
      base_log_score = sampler_draw$base_log_score,
      final_log_score = sampler_draw$final_log_score,
      tail_penalty = sampler_draw$tail_penalty
    ),
    tail_penalty = tail_penalty_info,
	corr_by_horizon = corr_by_horizon,
    return_dists = return_dists,
    level_dists = level_dists,

    classifier = classifier,
    classifier_diagnostics = classifier$diagnostics,
    diagnostics = diagnostics,
    backtest = backtest,

    valid_acc = valid_acc,
    valid_bal_acc = valid_bal_acc,
    valid_joint_acc = valid_joint_acc

  )

  class(out) <- c(
    "craft_fit",
    "regime_forecast_v5",
    "regime_forecast_transition",
    "regime_forecast"
  )

  out
}
craft_predict <- function(
  object,
  newdata = NULL,
  n_draws = NULL,
  seed = NULL,
  type = c("all", "trajectory_draws", "distributions", "probabilities", "labels")
) {
  type <- match.arg(type)

  if (!inherits(object, "regime_forecast_v5")) {
    stop("`object` must be a fitted CRAFT object returned by craft_fit().")
  }

  asset_names <- names(object$return_dists)
  if (is.null(asset_names) || !length(asset_names)) {
    asset_names <- unique(sub("_cum_lead_[0-9]+$", "", colnames(object$trajectory_draws)))
  }

  if (is.null(n_draws)) {
    n_draws <- if (is.data.frame(object$trajectory_draws)) nrow(object$trajectory_draws) else 1000L
  }
  n_draws <- as.integer(n_draws)[1L]
  if (!is.finite(n_draws) || n_draws < 1L) {
    stop("`n_draws` must be a positive integer.")
  }

  if (is.null(newdata)) {
    latest_regime_labels <- object$latest_regime_labels
    level_values <- object$last_level
    if (is.null(level_values)) {
      level_values <- rep(NA_real_, length(asset_names))
    }
  } else {
    if (!is.matrix(newdata) && !is.data.frame(newdata)) {
      stop("`newdata` must be NULL, a matrix, or a data.frame.")
    }

    newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
    newdata[] <- lapply(newdata, as.numeric)

    missing_cols <- setdiff(asset_names, names(newdata))
    if (length(missing_cols)) {
      stop("`newdata` is missing columns: ", paste(missing_cols, collapse = ", "))
    }

    newdata <- newdata[, asset_names, drop = FALSE]

    if (nrow(newdata) <= object$window) {
      stop("`newdata` must contain more rows than the fitted trajectory window.")
    }

    trajectory_model <- trajectory_embedding(newdata, object$window)
    latest_regime_labels_raw <- object$past_regimes$transform_newdata(trajectory_model$latest_trajectory)
    latest_regime_labels <- latest_regime_labels_raw$regimes
    level_values <- as.numeric(tail(newdata, 1))
  }

  past_labels <- as.data.frame(object$past_labels, stringsAsFactors = FALSE)
  latest_regime_labels <- as.data.frame(latest_regime_labels, stringsAsFactors = FALSE)

  past_names <- names(past_labels)
  latest_names <- names(latest_regime_labels)

  if (!setequal(past_names, latest_names)) {
    stop(
      "Fitted past labels and prediction labels do not contain the same regime columns."
    )
  }

  latest_regime_labels <- latest_regime_labels[, past_names, drop = FALSE]

  past_labels[] <- lapply(past_labels, function(z) {
    factor(z, levels = sort(unique(z[!is.na(z)])))
  })

  latest_unseen_counts <- vapply(seq_along(latest_regime_labels), function(j) {
    z <- as.character(latest_regime_labels[[j]])
    sum(!is.na(z) & !z %in% levels(past_labels[[j]]))
  }, integer(1))
  names(latest_unseen_counts) <- past_names

  latest_regime_labels[] <- lapply(seq_along(latest_regime_labels), function(j) {
    factor(latest_regime_labels[[j]], levels = levels(past_labels[[j]]))
  })

  latest_regime_labels <- as.data.frame(latest_regime_labels, stringsAsFactors = FALSE)
  names(latest_regime_labels) <- past_names

  latest_regime_probabilities <- object$classifier$predict(latest_regime_labels, type = "prob")
  latest_regime_probabilities <- .align_named_probs(
    x = latest_regime_probabilities,
    target_names = object$targets,
    object_name = "latest_regime_probabilities"
  )

  if (type == "labels") {
    return(latest_regime_labels)
  }

  if (type == "probabilities") {
    return(latest_regime_probabilities)
  }

  tail_penalty_lambda <- 0
  if (!is.null(object$diagnostics$tail_penalty$lambda)) {
    tail_penalty_lambda <- object$diagnostics$tail_penalty$lambda
  }

  sampler_draw <- .regime_sampler_tail_penalized(
    trafo_df = object$trajectories$future_trajectories,
    cluster_df = object$future_labels,
    prob_vectors = latest_regime_probabilities,
    weights = object$sampler_weights,
    tail_penalty = object$tail_penalty$scores,
    tail_penalty_lambda = tail_penalty_lambda,
    size = n_draws,
    seed = seed
  )

  trajectory_draws <- as.data.frame(sampler_draw$draws)
  trajectory_alignment <- .align_trajectory_columns(
    trajectory_draws = trajectory_draws,
    asset_names = asset_names,
    window = object$window
  )
  trajectory_draws <- trajectory_alignment$trajectory_draws

  if (type == "trajectory_draws") {
    return(trajectory_draws)
  }

  corr_by_horizon <- lapply(seq_len(object$window), function(h) {
    cols <- paste0(asset_names, "_cum_lead_", h)
    stats::cor(trajectory_draws[, cols, drop = FALSE])
  })
  names(corr_by_horizon) <- paste0("t", seq_len(object$window))
  corr_by_horizon <- lapply(corr_by_horizon, function(mat) {
    rownames(mat) <- colnames(mat) <- asset_names
    mat
  })

  all_horizons <- is.null(object$horizon)

  if (isTRUE(all_horizons)) {
    split_draws <- lapply(asset_names, function(asset) {
      cols <- paste0(asset, "_cum_lead_", seq_len(object$window))
      trajectory_draws[, cols, drop = FALSE]
    })
    names(split_draws) <- asset_names
  } else {
    split_draws <- lapply(asset_names, function(asset) {
      col <- paste0(asset, "_cum_lead_", object$horizon)
      trajectory_draws[, col, drop = FALSE]
    })
    names(split_draws) <- asset_names
  }

  return_dists <- lapply(seq_along(asset_names), function(j) {
    trajectory_forecast(split_draws[[j]], last_level = NULL, recon_mode = NULL)
  })
  level_dists <- lapply(seq_along(asset_names), function(j) {
    trajectory_forecast(
      split_draws[[j]],
      last_level = level_values[j],
      recon_mode = "scaled_growth"
    )
  })
  names(return_dists) <- asset_names
  names(level_dists) <- asset_names

  if (!isTRUE(all_horizons)) {
    return_dists <- lapply(return_dists, function(x) {
      names(x) <- paste0("t", object$horizon)
      x
    })
    level_dists <- lapply(level_dists, function(x) {
      names(x) <- paste0("t", object$horizon)
      x
    })
  }

  prediction <- list(
    latest_regime_labels = latest_regime_labels,
    latest_regime_probabilities = latest_regime_probabilities,
    trajectory_draws = trajectory_draws,
    sampler_trace = list(
      index = sampler_draw$index,
      probability = sampler_draw$sampling_prob,
      base_log_score = sampler_draw$base_log_score,
      final_log_score = sampler_draw$final_log_score,
      tail_penalty = sampler_draw$tail_penalty
    ),
    corr_by_horizon = corr_by_horizon,
    return_dists = return_dists,
    level_dists = level_dists,
    diagnostics = list(
      latest_unseen_counts = latest_unseen_counts,
      trajectory_columns_used_names = trajectory_alignment$used_names,
      expected_trajectory_columns = trajectory_alignment$expected_names
    )
  )

  class(prediction) <- c("craft_prediction", "regime_forecast_prediction")

  if (type == "distributions") {
    return(prediction[c("return_dists", "level_dists")])
  }

  prediction
}
