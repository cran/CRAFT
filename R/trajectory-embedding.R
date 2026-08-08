trajectory_embedding <- function(ts_set, trajectory_window) {
  feature_names <- colnames(ts_set)

  past_parts <- lapply(feature_names, function(feature_name) {
    trajectory_reframer(
      ts_set[[feature_name]],
      trajectory_window,
      "backward",
      feature_name
    )$trajectory
  })

  future_parts <- lapply(feature_names, function(feature_name) {
    trajectory_reframer(
      ts_set[[feature_name]],
      trajectory_window,
      "forward",
      feature_name
    )$trajectory
  })

  past_trajectories <- do.call(cbind, past_parts)
  future_trajectories <- do.call(cbind, future_parts)
  latest_trajectory <- tail(past_trajectories, 1)

  na_idx <- apply(past_trajectories, 1, anyNA) |
    apply(future_trajectories, 1, anyNA)

  past_trajectories <- past_trajectories[!na_idx, , drop = FALSE]
  future_trajectories <- future_trajectories[!na_idx, , drop = FALSE]

  list(
    past_trajectories = past_trajectories,
    future_trajectories = future_trajectories,
    latest_trajectory = latest_trajectory,
    time_index = !na_idx
  )
}

trajectory_reframer <- function(ts, k, mode = c("backward", "forward"), ts_name = NULL) {
  x <- as.numeric(ts)
  n <- length(x)
  if (is.null(ts_name)) {
    ts_name <- deparse(substitute(ts))
  }

  if (n < 2) {
    stop("ts must have length >= 2.")
  }
  if (k < 1 || k != as.integer(k)) {
    stop("k must be a positive integer.")
  }
  if (k >= n) {
    stop("k must be smaller than length(ts).")
  }

  mode <- match.arg(mode)
  trajectory <- matrix(NA_real_, nrow = n, ncol = k)

  if (mode == "backward") {
    for (t in seq_len(n)) {
      for (h in seq_len(k)) {
        base_idx <- t - h
        if (base_idx >= 1) {
          base <- x[base_idx]
          if (!is.na(base) && base != 0) {
            trajectory[t, h] <- x[t] / base - 1
          } else {
            trajectory[t, h] <- NA_real_
          }
        }
      }
    }
    colnames(trajectory) <- paste0(ts_name, "_cum_lag_", seq_len(k))
  } else {
    for (t in seq_len(n)) {
      base <- x[t]
      if (!is.na(base) && base != 0) {
        for (h in seq_len(k)) {
          fut_idx <- t + h
          if (fut_idx <= n) {
            trajectory[t, h] <- x[fut_idx] / base - 1
          }
        }
      }
    }
    colnames(trajectory) <- paste0(ts_name, "_cum_lead_", seq_len(k))
  }

  out <- list(
    trajectory = trajectory,
    ts = x,
    k = k,
    mode = mode
  )
  class(out) <- c("trajectory_reframer", class(out))
  out
}
