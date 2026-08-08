time_volume_weights <- function(
  volume_df,
  date_col = NULL,
  alpha = 1,
  normalize = TRUE,
  eps = 1e-8,
  winsorize = FALSE,
  winsor_q = 0.99
) {
  # Remove date column if provided
  if (!is.null(date_col)) {
    dates <- volume_df[[date_col]]
    x <- volume_df[, setdiff(names(volume_df), date_col), drop = FALSE]
  } else {
    dates <- NULL
    x <- volume_df
  }

  # Keep numeric columns only
  x <- x[, sapply(x, is.numeric), drop = FALSE]

  if (ncol(x) == 0) {
    stop("No numeric volume columns found.")
  }

  # Clean invalid values
  x[!is.finite(as.matrix(x))] <- NA
  x[x < 0] <- NA

  # Aggregate volumes by time point
  total_volume_t <- rowSums(x, na.rm = TRUE)

  # Optional winsorization to reduce extreme volume spikes
  if (winsorize) {
    cap <- quantile(total_volume_t, winsor_q, na.rm = TRUE)
    total_volume_t <- pmin(total_volume_t, cap)
  }

  # Optional compression/amplification
  # alpha < 1 compresses extreme high-volume days
  # alpha = 1 keeps raw proportional volume
  # alpha > 1 emphasizes high-volume days
  score <- (total_volume_t + eps)^alpha

  # Normalize into weights summing to 1
  if (normalize) {
    weights <- score / sum(score, na.rm = TRUE)
  } else {
    weights <- score
  }

  # Return named vector if dates exist
  if (!is.null(dates)) {
    names(weights) <- as.character(dates)
  }

  weights
}
