denoiser <- function(
  x,
  target = NULL,
  min_non_missing = 0.70,
  min_unique = 3L,
  min_mad = 1e-8,
  max_cor = 0.95,
  keep_quantile = 0.70,
  min_score = NULL,
  min_features = 1L,
  max_features = NULL,
  top_n = NULL,
  supervised_weight = if (!is.null(target)) 0.60 else 0,
  pca_weight = 0.45,
  temporal_weight = 0.25,
  snr_weight = 0.15,
  info_weight = 0.15,
  lag_max = 5L,
  impute = c("median", "zero"),
  return_report = FALSE,
  verbose = TRUE
) {
  impute <- match.arg(impute)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data.frame.")
  }

  input_is_matrix <- is.matrix(x)
  df <- as.data.frame(x, stringsAsFactors = FALSE)

  if (is.null(colnames(df)) || any(colnames(df) == "")) {
    colnames(df) <- paste0("feature_", seq_len(ncol(df)))
  }

  y <- NULL

  if (!is.null(target)) {
    if (is.character(target) && length(target) == 1L && target %in% colnames(df)) {
      y <- df[[target]]
      df[[target]] <- NULL
    } else {
      y <- target
    }

    if (length(y) != nrow(df)) {
      stop("target must have the same number of rows as x.")
    }
  }

  numeric_cols <- vapply(
    df,
    function(z) is.numeric(z) || is.integer(z) || is.logical(z),
    logical(1)
  )

  dropped_non_numeric <- names(df)[!numeric_cols]
  df_num <- df[, numeric_cols, drop = FALSE]

  if (!ncol(df_num)) {
    out <- if (input_is_matrix) {
      matrix(nrow = nrow(df), ncol = 0L)
    } else {
      df_num
    }

    report <- data.frame(
      feature = dropped_non_numeric,
      selected = FALSE,
      score = 0,
      reason = "non_numeric",
      stringsAsFactors = FALSE
    )

    attr(out, "denoiser_report") <- report
    attr(out, "selected_features") <- character()

    if (isTRUE(return_report)) {
      return(list(data = out, features = character(), report = report))
    }

    return(out)
  }

  df_num[] <- lapply(df_num, function(z) as.numeric(z))

  n <- nrow(df_num)
  p <- ncol(df_num)
  feature_names <- colnames(df_num)

  min_non_missing <- min(max(as.numeric(min_non_missing)[1L], 0), 1)
  min_unique <- max(1L, as.integer(min_unique))
  min_mad <- max(as.numeric(min_mad)[1L], 0)
  max_cor <- min(max(as.numeric(max_cor)[1L], 0), 1)
  keep_quantile <- min(max(as.numeric(keep_quantile)[1L], 0), 1)
  min_features <- max(0L, as.integer(min_features))
  lag_max <- max(1L, as.integer(lag_max))

  if (!is.null(max_features)) {
    max_features <- max(1L, as.integer(max_features))
  }

  if (!is.null(top_n)) {
    top_n <- max(1L, as.integer(top_n))
  }

  supervised_weight <- min(max(as.numeric(supervised_weight)[1L], 0), 1)

  .finite <- function(z) {
    z[is.finite(z)]
  }

  .rank01 <- function(z) {
    z <- as.numeric(z)
    out <- rep(0, length(z))
    ok <- is.finite(z)

    if (sum(ok) == 0L) {
      return(out)
    }

    lo <- min(z[ok], na.rm = TRUE)
    hi <- max(z[ok], na.rm = TRUE)

    if (!is.finite(lo) || !is.finite(hi) || abs(hi - lo) <= 1e-12) {
      out[ok] <- 0.5
      return(out)
    }

    out[ok] <- (z[ok] - lo) / (hi - lo)
    out
  }

  .safe_cor <- function(a, b, method = "spearman") {
    a <- as.numeric(a)
    b <- as.numeric(b)
    ok <- is.finite(a) & is.finite(b)

    if (sum(ok) < 4L) {
      return(NA_real_)
    }

    if (stats::sd(a[ok]) <= 1e-12 || stats::sd(b[ok]) <= 1e-12) {
      return(NA_real_)
    }

    suppressWarnings(stats::cor(a[ok], b[ok], method = method))
  }

  .eta_squared <- function(a, group) {
    a <- as.numeric(a)
    group <- as.factor(group)
    ok <- is.finite(a) & !is.na(group)

    if (sum(ok) < 4L || length(unique(group[ok])) < 2L) {
      return(NA_real_)
    }

    a <- a[ok]
    group <- group[ok]
    grand <- mean(a, na.rm = TRUE)
    total <- sum((a - grand)^2, na.rm = TRUE)

    if (!is.finite(total) || total <= 1e-12) {
      return(NA_real_)
    }

    groups <- split(a, group)
    between <- sum(vapply(
      groups,
      function(z) length(z) * (mean(z, na.rm = TRUE) - grand)^2,
      numeric(1)
    ))

    between / total
  }

  .temporal_structure <- function(z, max_lag) {
    z <- as.numeric(z)
    ok <- is.finite(z)

    if (sum(ok) < max(8L, max_lag + 3L)) {
      return(NA_real_)
    }

    vals <- numeric()
    max_lag <- min(max_lag, length(z) - 2L)

    for (lag in seq_len(max_lag)) {
      a <- z[(lag + 1L):length(z)]
      b <- z[seq_len(length(z) - lag)]
      c <- .safe_cor(a, b, method = "spearman")

      if (is.finite(c)) {
        vals <- c(vals, abs(c))
      }
    }

    if (!length(vals)) {
      return(NA_real_)
    }

    max(vals, na.rm = TRUE)
  }

  .fill_one <- function(z) {
    z <- as.numeric(z)
    ok <- is.finite(z)

    fill <- if (impute == "median" && any(ok)) {
      stats::median(z[ok], na.rm = TRUE)
    } else {
      0
    }

    z[!ok] <- fill
    z
  }

  missing_rate <- vapply(df_num, function(z) mean(!is.finite(z)), numeric(1))
  unique_n <- vapply(df_num, function(z) length(unique(.finite(z))), numeric(1))
  sd_raw <- vapply(df_num, function(z) stats::sd(.finite(z)), numeric(1))
  mad_raw <- vapply(
    df_num,
    function(z) stats::mad(.finite(z), constant = 1.4826, na.rm = TRUE),
    numeric(1)
  )

  sd_raw[!is.finite(sd_raw)] <- 0
  mad_raw[!is.finite(mad_raw)] <- 0

  basic_keep <- missing_rate <= (1 - min_non_missing) &
    unique_n >= min_unique &
    mad_raw >= min_mad &
    sd_raw > 0

  X_imp <- as.matrix(as.data.frame(lapply(df_num, .fill_one)))
  colnames(X_imp) <- feature_names

  center <- apply(X_imp, 2, stats::median, na.rm = TRUE)
  scale <- apply(X_imp, 2, stats::mad, constant = 1.4826, na.rm = TRUE)
  fallback_scale <- apply(X_imp, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 1e-12] <- fallback_scale[!is.finite(scale) | scale <= 1e-12]
  scale[!is.finite(scale) | scale <= 1e-12] <- 1

  Z <- sweep(sweep(X_imp, 2, center, "-"), 2, scale, "/")
  Z[!is.finite(Z)] <- 0

  pca_score <- rep(0, p)
  names(pca_score) <- feature_names

  cand_idx <- which(basic_keep)

  if (length(cand_idx) >= 2L && n >= 8L) {
    Zc <- Z[, cand_idx, drop = FALSE]
    C <- suppressWarnings(stats::cor(Zc, use = "pairwise.complete.obs"))
    C[!is.finite(C)] <- 0
    C <- (C + t(C)) / 2
    diag(C) <- 1

    eg <- eigen(C, symmetric = TRUE)
    eig <- pmax(eg$values, 0)
    q <- ncol(Zc) / max(1L, nrow(Zc))
    lambda_plus <- (1 + sqrt(q))^2
    signal_pc <- which(eig > lambda_plus)

    if (length(signal_pc)) {
      loadings <- eg$vectors[, signal_pc, drop = FALSE]
      weighted_loading <- loadings^2 %*% matrix(eig[signal_pc], ncol = 1L)
      pca_score[cand_idx] <- as.numeric(weighted_loading)
    }
  } else if (length(cand_idx) == 1L) {
    pca_score[cand_idx] <- 1
  }

  temporal_score <- vapply(
    seq_len(p),
    function(j) .temporal_structure(Z[, j], lag_max),
    numeric(1)
  )

  snr_score <- vapply(
    df_num,
    function(z) {
      z <- .finite(z)

      if (length(z) < 4L) {
        return(NA_real_)
      }

      s <- stats::sd(z, na.rm = TRUE)

      if (!is.finite(s) || s <= 1e-12) {
        return(NA_real_)
      }

      abs(mean(z, na.rm = TRUE)) / s
    },
    numeric(1)
  )

  info_score <- 0.5 * (1 - missing_rate) +
    0.5 * pmin(unique_n / max(1L, n), 1)

  supervised_score <- rep(NA_real_, p)
  names(supervised_score) <- feature_names

  if (!is.null(y)) {
    y_is_numeric <- is.numeric(y) || is.integer(y)
    y_unique <- length(unique(y[!is.na(y)]))

    supervised_score <- vapply(
      df_num,
      function(z) {
        if (y_is_numeric && y_unique > 10L) {
          abs(.safe_cor(z, y, method = "spearman"))
        } else {
          .eta_squared(z, y)
        }
      },
      numeric(1)
    )
  }

  pca_score01 <- .rank01(pca_score)
  temporal_score01 <- .rank01(temporal_score)
  snr_score01 <- .rank01(snr_score)
  info_score01 <- .rank01(info_score)
  supervised_score01 <- .rank01(supervised_score)

  unsup_weights <- c(
    pca = max(as.numeric(pca_weight)[1L], 0),
    temporal = max(as.numeric(temporal_weight)[1L], 0),
    snr = max(as.numeric(snr_weight)[1L], 0),
    info = max(as.numeric(info_weight)[1L], 0)
  )

  if (sum(unsup_weights) <= 0) {
    unsup_weights[] <- 1
  }

  unsup_weights <- unsup_weights / sum(unsup_weights)

  unsup_score <- unsup_weights["pca"] * pca_score01 +
    unsup_weights["temporal"] * temporal_score01 +
    unsup_weights["snr"] * snr_score01 +
    unsup_weights["info"] * info_score01

  score <- if (!is.null(y) && supervised_weight > 0) {
    (1 - supervised_weight) * unsup_score + supervised_weight * supervised_score01
  } else {
    unsup_score
  }

  score[!basic_keep] <- 0
  score[!is.finite(score)] <- 0
  names(score) <- feature_names

  selected <- rep(FALSE, p)
  names(selected) <- feature_names

  if (length(cand_idx)) {
    candidate_scores <- score[cand_idx]

    if (!is.null(top_n)) {
      keep_n <- min(top_n, length(cand_idx))
      keep_idx <- cand_idx[order(candidate_scores, decreasing = TRUE)][seq_len(keep_n)]
    } else {
      cutoff <- if (!is.null(min_score)) {
        as.numeric(min_score)[1L]
      } else {
        stats::quantile(candidate_scores, probs = keep_quantile, na.rm = TRUE, names = FALSE)
      }

      keep_idx <- cand_idx[candidate_scores >= cutoff]
    }

    if (length(keep_idx) < min(min_features, length(cand_idx))) {
      fallback_n <- min(min_features, length(cand_idx))
      keep_idx <- cand_idx[order(candidate_scores, decreasing = TRUE)][seq_len(fallback_n)]
    }

    if (!is.null(max_features) && length(keep_idx) > max_features) {
      keep_idx <- keep_idx[order(score[keep_idx], decreasing = TRUE)][seq_len(max_features)]
    }

    selected[keep_idx] <- TRUE
  }

  redundant <- rep(FALSE, p)
  names(redundant) <- feature_names

  selected_idx <- which(selected)

  if (length(selected_idx) > 1L && max_cor < 1) {
    ord <- selected_idx[order(score[selected_idx], decreasing = TRUE)]
    kept <- integer()

    for (idx in ord) {
      if (!length(kept)) {
        kept <- c(kept, idx)
        next
      }

      cors <- vapply(
        kept,
        function(k) abs(.safe_cor(Z[, idx], Z[, k], method = "spearman")),
        numeric(1)
      )

      cors[!is.finite(cors)] <- 0

      if (max(cors, na.rm = TRUE) >= max_cor) {
        redundant[idx] <- TRUE
      } else {
        kept <- c(kept, idx)
      }
    }

    selected[] <- FALSE
    selected[kept] <- TRUE
  }

  selected_idx <- which(selected)
  target_min <- min(min_features, length(cand_idx))

  if (length(selected_idx) < target_min && length(cand_idx)) {
    add_order <- cand_idx[order(score[cand_idx], decreasing = TRUE)]
    add_order <- setdiff(add_order, selected_idx)

    for (idx in add_order) {
      selected_idx <- which(selected)

      if (length(selected_idx) >= target_min) {
        break
      }

      if (max_cor < 1 && length(selected_idx)) {
        cors <- vapply(
          selected_idx,
          function(k) abs(.safe_cor(Z[, idx], Z[, k], method = "spearman")),
          numeric(1)
        )

        cors[!is.finite(cors)] <- 0

        if (max(cors, na.rm = TRUE) >= max_cor) {
          redundant[idx] <- TRUE
          next
        }
      }

      selected[idx] <- TRUE
    }
  }

  reason <- rep("low_signal_score", p)
  reason[missing_rate > (1 - min_non_missing)] <- "too_many_missing"
  reason[unique_n < min_unique] <- "too_few_unique_values"
  reason[mad_raw < min_mad | sd_raw <= 0] <- "near_constant"
  reason[redundant] <- "redundant_high_correlation"
  reason[selected] <- "selected"

  report <- data.frame(
    feature = feature_names,
    selected = as.logical(selected),
    score = as.numeric(score),
    pca_score = as.numeric(pca_score01),
    temporal_score = as.numeric(temporal_score01),
    snr_score = as.numeric(snr_score01),
    info_score = as.numeric(info_score01),
    supervised_score = as.numeric(supervised_score01),
    missing_rate = as.numeric(missing_rate),
    unique_n = as.integer(unique_n),
    sd = as.numeric(sd_raw),
    mad = as.numeric(mad_raw),
    reason = reason,
    stringsAsFactors = FALSE
  )

  if (length(dropped_non_numeric)) {
    non_num_report <- data.frame(
      feature = dropped_non_numeric,
      selected = FALSE,
      score = 0,
      pca_score = NA_real_,
      temporal_score = NA_real_,
      snr_score = NA_real_,
      info_score = NA_real_,
      supervised_score = NA_real_,
      missing_rate = NA_real_,
      unique_n = NA_integer_,
      sd = NA_real_,
      mad = NA_real_,
      reason = "non_numeric",
      stringsAsFactors = FALSE
    )

    report <- rbind(report, non_num_report)
  }

  report <- report[order(report$selected, report$score, decreasing = TRUE), , drop = FALSE]
  rownames(report) <- NULL

  selected_features <- feature_names[selected]

  if (input_is_matrix) {
    out <- as.matrix(df_num[, selected_features, drop = FALSE])
  } else {
    out <- df_num[, selected_features, drop = FALSE]
  }

  attr(out, "denoiser_report") <- report
  attr(out, "selected_features") <- selected_features
  attr(out, "removed_features") <- setdiff(feature_names, selected_features)

  if (isTRUE(verbose)) {
    message(
      "denoiser kept ", length(selected_features), "/", length(feature_names),
      " numeric features."
    )
  }

  if (isTRUE(return_report)) {
    return(list(
      data = out,
      features = selected_features,
      report = report,
      settings = list(
        min_non_missing = min_non_missing,
        min_unique = min_unique,
        min_mad = min_mad,
        max_cor = max_cor,
        keep_quantile = keep_quantile,
        min_score = min_score,
        min_features = min_features,
        max_features = max_features,
        top_n = top_n,
        supervised_weight = supervised_weight,
        return_mode = if (!is.null(y)) "supervised_plus_unsupervised" else "unsupervised"
      )
    ))
  }

  out
}
