svd_changepoint_regimes <- function(
  X,
  series_cols = NULL,
  denoise = TRUE,
  denoiser_args = NULL,

  reduction_criteria = c(
    "variance", "elbow", "kaiser", "bic", "fixed",
    "regime_strength", "hybrid"
  ),

  n_dim = NULL,
  var_threshold = 0.95,
  max_dim = NULL,

  center = TRUE,
  scale = TRUE,
  na_action = c("impute", "stop"),

  cpt_method = "PELT",
  cpt_penalty = "MBIC",
  minseglen = 5,

  filter_weak_dims = FALSE,
  min_regime_strength = 0.10,
  max_selected_dims = NULL,

  max_regimes_per_dim = 6,
  min_regime_size = 30,
  merge_similar_regimes = TRUE,
  merge_threshold = 0.35,

  stability_bootstrap = FALSE,
  n_bootstrap = 100,
  bootstrap_noise = 0.05,
  cpt_tolerance = NULL,
  min_cpt_stability = 0.60,

  newdata_temperature = 1,

  regime_prefix = "regime",
  return_scores = TRUE,
  return_probabilities = TRUE,

  plot = FALSE,
  plot_type = c("scores", "heatmap", "summary"),

  seed = NULL
) {
  reduction_criteria <- match.arg(reduction_criteria)
  na_action <- match.arg(na_action)
  plot_type <- match.arg(plot_type)

  if (!requireNamespace("changepoint", quietly = TRUE)) {
    stop("Package 'changepoint' is required. Install it with install.packages('changepoint').")
  }

  if (!is.null(seed)) set.seed(seed)

  X <- as.data.frame(X)

  # ------------------------------------------------------------
  # 1. Select numeric columns
  # ------------------------------------------------------------

  if (is.null(series_cols)) {
    series_cols <- names(X)[vapply(X, is.numeric, logical(1))]
  }

  original_series_cols <- series_cols
  denoiser_report <- NULL
  denoiser_settings <- NULL

  if (length(series_cols) < 2) {
    stop("At least two numeric time-series columns are required.")
  }

  if (isTRUE(denoise)) {
    if (is.null(denoiser_args)) {
      denoiser_args <- list()
    }

    if (!is.list(denoiser_args)) {
      stop("denoiser_args must be NULL or a list of arguments passed to denoiser().")
    }

    denoiser_args$x <- X[, series_cols, drop = FALSE]
    denoiser_args$return_report <- TRUE

    if (is.null(denoiser_args$verbose)) {
      denoiser_args$verbose <- FALSE
    }

    if (is.null(denoiser_args$min_features)) {
      denoiser_args$min_features <- 2L
    }

    denoised <- do.call(denoiser, denoiser_args)

    X_denoised <- denoised$data
    series_cols <- colnames(X_denoised)
    denoiser_report <- denoised$report
    denoiser_settings <- denoised$settings

    if (length(series_cols) < 2) {
      stop("Denoising left fewer than two features. Relax denoiser_args or set denoise = FALSE.")
    }

    X_mat <- as.matrix(X_denoised)
  } else {
    X_mat <- as.matrix(X[, series_cols, drop = FALSE])
  }

  if (!all(vapply(as.data.frame(X_mat), is.numeric, logical(1)))) {
    stop("All selected series columns must be numeric.")
  }

  if (nrow(X_mat) < 2) {
    stop("X must contain at least two rows.")
  }

  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # Store imputation values from training data
  impute_values <- colMeans(X_mat, na.rm = TRUE)
  impute_values[!is.finite(impute_values)] <- 0

  if (anyNA(X_mat)) {
    if (na_action == "stop") {
      stop("Missing values found. Use na_action = 'impute' or clean the data first.")
    }

    for (j in seq_len(ncol(X_mat))) {
      X_mat[is.na(X_mat[, j]), j] <- impute_values[j]
    }
  }

  if (is.null(cpt_tolerance)) {
    cpt_tolerance <- max(2, minseglen)
  }

  # ------------------------------------------------------------
  # 2. Scale data
  # ------------------------------------------------------------

  center_vec <- if (center) colMeans(X_mat) else rep(0, p)

  scale_vec <- if (scale) {
    s <- apply(X_mat, 2, sd)
    s[!is.finite(s) | s == 0] <- 1
    s
  } else {
    rep(1, p)
  }

  X_scaled <- sweep(X_mat, 2, center_vec, "-")
  X_scaled <- sweep(X_scaled, 2, scale_vec, "/")

  # ------------------------------------------------------------
  # 3. SVD
  # ------------------------------------------------------------

  sv <- svd(X_scaled)

  singular_values <- sv$d
  rank_max <- length(singular_values)

  if (is.null(max_dim)) {
    max_dim <- rank_max
  } else {
    max_dim <- min(max_dim, rank_max)
  }

  eig <- singular_values^2

  if (sum(eig) <= .Machine$double.eps) {
    explained_var <- rep(0, length(eig))
    explained_var[1] <- 1
  } else {
    explained_var <- eig / sum(eig)
  }

  cum_explained_var <- cumsum(explained_var)

  V_all <- sv$v

  # Deterministic sign correction
  for (j in seq_len(ncol(V_all))) {
    if (sum(V_all[, j], na.rm = TRUE) < 0) {
      V_all[, j] <- -V_all[, j]
    }
  }

  # ------------------------------------------------------------
  # 4. Initial SVD dimension selection
  # ------------------------------------------------------------

  select_initial_k <- function() {
    if (reduction_criteria == "fixed") {
      if (is.null(n_dim)) {
        stop("When reduction_criteria = 'fixed', provide n_dim.")
      }

      return(max(1, min(n_dim, max_dim)))
    }

    if (reduction_criteria == "variance") {
      k <- which(cum_explained_var >= var_threshold)[1]
      if (is.na(k)) k <- max_dim
      return(max(1, min(k, max_dim)))
    }

    if (reduction_criteria == "kaiser") {
      k <- sum(eig > mean(eig))
      return(max(1, min(k, max_dim)))
    }

    if (reduction_criteria == "elbow") {
      if (max_dim <= 2) return(1)

      y <- cum_explained_var[seq_len(max_dim)]
      x <- seq_along(y)

      x1 <- x[1]
      y1 <- y[1]
      x2 <- x[length(x)]
      y2 <- y[length(y)]

      dist <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
        sqrt((y2 - y1)^2 + (x2 - x1)^2)

      k <- which.max(dist)

      return(max(1, min(k, max_dim)))
    }

    if (reduction_criteria == "bic") {
      bic_values <- rep(NA_real_, max_dim)

      for (k in seq_len(max_dim)) {
        rss <- if (k < rank_max) {
          sum(eig[(k + 1):rank_max])
        } else {
          0
        }

        rss <- max(rss, .Machine$double.eps)

        n_params <- k * (n + p - k)

        bic_values[k] <- n * p * log(rss / (n * p)) +
          log(n * p) * n_params
      }

      k <- which.min(bic_values)

      return(max(1, min(k, max_dim)))
    }

    if (reduction_criteria %in% c("regime_strength", "hybrid")) {
      return(max_dim)
    }
  }

  k_initial <- select_initial_k()
  initial_dims <- seq_len(k_initial)

  scores_initial <- X_scaled %*% V_all[, initial_dims, drop = FALSE]
  colnames(scores_initial) <- paste0("svd_dim_", initial_dims)

  # ------------------------------------------------------------
  # 5. Helpers
  # ------------------------------------------------------------

  run_cpt <- function(z) {
    z <- as.numeric(z)

    if (length(z) < 2 * minseglen || sd(z) == 0 || all(z == z[1])) {
      return(NULL)
    }

    tryCatch(
      changepoint::cpt.meanvar(
        z,
        method = cpt_method,
        penalty = cpt_penalty,
        minseglen = minseglen
      ),
      error = function(e) NULL
    )
  }

  make_segments <- function(cpts, n_obs) {
    cpts <- sort(unique(as.integer(cpts)))
    cpts <- cpts[cpts > 0 & cpts < n_obs]

    starts <- c(1, cpts + 1)
    ends <- c(cpts, n_obs)

    data.frame(
      regime = seq_along(starts),
      start = starts,
      end = ends,
      stringsAsFactors = FALSE
    )
  }

  segments_to_regime <- function(segments, n_obs) {
    r <- integer(n_obs)

    for (i in seq_len(nrow(segments))) {
      r[segments$start[i]:segments$end[i]] <- i
    }

    r
  }

  profile_from_segments <- function(z, segments, dim_id) {
    prof <- do.call(
      rbind,
      lapply(seq_len(nrow(segments)), function(i) {
        idx <- segments$start[i]:segments$end[i]
        zz <- z[idx]

        data.frame(
          dim = dim_id,
          regime = i,
          start = segments$start[i],
          end = segments$end[i],
          n = length(zz),
          mean = mean(zz),
          sd = sd(zz),
          stringsAsFactors = FALSE
        )
      })
    )

    global_sd <- sd(z)
    if (!is.finite(global_sd) || global_sd == 0) global_sd <- 1

    prof$sd[!is.finite(prof$sd) | prof$sd == 0] <- global_sd

    prof
  }

  segment_distance <- function(z, segments, i) {
    prof <- profile_from_segments(z, segments, dim_id = 1)

    m1 <- prof$mean[i]
    m2 <- prof$mean[i + 1]

    s1 <- prof$sd[i]
    s2 <- prof$sd[i + 1]

    global_sd <- sd(z)
    if (!is.finite(global_sd) || global_sd == 0) global_sd <- 1

    mean_dist <- abs(m1 - m2) / global_sd
    sd_dist <- abs(log(s1 / s2))

    sqrt(mean_dist^2 + sd_dist^2)
  }

  merge_pair <- function(segments, i) {
    segments$end[i] <- segments$end[i + 1]
    segments <- segments[-(i + 1), , drop = FALSE]
    segments$regime <- seq_len(nrow(segments))
    rownames(segments) <- NULL
    segments
  }

  postprocess_segments <- function(z, segments) {
    if (nrow(segments) <= 1) {
      return(segments)
    }

    # Merge short regimes
    repeat {
      lens <- segments$end - segments$start + 1

      if (nrow(segments) <= 1 || min(lens) >= min_regime_size) {
        break
      }

      i_short <- which.min(lens)

      if (i_short == 1) {
        merge_i <- 1
      } else if (i_short == nrow(segments)) {
        merge_i <- nrow(segments) - 1
      } else {
        d_left <- segment_distance(z, segments, i_short - 1)
        d_right <- segment_distance(z, segments, i_short)

        merge_i <- if (d_left <= d_right) i_short - 1 else i_short
      }

      segments <- merge_pair(segments, merge_i)
    }

    # Merge similar adjacent regimes
    if (merge_similar_regimes && nrow(segments) > 1) {
      repeat {
        if (nrow(segments) <= 1) break

        d <- sapply(seq_len(nrow(segments) - 1), function(i) {
          segment_distance(z, segments, i)
        })

        if (length(d) == 0 || min(d) > merge_threshold) break

        merge_i <- which.min(d)
        segments <- merge_pair(segments, merge_i)
      }
    }

    # Enforce max number of regimes
    if (!is.null(max_regimes_per_dim)) {
      repeat {
        if (nrow(segments) <= max_regimes_per_dim || nrow(segments) <= 1) {
          break
        }

        d <- sapply(seq_len(nrow(segments) - 1), function(i) {
          segment_distance(z, segments, i)
        })

        merge_i <- which.min(d)
        segments <- merge_pair(segments, merge_i)
      }
    }

    segments
  }

  calc_regime_strength <- function(z, profile) {
    z <- as.numeric(z)

    total <- sum((z - mean(z))^2)

    if (!is.finite(total) || total <= .Machine$double.eps) {
      return(0)
    }

    mean_between <- sum(profile$n * (profile$mean - mean(z))^2)
    mean_strength <- mean_between / total

    global_sd <- sd(z)
    if (!is.finite(global_sd) || global_sd == 0) global_sd <- 1

    log_sd_dev <- log(profile$sd / global_sd)
    var_strength <- 1 - exp(-weighted.mean(log_sd_dev^2, w = profile$n))

    strength <- 0.75 * mean_strength + 0.25 * var_strength

    max(0, min(1, strength))
  }

  stable_changepoints <- function(z, base_cpts) {
    if (!stability_bootstrap || length(base_cpts) == 0) {
      return(list(
        cpts = base_cpts,
        stability = rep(NA_real_, length(base_cpts))
      ))
    }

    z <- as.numeric(z)

    global_sd <- sd(z)
    if (!is.finite(global_sd) || global_sd == 0) {
      global_sd <- 1
    }

    counts <- rep(0, length(base_cpts))
    valid_b <- 0

    for (b in seq_len(n_bootstrap)) {
      z_b <- z + rnorm(length(z), mean = 0, sd = global_sd * bootstrap_noise)

      fit_b <- run_cpt(z_b)

      if (is.null(fit_b)) next

      cpts_b <- changepoint::cpts(fit_b)
      cpts_b <- cpts_b[cpts_b < length(z)]

      valid_b <- valid_b + 1

      for (i in seq_along(base_cpts)) {
        if (any(abs(cpts_b - base_cpts[i]) <= cpt_tolerance)) {
          counts[i] <- counts[i] + 1
        }
      }
    }

    stability <- if (valid_b == 0) {
      rep(0, length(base_cpts))
    } else {
      counts / valid_b
    }

    keep <- stability >= min_cpt_stability

    list(
      cpts = base_cpts[keep],
      stability = stability
    )
  }

  fit_one_dimension <- function(z, dim_id) {
    z <- as.numeric(z)

    fit <- run_cpt(z)

    if (is.null(fit)) {
      raw_cpts <- integer(0)
    } else {
      raw_cpts <- changepoint::cpts(fit)
      raw_cpts <- raw_cpts[raw_cpts < length(z)]
    }

    stable <- stable_changepoints(z, raw_cpts)
    stable_cpts <- stable$cpts

    segments <- make_segments(stable_cpts, length(z))
    segments <- postprocess_segments(z, segments)

    regime <- segments_to_regime(segments, length(z))
    profile <- profile_from_segments(z, segments, dim_id = dim_id)

    final_cpts <- segments$end[-nrow(segments)]

    strength <- calc_regime_strength(z, profile)

    list(
      fit = fit,
      raw_cpts = raw_cpts,
      raw_cpt_stability = stable$stability,
      cpts = final_cpts,
      segments = segments,
      regime = regime,
      profile = profile,
      regime_strength = strength
    )
  }

  # ------------------------------------------------------------
  # 6. Fit changepoints on candidate dimensions
  # ------------------------------------------------------------

  cpt_results_all <- lapply(seq_along(initial_dims), function(i) {
    dim_id <- initial_dims[i]
    fit_one_dimension(scores_initial[, i], dim_id = dim_id)
  })

  names(cpt_results_all) <- paste0("dim_", initial_dims)

  diagnostics <- do.call(
    rbind,
    lapply(seq_along(initial_dims), function(i) {
      dim_id <- initial_dims[i]
      res <- cpt_results_all[[i]]
      prof <- res$profile

      data.frame(
        dim = dim_id,
        explained_variance = explained_var[dim_id],
        cumulative_explained_variance = cum_explained_var[dim_id],
        n_regimes = nrow(prof),
        n_changepoints = length(res$cpts),
        min_regime_length = min(prof$n),
        mean_regime_length = mean(prof$n),
        regime_strength = res$regime_strength,
        hybrid_score = explained_var[dim_id] * res$regime_strength,
        stringsAsFactors = FALSE
      )
    })
  )

  # ------------------------------------------------------------
  # 7. Final dimension selection
  # ------------------------------------------------------------

  if (reduction_criteria == "regime_strength") {
    selected_dims <- diagnostics$dim[
      diagnostics$regime_strength >= min_regime_strength
    ]

    if (length(selected_dims) == 0) {
      selected_dims <- diagnostics$dim[which.max(diagnostics$regime_strength)]
    }

  } else if (reduction_criteria == "hybrid") {
    k_var <- which(cum_explained_var >= var_threshold)[1]
    if (is.na(k_var)) k_var <- max_dim

    valid <- diagnostics$dim[
      diagnostics$regime_strength >= min_regime_strength &
        diagnostics$dim <= k_var
    ]

    if (length(valid) == 0) {
      valid <- diagnostics$dim[which.max(diagnostics$hybrid_score)]
    }

    selected_dims <- valid

  } else {
    selected_dims <- initial_dims

    if (filter_weak_dims) {
      valid <- diagnostics$dim[
        diagnostics$regime_strength >= min_regime_strength
      ]

      selected_dims <- intersect(selected_dims, valid)

      if (length(selected_dims) == 0) {
        selected_dims <- diagnostics$dim[which.max(diagnostics$regime_strength)]
      }
    }
  }

  if (!is.null(max_selected_dims) && length(selected_dims) > max_selected_dims) {
    score_order <- diagnostics$dim[
      order(diagnostics$hybrid_score, decreasing = TRUE)
    ]

    selected_dims <- score_order[score_order %in% selected_dims][seq_len(max_selected_dims)]
    selected_dims <- sort(selected_dims)
  }

  selected_dims <- sort(unique(selected_dims))

  diagnostics$selected <- diagnostics$dim %in% selected_dims

  selected_idx <- match(selected_dims, initial_dims)

  scores <- scores_initial[, selected_idx, drop = FALSE]
  colnames(scores) <- paste0("svd_dim_", selected_dims)

  cpt_results <- cpt_results_all[paste0("dim_", selected_dims)]

  regimes <- as.data.frame(
    setNames(
      lapply(seq_along(selected_dims), function(i) {
        cpt_results[[i]]$regime
      }),
      paste0(regime_prefix, "_dim_", selected_dims)
    )
  )

  profiles <- do.call(
    rbind,
    lapply(cpt_results, function(x) x$profile)
  )

  rownames(profiles) <- NULL

  V_selected <- V_all[, selected_dims, drop = FALSE]
  colnames(V_selected) <- paste0("svd_dim_", selected_dims)
  rownames(V_selected) <- series_cols

  # ------------------------------------------------------------
  # 8. New-data transform - FIXED
  # ------------------------------------------------------------

  return_probabilities_default <- return_probabilities

  transform_newdata <- function(
    newdata,
    return_scores = FALSE,
    return_probabilities = return_probabilities_default
  ) {
    newdata <- as.data.frame(newdata)

    missing_cols <- setdiff(series_cols, names(newdata))

    if (length(missing_cols) > 0) {
      stop("newdata is missing columns: ", paste(missing_cols, collapse = ", "))
    }

    X_new <- as.matrix(newdata[, series_cols, drop = FALSE])

    if (anyNA(X_new)) {
      if (na_action == "stop") {
        stop("Missing values found in newdata.")
      }

      for (j in seq_len(ncol(X_new))) {
        X_new[is.na(X_new[, j]), j] <- impute_values[j]
      }
    }

    X_new_scaled <- sweep(X_new, 2, center_vec, "-")
    X_new_scaled <- sweep(X_new_scaled, 2, scale_vec, "/")

    new_scores <- X_new_scaled %*% V_selected
    colnames(new_scores) <- paste0("svd_dim_", selected_dims)

    assign_one_dim <- function(z_new, profile_dim, dim_id) {
      regimes_available <- profile_dim$regime

      prob_mat <- matrix(
        NA_real_,
        nrow = length(z_new),
        ncol = length(regimes_available)
      )

      colnames(prob_mat) <- paste0("p_dim_", dim_id, "_R", regimes_available)

      assigned <- integer(length(z_new))

      for (i in seq_along(z_new)) {
        d <- abs(z_new[i] - profile_dim$mean) / profile_dim$sd

        w <- exp(
          -0.5 * d^2 / max(newdata_temperature, .Machine$double.eps)
        )

        if (sum(w) <= .Machine$double.eps || any(!is.finite(w))) {
          w <- rep(1 / length(w), length(w))
        } else {
          w <- w / sum(w)
        }

        prob_mat[i, ] <- w
        assigned[i] <- regimes_available[which.max(w)]
      }

      list(
        assigned = assigned,
        probabilities = as.data.frame(prob_mat)
      )
    }

    assigned_list <- vector("list", length(selected_dims))
    prob_list <- vector("list", length(selected_dims))

    for (i in seq_along(selected_dims)) {
      dim_id <- selected_dims[i]

      profile_j <- profiles[profiles$dim == dim_id, , drop = FALSE]

      out_j <- assign_one_dim(
        z_new = new_scores[, i],
        profile_dim = profile_j,
        dim_id = dim_id
      )

      assigned_list[[i]] <- out_j$assigned
      prob_list[[i]] <- out_j$probabilities
    }

    new_regimes <- as.data.frame(
      setNames(
        assigned_list,
        paste0(regime_prefix, "_dim_", selected_dims)
      )
    )

    out <- list(regimes = new_regimes)

    if (return_scores) {
      out$scores <- as.data.frame(new_scores)
    }

    if (return_probabilities) {
      out$probabilities <- do.call(cbind, prob_list)
    }

    if (length(out) == 1) {
      return(out$regimes)
    }

    out
  }

  # ------------------------------------------------------------
  # 9. Plotting
  # ------------------------------------------------------------

  plot_regimes <- function(
    dims = NULL,
    type = c("scores", "heatmap", "summary"),
    show_changepoints = TRUE,
    label_regimes = TRUE,
    main = "SVD Changepoint Regimes",
    line_lwd = 1.1,
    point_cex = 0.45,
    alpha = 0.22,
    regime_label_cex = 0.75
  ) {
    type <- match.arg(type)

    if (is.null(dims)) {
      dims <- selected_dims
    }

    dims <- intersect(dims, selected_dims)

    if (length(dims) == 0) {
      stop("No valid selected dimensions to plot.")
    }

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    if (type == "scores") {
      par(
        mfrow = c(length(dims), 1),
        mar = c(3.5, 4.2, 2.5, 1),
        oma = c(0, 0, 2, 0)
      )

      for (dim_id in dims) {
        score_col <- paste0("svd_dim_", dim_id)
        z <- as.numeric(scores[, score_col])

        res <- cpt_results[[paste0("dim_", dim_id)]]

        r <- res$regime
        cpts <- res$cpts

        reg_ids <- sort(unique(r))

        cols <- grDevices::hcl.colors(
          n = max(length(reg_ids), 3),
          palette = "Pastel 1"
        )[seq_along(reg_ids)]

        names(cols) <- as.character(reg_ids)

        y_range <- range(z, na.rm = TRUE)
        y_pad <- diff(y_range) * 0.08

        if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

        plot(
          z,
          type = "n",
          xlab = "Time index",
          ylab = paste0("SVD dim ", dim_id),
          main = paste0(
            "Dimension ", dim_id,
            " | regimes: ", length(reg_ids),
            " | strength: ", round(res$regime_strength, 3)
          ),
          ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)
        )

        usr <- par("usr")

        run <- rle(r)
        ends <- cumsum(run$lengths)
        starts <- c(1, head(ends, -1) + 1)

        for (i in seq_along(run$values)) {
          g <- as.character(run$values[i])

          rect(
            xleft = starts[i],
            ybottom = usr[3],
            xright = ends[i],
            ytop = usr[4],
            col = grDevices::adjustcolor(cols[g], alpha.f = alpha),
            border = NA
          )

          if (label_regimes) {
            text(
              x = mean(c(starts[i], ends[i])),
              y = usr[4] - 0.06 * diff(usr[3:4]),
              labels = paste0("R", run$values[i]),
              cex = regime_label_cex,
              col = "gray25"
            )
          }
        }

        grid(col = "gray88", lty = 1)

        lines(z, lwd = line_lwd, col = "black")
        points(z, pch = 16, cex = point_cex, col = "black")

        if (show_changepoints && length(cpts) > 0) {
          abline(v = cpts, lty = 2, lwd = 1, col = "gray25")
        }

        box()
      }

      mtext(main, outer = TRUE, cex = 1.1, font = 2)
    }

    if (type == "heatmap") {
      reg_mat <- as.matrix(
        regimes[, paste0(regime_prefix, "_dim_", dims), drop = FALSE]
      )

      max_r <- max(reg_mat, na.rm = TRUE)

      cols <- grDevices::hcl.colors(
        n = max(max_r, 3),
        palette = "Pastel 1"
      )

      par(mar = c(4, 5, 3, 1))

      image(
        x = seq_len(nrow(reg_mat)),
        y = seq_along(dims),
        z = reg_mat,
        col = cols,
        axes = FALSE,
        xlab = "Time index",
        ylab = "Selected SVD dimension",
        main = "Regime Heatmap"
      )

      axis(1)
      axis(2, at = seq_along(dims), labels = paste0("Dim ", dims), las = 2)
      box()
    }

    if (type == "summary") {
      d <- diagnostics[diagnostics$dim %in% selected_dims, , drop = FALSE]

      par(mfrow = c(3, 1), mar = c(4, 4.2, 2.5, 1))

      barplot(
        d$n_regimes,
        names.arg = paste0("D", d$dim),
        main = "Number of regimes by selected dimension",
        ylab = "N regimes"
      )

      barplot(
        d$regime_strength,
        names.arg = paste0("D", d$dim),
        main = "Regime strength by selected dimension",
        ylab = "Strength"
      )

      barplot(
        d$min_regime_length,
        names.arg = paste0("D", d$dim),
        main = "Minimum regime length by selected dimension",
        ylab = "Min length"
      )
    }

    invisible(NULL)
  }

  # ------------------------------------------------------------
  # 10. Output
  # ------------------------------------------------------------

  out <- list(
    regimes = regimes,
    selected_dims = selected_dims,
    n_dim = length(selected_dims),
    initial_dims = initial_dims,
    reduction_criteria = reduction_criteria,

    explained_variance = explained_var,
    cumulative_explained_variance = cum_explained_var,
    selected_cumulative_variance = max(cum_explained_var[selected_dims]),

    diagnostics = diagnostics,
    profiles = profiles,
    changepoint_results = cpt_results,
    denoiser = list(
      enabled = isTRUE(denoise),
      original_series_cols = original_series_cols,
      selected_features = series_cols,
      removed_features = setdiff(original_series_cols, series_cols),
      report = denoiser_report,
      settings = denoiser_settings
    ),

    svd = list(
      V = V_selected,
      selected_dims = selected_dims,
      singular_values = singular_values,
      center = center_vec,
      scale = scale_vec,
      impute_values = impute_values,
      series_cols = series_cols
    ),

    transform_newdata = transform_newdata,
    plot_regimes = plot_regimes
  )

  if (return_scores) {
    out$scores <- as.data.frame(scores)
  }

  class(out) <- "svd_changepoint_regimes_v2"

  if (plot) {
    plot_regimes(type = plot_type)
  }

  out
}
