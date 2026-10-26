.craft_autotune_integer <- function(x, name, lower = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > .Machine$integer.max || x != floor(x)) {
    stop("`", name, "` must be a finite integer >= ", lower, ".")
  }
}

.craft_autotune_named_list <- function(x, name) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))) {
    stop("`", name, "` must be a list with unique, non-empty names.")
  }
}

.craft_autotune_validate_spec <- function(spec) {
  if (is.function(spec)) return(invisible(NULL))
  if (is.list(spec) && !is.null(spec$type)) {
    type <- match.arg(spec$type, c("choice", "integer", "numeric"))
    if (type == "choice") {
      if (!is.atomic(spec$values) || !length(spec$values) || anyNA(spec$values)) {
        stop("Choice specifications require non-empty atomic `values` without missing values.")
      }
    } else {
      for (nm in c("min", "max")) {
        x <- spec[[nm]]
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
          stop("Range specifications require finite scalar `min` and `max`.")
        }
      }
      if (spec$min > spec$max) stop("Search range `min` cannot exceed `max`.")
      if (type == "integer" && (spec$min != floor(spec$min) ||
          spec$max != floor(spec$max) || spec$min < -.Machine$integer.max ||
          spec$max > .Machine$integer.max)) stop("Integer ranges require integer bounds.")
      if (!is.null(spec$log) && (!is.logical(spec$log) ||
          length(spec$log) != 1L || is.na(spec$log))) stop("`log` must be TRUE or FALSE.")
      if (isTRUE(spec$log) && spec$min <= 0) stop("Log ranges require positive bounds.")
    }
  } else if (is.atomic(spec) && (!length(spec) || anyNA(spec))) {
    stop("Atomic search specifications must be non-empty and cannot contain missing values.")
  }
  invisible(NULL)
}

# ============================================================
# Random search wrapper for craft_fit()
# ============================================================

.craft_autotune_default_search_space <- function() {
  list(
    window = function() sample(c(10L, 15L, 20L, 30L), 1L),
    n_factors = function() sample(2L:5L, 1L),
    var_threshold = function() stats::runif(1L, 0.85, 0.98),
    factor_selection = function() sample(c("fixed", "variance"), 1L),
    max_regimes_per_factor = function() sample(3L:6L, 1L),
    min_segment = function() sample(c(10L, 15L, 20L, 30L, 40L), 1L),
    transition_alpha = function() exp(stats::runif(1L, log(0.25), log(5))),
    transition_min_count = function() sample(1L:5L, 1L),
    sampler_weight = function() sample(c("uniform", "weighted_bal_acc", "bal_acc"), 1L),
    tail_penalty_lambda = function() sample(c(0, exp(stats::runif(1L, log(0.02), log(0.75)))), 1L),
    tail_penalty_threshold = function() stats::runif(1L, 2.5, 4.5),
    tail_penalty_power = function() sample(c(1, 1.5, 2), 1L),
    tail_penalty_aggregation = function() sample(c("mean", "max"), 1L)
  )
}

.craft_autotune_draw_value <- function(spec) {
  if (is.function(spec)) {
    return(spec())
  }

  if (is.list(spec) && !is.null(spec$type)) {
    type <- match.arg(spec$type, c("choice", "integer", "numeric"))

    if (type == "choice") {
      return(spec$values[sample.int(length(spec$values), 1L)])
    }

    if (type == "integer") {
      return(seq.int(spec$min, spec$max)[sample.int(spec$max - spec$min + 1L, 1L)])
    }

    if (type == "numeric") {
      lo <- as.numeric(spec$min)
      hi <- as.numeric(spec$max)

      if (isTRUE(spec$log)) {
        return(exp(stats::runif(1L, log(lo), log(hi))))
      }

      return(stats::runif(1L, lo, hi))
    }
  }

  if (is.atomic(spec) && length(spec) > 1L) {
    return(spec[sample.int(length(spec), 1L)])
  }

  spec
}

.craft_autotune_draw_params <- function(search_space) {
  out <- lapply(search_space, .craft_autotune_draw_value)
  names(out) <- names(search_space)
  out
}

.craft_autotune_param_row <- function(params) {
  if (!length(params)) return(data.frame(row.names = 1L))
  out <- lapply(params, function(x) {
    if (is.atomic(x) && length(x) == 1L) {
      return(x)
    }

    paste(utils::capture.output(utils::str(x, give.attr = FALSE)), collapse = " ")
  })

  as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

.craft_autotune_score_fit <- function(fit, score_metric = "return_rmse") {
  if (is.null(fit$backtest) || is.null(fit$backtest$summary$overall) || !nrow(fit$backtest$summary$overall)) {
    stop("Backtest summary is empty; set `n_testing` high enough and check fit failures.")
  }

  overall <- fit$backtest$summary$overall[1L, , drop = FALSE]

  score <- switch(
    score_metric,
    return_rmse = overall$return_rmse,
    return_mae = overall$return_mae,
    level_rmse = overall$level_rmse,
    level_mae = overall$level_mae,
    coverage_90_gap = abs(overall$coverage_90 - 0.90),
    direction_accuracy = -overall$direction_accuracy,
    composite = overall$return_rmse,
    stop("Unknown `score_metric`.")
  )

  list(score = as.numeric(score), overall = overall)
}

.craft_autotune_evaluate_one <- function(task) {
  tryCatch({
    fit_args <- utils::modifyList(task$base_args, task$params, keep.null = TRUE)
    fit_args$seed <- as.integer((as.double(task$seed) + task$model_id) %% .Machine$integer.max)
    fit <- do.call(craft_fit, fit_args)
    scored <- .craft_autotune_score_fit(fit, task$score_metric)
    if (length(scored$score) != 1L || !is.finite(scored$score)) {
      stop("The selected backtest score is not finite.")
    }
    if (nrow(fit$backtest$failures)) stop("One or more backtest origins failed.")
    param_row <- .craft_autotune_param_row(task$params)

    result <- cbind(
      data.frame(
        model_id = task$model_id,
        ok = TRUE,
        score = scored$score,
        error = NA_character_,
        stringsAsFactors = FALSE
      ),
      scored$overall,
      param_row
    )

    list(
      ok = TRUE,
      model_id = task$model_id,
      params = task$params,
      result = result,
      fit = fit,
      error = NA_character_
    )
  }, error = function(e) {
    param_row <- .craft_autotune_param_row(task$params)
    result <- cbind(
      data.frame(
        model_id = task$model_id,
        ok = FALSE,
        score = Inf,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      ),
      param_row
    )

    list(
      ok = FALSE,
      model_id = task$model_id,
      params = task$params,
      result = result,
      fit = NULL,
      error = conditionMessage(e)
    )
  })
}

.craft_autotune_bind_rows <- function(x) {
  x <- Filter(Negate(is.null), x)

  if (!length(x)) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(row) {
    missing <- setdiff(all_names, names(row))
    for (nm in missing) row[[nm]] <- NA
    row[, all_names, drop = FALSE]
  })

  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.craft_autotune_rank01 <- function(x) {
  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)

  if (!any(finite)) {
    return(out)
  }

  if (sum(finite) == 1L || diff(range(x[finite], na.rm = TRUE)) <= .Machine$double.eps) {
    out[finite] <- 0
    return(out)
  }

  out[finite] <- (rank(x[finite], ties.method = "average") - 1) / (sum(finite) - 1)
  out
}

.craft_autotune_composite_scores <- function(result_table, weights = NULL) {
  scores <- rep(Inf, nrow(result_table))

  if (!nrow(result_table)) {
    return(scores)
  }

  get_col <- function(name) {
    if (name %in% names(result_table)) as.numeric(result_table[[name]]) else rep(NA_real_, nrow(result_table))
  }

  penalties <- list(
    return_rmse = get_col("return_rmse"),
    return_mae = get_col("return_mae"),
    level_rmse = get_col("level_rmse"),
    level_mae = get_col("level_mae"),
    coverage_90_gap = abs(get_col("coverage_90") - 0.90),
    direction_error = 1 - get_col("direction_accuracy"),
    mean_p_actual_gap = abs(get_col("mean_p_actual") - 0.50)
  )

  all_penalty_names <- names(penalties)
  user_weights <- NULL

  if (!is.null(weights)) {
    if (is.null(names(weights))) {
      stop("`composite_weights` must be a named numeric vector or named list.")
    }

    user_weights <- unlist(weights, use.names = TRUE)
    unknown <- setdiff(names(user_weights), all_penalty_names)
    if (length(unknown)) {
      stop("Unknown `composite_weights`: ", paste(unknown, collapse = ", "))
    }
  }

  usable <- vapply(penalties, function(x) any(is.finite(x)), logical(1))
  penalties <- penalties[usable]

  if (!length(penalties)) {
    return(scores)
  }

  weights <- rep(1, length(penalties))
  names(weights) <- names(penalties)

  if (!is.null(user_weights)) {
    common <- intersect(names(user_weights), names(weights))
    weights[common] <- as.numeric(user_weights[common])
  }

  weights[!is.finite(weights) | weights < 0] <- 0

  if (!any(weights > 0)) {
    stop("At least one `composite_weights` value must be positive.")
  }

  component_scores <- as.data.frame(lapply(penalties, .craft_autotune_rank01), stringsAsFactors = FALSE)
  ok <- rep(FALSE, nrow(result_table))
  if ("ok" %in% names(result_table)) {
    ok <- !is.na(result_table$ok) & result_table$ok
  }

  weighted <- sweep(as.matrix(component_scores), 2L, weights, `*`)
  available_weight <- sweep(!is.na(as.matrix(component_scores)), 2L, weights, `*`)
  denom <- rowSums(available_weight, na.rm = TRUE)
  valid <- ok & denom > 0
  scores[valid] <- rowSums(weighted[valid, , drop = FALSE], na.rm = TRUE) / denom[valid]
  scores
}

#' Tune CRAFT hyperparameters by random search
#'
#' See the craft_autotune manual page for search specifications, scoring,
#' parallel execution, and the returned candidate diagnostics.
craft_autotune <- function(
  series,
  volumes = NULL,
  n_models = 30L,
  top_k = 5L,
  n_testing = 5L,
  search_space = NULL,
  score_metric = c(
    "return_rmse",
    "return_mae",
    "level_rmse",
    "level_mae",
    "coverage_90_gap",
    "direction_accuracy",
    "composite"
  ),
  seed = 42,
  fixed_args = list(),
  composite_weights = NULL,
  keep_all_fits = FALSE,
  search_parallel = FALSE,
  search_workers = NULL,
  verbose = TRUE,
  ...
) {
  score_metric <- match.arg(score_metric)
  .craft_autotune_integer(n_models, "n_models")
  .craft_autotune_integer(top_k, "top_k")
  .craft_autotune_integer(n_testing, "n_testing")
  .craft_autotune_integer(seed, "seed", lower = 0)
  for (nm in c("search_parallel", "keep_all_fits", "verbose")) {
    x <- get(nm)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop("`", nm, "` must be TRUE or FALSE.")
    }
  }
  if (!is.null(search_workers)) .craft_autotune_integer(search_workers, "search_workers")
  if (is.null(search_space)) search_space <- .craft_autotune_default_search_space()
  .craft_autotune_named_list(search_space, "search_space")
  .craft_autotune_named_list(fixed_args, "fixed_args")
  dots <- list(...)
  .craft_autotune_named_list(dots, "...")
  allowed <- setdiff(names(formals(craft_fit)), c("series", "volumes", "seed", "n_testing"))
  for (nm in c("search_space", "fixed_args", "dots")) {
    unknown <- setdiff(names(get(nm)), allowed)
    if (length(unknown)) stop("Unknown or reserved craft_fit arguments in `", nm,
                              "`: ", paste(unknown, collapse = ", "))
  }
  if (search_parallel && "backtest_parallel" %in% names(search_space)) {
    stop("Cannot tune `backtest_parallel` during parallel search.")
  }
  invisible(lapply(search_space, .craft_autotune_validate_spec))
  if (!is.null(composite_weights)) {
    if (is.list(composite_weights)) composite_weights <- unlist(composite_weights)
    valid_names <- c("return_rmse", "return_mae", "level_rmse", "level_mae",
                     "coverage_90_gap", "direction_error", "mean_p_actual_gap")
    if (!is.numeric(composite_weights) || !length(composite_weights) ||
        is.null(names(composite_weights)) || anyNA(names(composite_weights)) ||
        anyDuplicated(names(composite_weights)) ||
        any(!names(composite_weights) %in% valid_names) ||
        any(!is.finite(composite_weights) | composite_weights < 0) ||
        !any(composite_weights > 0)) {
      stop("`composite_weights` must have unique known names and finite non-negative weights, with at least one positive weight.")
    }
  }

  set.seed(seed)
  sampled_params <- lapply(seq_len(n_models), function(i) .craft_autotune_draw_params(search_space))

  base_args <- utils::modifyList(dots, fixed_args, keep.null = TRUE)
  base_args$series <- as.data.frame(series, stringsAsFactors = FALSE)
  base_args$volumes <- volumes
  base_args$n_testing <- n_testing

  if (is.null(base_args$backtest_parallel)) {
    base_args$backtest_parallel <- FALSE
  }

  if (is.null(base_args$return_train_probs)) {
    base_args$return_train_probs <- FALSE
  }

  if (is.null(base_args$verbose)) {
    base_args$verbose <- FALSE
  }

  if (search_parallel) base_args$backtest_parallel <- FALSE

  tasks <- lapply(seq_len(n_models), function(i) {
    list(
      model_id = i,
      params = sampled_params[[i]],
      base_args = base_args,
      score_metric = score_metric,
      seed = seed
    )
  })

  if (isTRUE(search_parallel) && n_models > 1L) {
    workers <- if (is.null(search_workers)) 2L else search_workers
    workers <- min(workers, n_models)

    if (verbose) message("Running random search in parallel with ", workers, " workers.")

    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    results <- parallel::parLapply(cl, tasks, .craft_autotune_evaluate_one)
  } else {
    if (verbose) message("Running random search serially over ", n_models, " models.")
    results <- lapply(tasks, .craft_autotune_evaluate_one)
  }

  result_table <- .craft_autotune_bind_rows(lapply(results, `[[`, "result"))

  if (score_metric == "composite") {
    result_table$score <- .craft_autotune_composite_scores(result_table, composite_weights)
  }

  score_by_model_id <- stats::setNames(result_table$score, result_table$model_id)
  ok_idx <- which(vapply(results, function(x) {
    score <- score_by_model_id[as.character(x$model_id)]
    isTRUE(x$ok) && length(score) == 1L && is.finite(score)
  }, logical(1)))
  ok_idx <- ok_idx[order(vapply(results[ok_idx], function(x) score_by_model_id[as.character(x$model_id)], numeric(1)))]
  top_idx <- head(ok_idx, min(top_k, length(ok_idx)))

  top_models <- lapply(top_idx, function(i) {
    list(
      model_id = results[[i]]$model_id,
      params = results[[i]]$params,
      score = unname(score_by_model_id[as.character(results[[i]]$model_id)]),
      fit = results[[i]]$fit
    )
  })

  all_fits <- if (isTRUE(keep_all_fits)) lapply(results, `[[`, "fit") else NULL

  out <- list(
    call = match.call(),
    score_metric = score_metric,
    search_space = search_space,
    results = result_table,
    top_models = top_models,
    failures = result_table[!result_table$ok, , drop = FALSE],
    all_fits = all_fits,
    settings = list(
      n_models = n_models,
      top_k = top_k,
      n_testing = n_testing,
      seed = seed,
      composite_weights = composite_weights,
      search_parallel = search_parallel,
      search_workers = search_workers,
      keep_all_fits = keep_all_fits
    )
  )

  if (!length(top_models)) warning("No candidate produced a finite score; inspect $failures.")
  class(out) <- c("craft_autotune", class(out))
  out
}



