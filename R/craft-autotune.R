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


# Caches belong to a single search and a single immutable input history.
# Prefix length is therefore sufficient to identify a training data set.
.craft_fast_cache <- function(size) {
  x <- new.env(parent = emptyenv())
  x$values <- new.env(parent = emptyenv())
  x$keys <- list()
  x$next_id <- 0
  x$order <- character()
  x$size <- size
  x$hits <- 0L
  x$misses <- 0L
  x
}

.craft_fast_cached <- function(cache, key, compute) {
  if (is.null(cache)) return(compute())
  serialized <- serialize(key, NULL, version = 2)
  matched <- which(vapply(cache$keys, identical, logical(1), serialized))
  if (length(matched)) {
    key <- names(cache$keys)[matched[1L]]
    cache$hits <- cache$hits + 1L
    cache$order <- c(setdiff(cache$order, key), key)
    return(get(key, envir = cache$values, inherits = FALSE))
  }
  cache$misses <- cache$misses + 1L
  value <- compute()
  if (length(cache$order) >= cache$size) {
    rm(list = cache$order[1L], envir = cache$values)
    cache$keys[[cache$order[1L]]] <- NULL
    cache$order <- cache$order[-1L]
  }
  cache$next_id <- cache$next_id + 1
  key <- as.character(cache$next_id)
  cache$keys[[key]] <- serialized
  assign(key, value, envir = cache$values)
  cache$order <- c(cache$order, key)
  value
}

.craft_fast_regimes <- function(cache, key, args) {
  # Stochastic/bootstrap and plotting calls must execute on every invocation.
  if (isTRUE(args$stability_bootstrap) || !is.null(args$seed) || isTRUE(args$plot)) {
    return(do.call(svd_changepoint_regimes, args))
  }
  settings <- args[setdiff(names(args), "X")]
  .craft_fast_cached(cache, c(key, list(settings)),
                     function() do.call(svd_changepoint_regimes, args))
}

.craft_fast_worker_state <- new.env(parent = emptyenv())

.craft_fast_worker_init <- function(context) {
  state <- new.env(parent = emptyenv())
  state$context <- context
  state$cache <- if (context$cache) .craft_fast_cache(context$cache_size) else NULL
  .craft_fast_worker_state$state <- state
  invisible(NULL)
}

.craft_fast_worker <- function(task, state = NULL) {
  if (is.null(state)) state <- .craft_fast_worker_state$state
  context <- state$context
  before <- if (is.null(state$cache)) c(0L, 0L) else c(state$cache$hits, state$cache$misses)
  started <- proc.time()[[3L]]
  origins <- lapply(task$points, function(point) {
    args <- task$args
    args$seed <- task$seed
    # The point ID is stable across stages and scheduling assignments.
    set.seed(as.integer((as.double(task$seed) + point) %% .Machine$integer.max))
    out <- .v5_backtest_one_origin(list(
      series = context$series, volumes = context$volumes,
      origin = context$origins[point], point_id = point,
      horizons = .v5_backtest_horizons(args$window, args$horizon),
      fit_args = args, score_only = TRUE, cache = state$cache
    ))
    if (out$ok) {
      # No regime models, sampled trajectories, distributions, or closures cross
      # the worker boundary. These columns suffice for all standard summaries.
      keep <- c("point_id", "origin_index", "asset", "horizon",
                "return_abs_error", "return_sq_error", "level_abs_error",
                "level_sq_error", "covered_90", "direction_correct", "p_actual")
      out$predictions <- out$predictions[, keep, drop = FALSE]
    }
    out
  })
  after <- if (is.null(state$cache)) c(0L, 0L) else c(state$cache$hits, state$cache$misses)
  list(model_id = task$model_id, points = task$points, origins = origins,
       performance = data.frame(model_id = task$model_id,
         origins_computed = length(task$points),
         elapsed = proc.time()[[3L]] - started,
         cache_hits = after[1L] - before[1L],
         cache_misses = after[2L] - before[2L]))
}

.craft_fast_backtest <- function(point_results, origins, horizons) {
  predictions <- .v5_bind_rows(lapply(point_results, function(x) {
    if (isTRUE(x$ok)) x$predictions else NULL
  }))
  failures <- .v5_bind_rows(lapply(point_results, function(x) {
    if (isTRUE(x$ok)) return(NULL)
    data.frame(point_id = x$point_id, origin_index = x$origin_index, error = x$error)
  }))
  list(settings = list(actual_n_testing = length(point_results),
         origins = origins, horizons = horizons, scoring_only = TRUE),
       predictions = predictions, failures = failures,
       summary = .v5_backtest_summary(predictions))
}

.craft_fast_validate <- function(n_models, top_k, n_testing, seed, search_workers,
                                 cache_size, survivor_fraction, flags,
                                 search_space, fixed_args, dots, weights) {
  for (nm in c("n_models", "top_k", "n_testing", "cache_size")) {
    .craft_autotune_integer(get(nm), nm)
  }
  .craft_autotune_integer(seed, "seed", lower = 0)
  if (!is.null(search_workers)) .craft_autotune_integer(search_workers, "search_workers")
  for (nm in names(flags)) {
    x <- flags[[nm]]
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("`", nm, "` must be TRUE or FALSE.")
  }
  if (!is.numeric(survivor_fraction) || length(survivor_fraction) != 1L ||
      !is.finite(survivor_fraction) || survivor_fraction <= 0 || survivor_fraction > 1) {
    stop("`survivor_fraction` must be in (0, 1].")
  }
  allowed <- setdiff(names(formals(craft_fit)), c("series", "volumes", "seed", "n_testing"))
  for (nm in c("search_space", "fixed_args", "dots")) {
    x <- get(nm)
    .craft_autotune_named_list(x, nm)
    unknown <- setdiff(names(x), allowed)
    if (length(unknown)) stop("Unknown or reserved craft_fit arguments in `", nm,
                              "`: ", paste(unknown, collapse = ", "))
  }
  if (any(c("backtest_parallel", "backtest_workers", "verbose", "return_train_probs") %in%
          names(search_space))) stop("Execution/diagnostic controls cannot be tuned in fast search.")
  invisible(lapply(search_space, .craft_autotune_validate_spec))
  if (!is.null(weights)) {
    weights <- unlist(weights)
    known <- c("return_rmse", "return_mae", "level_rmse", "level_mae",
               "coverage_90_gap", "direction_error", "mean_p_actual_gap")
    if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) ||
        anyNA(names(weights)) || anyDuplicated(names(weights)) ||
        any(!names(weights) %in% known) || any(!is.finite(weights) | weights < 0) ||
        !any(weights > 0)) stop("Invalid `composite_weights`.")
  }
}

#' Fast staged random-search tuning for CRAFT
#'
#' See the craft_autotune manual page for stage budgets, cache scope,
#' compact scoring results, and deferred fitting of selected candidates.
craft_autotune <- function(
  series, volumes = NULL, n_models = 30L, top_k = 5L, n_testing = 5L,
  search_space = NULL,
  score_metric = c("return_rmse", "return_mae", "level_rmse", "level_mae",
                   "coverage_90_gap", "direction_accuracy", "composite"),
  seed = 42, fixed_args = list(), composite_weights = NULL,
  keep_all_fits = FALSE, search_parallel = FALSE, search_workers = NULL,
  stage_budgets = NULL, survivor_fraction = 0.3,
  cache = TRUE, cache_size = 64L, verbose = TRUE, ...
) {
  score_metric <- match.arg(score_metric)
  if (is.null(search_space)) search_space <- .craft_autotune_default_search_space()
  dots <- list(...)
  .craft_fast_validate(n_models, top_k, n_testing, seed, search_workers,
    cache_size, survivor_fraction,
    list(keep_all_fits = keep_all_fits, search_parallel = search_parallel,
         cache = cache, verbose = verbose), search_space, fixed_args, dots,
    composite_weights)
  if (is.null(stage_budgets)) stage_budgets <- unique(c(min(2L, n_testing), n_testing))
  if (!is.numeric(stage_budgets) || !length(stage_budgets) ||
      any(!is.finite(stage_budgets)) || any(stage_budgets < 1) ||
      any(stage_budgets != floor(stage_budgets)) || any(diff(stage_budgets) <= 0) ||
      tail(stage_budgets, 1L) != n_testing) {
    stop("`stage_budgets` must be increasing positive integers ending at `n_testing`.")
  }
  if (!is.matrix(series) && !is.data.frame(series)) stop("`series` must be a matrix or data frame.")
  series <- as.data.frame(series, stringsAsFactors = FALSE)
  if (!nrow(series) || !ncol(series) || !all(vapply(series, is.numeric, logical(1)))) {
    stop("`series` must contain numeric time-series columns and observations.")
  }
  if (!is.null(volumes)) {
    if ((!is.matrix(volumes) && !is.data.frame(volumes)) ||
        !identical(dim(volumes), dim(series))) stop("`volumes` must match the dimensions of `series`.")
    volumes <- as.data.frame(volumes)
  }
  default_exprs <- formals(craft_fit)
  default_exprs$series <- NULL
  defaults <- lapply(default_exprs, eval, envir = environment(craft_fit))
  base_args <- utils::modifyList(defaults, utils::modifyList(dots, fixed_args,
                                keep.null = TRUE), keep.null = TRUE)
  base_args[c("series", "volumes")] <- NULL
  base_args$n_testing <- 0L
  base_args$backtest_parallel <- FALSE
  base_args$backtest_workers <- 1L
  base_args$return_train_probs <- FALSE
  base_args$verbose <- FALSE
  # Resolve match.arg defaults before serializing them as explicit arguments.
  for (nm in c("valid_split", "sampler_weight", "tail_penalty_aggregation", "tail_penalty_scope")) {
    if (!nm %in% c(names(dots), names(fixed_args))) base_args[[nm]] <- base_args[[nm]][1L]
  }
  set.seed(seed)
  sampled_params <- lapply(seq_len(n_models), function(i) .craft_autotune_draw_params(search_space))
  args <- lapply(sampled_params, function(p) utils::modifyList(base_args, p, keep.null = TRUE))
  seeds <- as.integer((as.double(seed) + seq_len(n_models)) %% .Machine$integer.max)
  errors <- rep(NA_character_, n_models)
  minimum <- maximum <- rep(NA_real_, n_models)
  for (i in seq_len(n_models)) {
    tryCatch({
      a <- args[[i]]
      .craft_autotune_integer(a$window, "window")
      .craft_autotune_integer(a$min_segment, "min_segment")
      .craft_autotune_integer(a$n_draws, "n_draws")
      if (a$window >= nrow(series)) stop("`window` must be shorter than the series.")
      if (!is.null(a$horizon)) {
        .craft_autotune_integer(a$horizon, "horizon")
        if (a$horizon > a$window) stop("`horizon` cannot exceed `window`.")
      }
      if (!is.null(a$backtest_min_train)) .craft_autotune_integer(a$backtest_min_train, "backtest_min_train")
      lo <- if (is.null(a$backtest_min_train)) {
        max(a$window + 2L, a$window + 2L * a$min_segment, 3L * a$window)
      } else max(a$window + 1L, a$backtest_min_train)
      hi <- nrow(series) - max(.v5_backtest_horizons(a$window, a$horizon))
      if (lo > hi) stop("No feasible backtest origins for this candidate.")
      minimum[i] <- lo
      maximum[i] <- hi
    }, error = function(e) { errors[i] <<- conditionMessage(e) })
  }
  valid <- which(is.na(errors))
  origins <- integer()
  if (length(valid)) {
    lo <- max(minimum[valid])
    hi <- min(maximum[valid])
    if (lo > hi) stop("Candidates have no common backtest origins; narrow the search space.")
    origins <- if (n_testing == 1L) as.integer(hi) else {
      unique(as.integer(round(seq(lo, hi, length.out = n_testing))))
    }
  }
  # First stage spans the full origin grid; later stages add rather than refit.
  priority <- if (length(origins)) unique(c(1L, length(origins), seq_along(origins))) else integer()
  budgets <- if (length(origins)) unique(pmin(stage_budgets, length(origins))) else integer()
  context <- list(series = series, volumes = volumes, origins = origins,
                  cache = cache, cache_size = cache_size)
  state <- new.env(parent = emptyenv())
  state$context <- context
  state$cache <- if (cache) .craft_fast_cache(cache_size) else NULL
  cl <- NULL
  workers <- 1L
  if (search_parallel && length(valid) > 1L) {
    workers <- min(if (is.null(search_workers)) 2L else search_workers, length(valid))
    cl <- parallel::makeCluster(workers)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
    # Export the immutable history once, not once per candidate or stage.
    parallel::clusterCall(cl, .craft_fast_worker_init, context)
  }
  point_results <- lapply(seq_len(n_models), function(i) vector("list", length(origins)))
  backtests <- vector("list", n_models)
  rows <- lapply(seq_len(n_models), function(i) {
    cbind(data.frame(model_id = i, ok = FALSE, score = Inf,
      error = errors[i], stage = 0L, n_testing = 0L,
      status = if (is.na(errors[i])) "pending" else "failed", selected = FALSE),
      .craft_autotune_param_row(sampled_params[[i]]))
  })
  stages <- performance <- list()
  active <- valid
  for (stage in seq_along(budgets)) {
    if (!length(active)) break
    points <- sort(priority[seq_len(budgets[stage])])
    if (verbose) message("Fast tuning stage ", stage, ": ", length(active),
                         " candidates, ", length(points), " origins each.")
    tasks <- lapply(active, function(i) list(model_id = i, args = args[[i]],
      seed = seeds[i], points = points[vapply(point_results[[i]][points], is.null, logical(1))]))
    evaluated <- if (is.null(cl)) lapply(tasks, .craft_fast_worker, state = state) else {
      parallel::parLapplyLB(cl, tasks, .craft_fast_worker, chunk.size = 1L)
    }
    for (item in evaluated) {
      i <- item$model_id
      point_results[[i]][item$points] <- item$origins
      perf <- item$performance
      perf$stage <- stage
      performance[[length(performance) + 1L]] <- perf
      bt <- .craft_fast_backtest(point_results[[i]][points], origins[points],
                                .v5_backtest_horizons(args[[i]]$window, args[[i]]$horizon))
      backtests[[i]] <- bt
      error <- tryCatch({
        if (nrow(bt$failures)) stop(paste(unique(bt$failures$error), collapse = "; "))
        scored <- .craft_autotune_score_fit(list(backtest = bt), score_metric)
        if (length(scored$score) != 1L || !is.finite(scored$score)) stop("Non-finite backtest score.")
        rows[[i]] <- cbind(data.frame(model_id = i, ok = TRUE,
          score = scored$score, error = NA_character_, stage = stage,
          n_testing = length(points), status = "evaluated", selected = FALSE),
          scored$overall, .craft_autotune_param_row(sampled_params[[i]]))
        NA_character_
      }, error = function(e) conditionMessage(e))
      if (!is.na(error)) {
        rows[[i]]$ok <- FALSE
        rows[[i]]$error <- error
        rows[[i]]$score <- Inf
        rows[[i]]$stage <- stage
        rows[[i]]$n_testing <- length(points)
        rows[[i]]$status <- "failed"
      }
    }
    table <- .craft_autotune_bind_rows(rows[active])
    if (score_metric == "composite") {
      table$score <- .craft_autotune_composite_scores(table, composite_weights)
      for (j in seq_len(nrow(table))) rows[[table$model_id[j]]]$score <- table$score[j]
    }
    eligible <- table$model_id[table$ok & is.finite(table$score)]
    ordered <- eligible[order(vapply(rows[eligible], function(x) x$score, numeric(1)), eligible)]
    final <- stage == length(budgets)
    survivors <- if (final) ordered else utils::head(ordered,
      max(top_k, ceiling(length(active) * survivor_fraction)))
    for (i in setdiff(ordered, survivors)) rows[[i]]$status <- "screened"
    for (i in survivors) rows[[i]]$status <- if (final) "finalist" else "promoted"
    stages[[length(stages) + 1L]] <- list(stage = stage, budget = length(points),
      origins = origins[points], results = table, promoted = survivors)
    active <- survivors
  }
  if (!is.null(cl)) { parallel::stopCluster(cl); cl <- NULL }
  finalists <- active
  fits <- vector("list", n_models)
  refit_failures <- list()
  selected <- integer()
  refit <- function(i) {
    a <- args[[i]]
    a$series <- series
    a$volumes <- volumes
    a$seed <- seeds[i]
    a$n_testing <- 0L
    set.seed(seeds[i])
    fit <- do.call(craft_fit, a)
    fit$backtest <- backtests[[i]]
    fit$diagnostics$backtest <- fit$backtest$summary
    fit$diagnostics$backtest_settings <- fit$backtest$settings
    fit
  }
  refit_ids <- if (keep_all_fits) {
    c(finalists, setdiff(which(vapply(rows, function(x) isTRUE(x$ok), logical(1))), finalists))
  } else finalists
  for (i in refit_ids) {
    if (!keep_all_fits && length(selected) >= top_k) break
    attempt <- tryCatch(list(fit = refit(i)), error = function(e) list(error = conditionMessage(e)))
    if (!is.null(attempt$error)) {
      refit_failures[[length(refit_failures) + 1L]] <- data.frame(model_id = i, error = attempt$error)
      rows[[i]]$status <- "refit_failed"
    } else {
      fits[[i]] <- attempt$fit
      if (i %in% finalists && length(selected) < top_k) selected <- c(selected, i)
    }
  }
  for (i in selected) rows[[i]]$selected <- TRUE
  result_table <- .craft_autotune_bind_rows(rows)
  top_models <- lapply(selected, function(i) list(model_id = i, params = sampled_params[[i]],
    score = rows[[i]]$score, fit = fits[[i]]))
  if (!length(top_models)) warning("No candidate produced a usable final fit; inspect $failures and $refit_failures.")
  out <- list(call = match.call(), score_metric = score_metric, search_space = search_space,
    results = result_table, top_models = top_models,
    failures = result_table[!result_table$ok, , drop = FALSE],
    refit_failures = .v5_bind_rows(refit_failures),
    all_fits = if (keep_all_fits) fits else NULL,
    stages = stages, performance = .v5_bind_rows(performance),
    settings = list(n_models = n_models, top_k = top_k, n_testing = n_testing,
      actual_n_testing = length(origins), origins = origins, seed = seed,
      stage_budgets = budgets, survivor_fraction = survivor_fraction,
      search_parallel = search_parallel, search_workers = workers,
      cache = cache, cache_size = cache_size, keep_all_fits = keep_all_fits,
      composite_weights = composite_weights))
  class(out) <- c("craft_autotune", "list")
  out
}
