fast_data <- function(n = 100) {
  set.seed(21)
  data.frame(a = cumprod(1 + rnorm(n, .001, .01)),
             b = cumprod(1 + rnorm(n, .0005, .012)))
}

fast_args <- function() list(window = 5, horizon = 5, n_draws = 40,
  n_factors = 1, min_segment = 5, max_regimes_per_factor = 3)

fast_search <- function(...) craft_autotune(fast_data(), n_models = 4,
  top_k = 1, n_testing = 3, search_space = list(transition_alpha = c(.5, 1, 2)),
  fixed_args = fast_args(), seed = 42, verbose = FALSE, ...)

test_that("the public tuner is the accelerated engine under one name", {
  expect_true("craft_autotune" %in% getNamespaceExports("CRAFT"))
  expect_false("craft_autotune_fast" %in% getNamespaceExports("CRAFT"))
  expect_false(exists("craft_autotune_fast", envir = asNamespace("CRAFT"), inherits = FALSE))
  expect_true(all(c("stage_budgets", "survivor_fraction", "cache", "cache_size") %in%
                    names(formals(craft_autotune))))
})

test_that("one-stage fast scoring reproduces a full-pipeline scoring reference", {
  old <- reference_autotune(fast_data(), n_models = 4, top_k = 1, n_testing = 3,
    search_space = list(transition_alpha = c(.5, 1, 2)), fixed_args = fast_args(),
    seed = 42, verbose = FALSE)
  fast <- fast_search(stage_budgets = 3)
  expect_s3_class(fast, "craft_autotune")
  expect_equal(fast$results[, names(old$results)], old$results)
  expect_equal(fast$top_models[[1]]$model_id, old$top_models[[1]]$model_id)
  expect_equal(fast$top_models[[1]]$fit$trajectory_draws,
               old$top_models[[1]]$fit$trajectory_draws)
  expect_equal(fast$top_models[[1]]$fit$backtest$summary,
               old$top_models[[1]]$fit$backtest$summary)
  expect_s3_class(craft_predict(fast$top_models[[1]]$fit, n_draws = 10), "craft_prediction")
})

test_that("score-only fits bypass distributions and diagnostics", {
  args <- fast_args()
  args$series <- fast_data()
  args$verbose <- FALSE
  full <- do.call(CRAFT:::.craft_fit_impl, args)
  testthat::local_mocked_bindings(
    trajectory_forecast = function(...) stop("distribution fitting was invoked"),
    .package = "CRAFT")
  args$.score_only <- TRUE
  scored <- do.call(CRAFT:::.craft_fit_impl, args)
  expect_equal(scored$trajectory_draws, full$trajectory_draws)
  expect_equal(scored$valid_acc, full$valid_acc)
  expect_null(scored$return_dists)
  expect_null(scored$corr_by_horizon)
})

test_that("stages reuse origins and distinguish screening from final scores", {
  staged <- fast_search()
  full <- fast_search(stage_budgets = 3)
  expect_length(staged$stages, 2)
  expect_equal(staged$stages[[1]]$budget, 2)
  expect_equal(staged$stages[[2]]$budget, 3)
  expect_equal(sum(staged$performance$origins_computed), sum(staged$results$n_testing))
  expect_lt(sum(staged$performance$origins_computed), 4 * 3)
  finalists <- staged$results$status == "finalist"
  expect_true(any(staged$results$status == "screened"))
  expect_equal(staged$results$score[finalists], full$results$score[finalists])
  expect_equal(staged$results$n_testing[staged$results$selected], 3)
  expect_equal(sum(staged$results$selected), 1)
  expect_equal(staged$results, fast_search()$results)
})

test_that("cache hits, eviction, and disabling caching preserve scores", {
  cached <- fast_search(stage_budgets = 3)
  plain <- fast_search(stage_budgets = 3, cache = FALSE)
  tiny <- fast_search(stage_budgets = 3, cache_size = 1)
  expect_gt(sum(cached$performance$cache_hits), 0)
  expect_equal(cached$results, plain$results)
  expect_equal(cached$results, tiny$results)
  expect_equal(sum(plain$performance$cache_hits), 0)
  store <- CRAFT:::.craft_fast_cache(1)
  counter <- 0L
  compute <- function() { counter <<- counter + 1L; counter }
  get_value <- CRAFT:::.craft_fast_cached
  expect_equal(get_value(store, list(1), compute), 1)
  expect_equal(get_value(store, list(1), compute), 1)
  expect_equal(get_value(store, list(2), compute), 2)
  expect_equal(get_value(store, list(1), compute), 3)
  expect_equal(get_value(store, list(raw(20000)), compute), 4)
  expect_equal(get_value(store, list(raw(20000)), compute), 4)
  expect_length(ls(store$values), 1)
})

test_that("stochastic regime calls bypass the deterministic cache", {
  calls <- 0L
  testthat::local_mocked_bindings(svd_changepoint_regimes = function(...) {
    calls <<- calls + 1L; calls
  }, .package = "CRAFT")
  cache <- CRAFT:::.craft_fast_cache(4)
  run <- CRAFT:::.craft_fast_regimes
  expect_equal(run(cache, list("a"), list(X = matrix(1, 2, 2))), 1)
  expect_equal(run(cache, list("a"), list(X = matrix(1, 2, 2))), 1)
  expect_equal(run(cache, list("a"), list(X = matrix(1, 2, 2), stability_bootstrap = TRUE)), 2)
  expect_equal(run(cache, list("a"), list(X = matrix(1, 2, 2), seed = 1)), 3)
})

test_that("only winners are fully fitted unless all fits are requested", {
  impl <- CRAFT:::.craft_fit_impl
  full_calls <- 0L
  testthat::local_mocked_bindings(.craft_fit_impl = function(...) {
    args <- list(...)
    if (!isTRUE(args$.score_only)) full_calls <<- full_calls + 1L
    do.call(impl, args)
  }, .package = "CRAFT")
  result <- fast_search(stage_budgets = 3)
  expect_equal(full_calls, 1)
  expect_null(result$all_fits)
  full_calls <- 0L
  result <- fast_search(keep_all_fits = TRUE)
  expect_equal(full_calls, sum(result$results$ok))
  expect_length(result$all_fits, 4)
  expect_true(all(vapply(result$all_fits, inherits, logical(1), "craft_fit")))
})

test_that("failed candidates and failed final refits are inspectable", {
  count <- 0L
  spec <- function() { count <<- count + 1L; c(5, 500)[count] }
  mixed <- craft_autotune(fast_data(), n_models = 2, top_k = 1,
    n_testing = 1, search_space = list(window = spec), fixed_args = fast_args(), verbose = FALSE)
  expect_equal(mixed$results$ok, c(TRUE, FALSE))
  expect_equal(mixed$failures$model_id, 2)
  expect_warning(all_bad <- craft_autotune(fast_data(), n_models = 1,
    n_testing = 1, search_space = list(window = 500), verbose = FALSE), "No candidate")
  expect_length(all_bad$top_models, 0)
  impl <- CRAFT:::.craft_fit_impl
  testthat::local_mocked_bindings(.craft_fit_impl = function(...) {
    args <- list(...)
    if (!isTRUE(args$.score_only)) stop("refit failure")
    do.call(impl, args)
  }, .package = "CRAFT")
  expect_warning(failed <- fast_search(stage_budgets = 3), "No candidate")
  expect_equal(nrow(failed$refit_failures), 4)
  expect_true(all(failed$results$status == "refit_failed"))
})

test_that("common origins and structural cache keys work across windows", {
  run <- function(cache) {
    counter <- 0L
    spec <- function() { counter <<- counter + 1L; c(5, 6, 5, 6)[counter] }
    fixed <- fast_args()
    fixed$backtest_min_train <- 30
    craft_autotune(fast_data(), n_models = 4, top_k = 1, n_testing = 2,
      search_space = list(window = spec), fixed_args = fixed, verbose = FALSE, cache = cache)
  }
  cached <- run(TRUE)
  expect_equal(cached$results, run(FALSE)$results)
  expect_equal(cached$settings$origins, c(30L, 95L))
})

test_that("fast controls and composite scores are validated", {
  expect_error(fast_search(stage_budgets = c(3, 2)), "stage_budgets")
  expect_error(fast_search(survivor_fraction = 0), "survivor_fraction")
  expect_error(fast_search(cache_size = 0), "cache_size")
  expect_error(fast_search(cache = NA), "cache")
  expect_error(fast_search(composite_weights = c(foo = 1)), "composite_weights")
  composite <- fast_search(score_metric = "composite")
  expect_true(all(is.finite(composite$results$score)))
  expect_length(composite$top_models, 1)
})

test_that("dynamic parallel search reproduces serial staged search", {
  skip_on_cran()
  serial <- fast_search()
  parallel <- fast_search(search_parallel = TRUE, search_workers = 2)
  expect_equal(parallel$results, serial$results)
  expect_equal(parallel$stages, serial$stages)
  expect_equal(parallel$top_models[[1]]$fit$trajectory_draws,
               serial$top_models[[1]]$fit$trajectory_draws)
})

test_that("the original parallel backtest still works after core extraction", {
  skip_on_cran()
  args <- fast_args()
  args$series <- fast_data()
  args$n_testing <- 2
  args$verbose <- FALSE
  serial <- do.call(craft_fit, c(args, list(backtest_parallel = FALSE)))
  parallel <- do.call(craft_fit, c(args, list(backtest_parallel = TRUE, backtest_workers = 2)))
  expect_equal(parallel$backtest$summary, serial$backtest$summary)
})

test_that("full horizon scoring is preserved without the cache", {
  fixed <- fast_args()
  fixed$horizon <- NULL
  args <- list(series = fast_data(), n_models = 2, top_k = 1, n_testing = 2,
               search_space = list(transition_alpha = c(.5, 1)),
               fixed_args = fixed, verbose = FALSE)
  old <- do.call(reference_autotune, args)
  fast <- do.call(craft_autotune, c(args, list(cache = FALSE)))
  expect_equal(fast$results$score, old$results$score)
})

test_that("craft_fit still evaluates arguments in the caller's frame", {
  fit <- local({
    local_history <- fast_data()
    local_window <- 5
    craft_fit(head(local_history, 80), window = local_window,
      n_draws = 20, n_factors = 1, min_segment = 5,
      max_regimes_per_factor = 3, verbose = FALSE)
  })
  expect_s3_class(fit, "craft_fit")
  expect_equal(fit$window, 5)
})
