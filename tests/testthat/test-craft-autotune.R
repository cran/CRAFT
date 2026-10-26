tuning_series <- function() {
  set.seed(21)
  data.frame(a = cumprod(1 + rnorm(80, 0.001, 0.01)),
             b = cumprod(1 + rnorm(80, 0.0005, 0.012)))
}

tuning_args <- function() {
  list(window = 5, horizon = 5, n_draws = 30, n_factors = 1,
       min_segment = 5, max_regimes_per_factor = 3)
}

test_that("search specifications preserve singleton numeric choices", {
  draw <- CRAFT:::.craft_autotune_draw_value
  expect_equal(draw(list(type = "choice", values = 5)), 5)
  expect_equal(draw(list(type = "integer", min = 5, max = 5)), 5)
  expect_equal(draw(5), 5)
  expect_equal(draw(function() 7), 7)
  expect_equal(draw(list(type = "numeric", min = 2, max = 2)), 2)
  expect_true(draw(list(type = "numeric", min = 0.1, max = 1, log = TRUE)) >= 0.1)
  expect_equal(draw(list(a = 1)), list(a = 1))
})

test_that("invalid tuning controls fail before fitting", {
  x <- tuning_series()
  for (v in list(0, NA_real_, Inf, 1.5, c(1, 2), numeric())) {
    expect_error(craft_autotune(x, n_models = v), "n_models")
  }
  expect_error(craft_autotune(x, n_testing = 0), "n_testing")
  expect_error(craft_autotune(x, seed = NULL), "seed")
  expect_error(craft_autotune(x, keep_all_fits = NA), "keep_all_fits")
  expect_error(craft_autotune(x, search_workers = 1.5), "search_workers")
  expect_error(craft_autotune(x, search_space = list(1)), "unique")
  expect_error(craft_autotune(x, search_space = list(foo = 1)), "Unknown")
  expect_error(craft_autotune(x, fixed_args = list(seed = 1)), "reserved")
  expect_error(craft_autotune(x, foo = 1), "Unknown")
  expect_error(craft_autotune(x, search_space = list(window = 1, window = 2)), "unique")
  expect_error(craft_autotune(x, search_space = list(window = list(type = "integer", min = 4, max = 2))), "min")
  expect_error(craft_autotune(x, search_space = list(window = list(type = "choice", values = numeric()))), "values")
  expect_error(craft_autotune(x, search_space = list(window = list(type = "numeric", min = 0, max = 2, log = TRUE))), "positive")
  expect_error(craft_autotune(x, composite_weights = c(foo = 1)), "composite_weights")
  expect_error(craft_autotune(x, composite_weights = c(return_rmse = -1)), "composite_weights")
  expect_error(craft_autotune(x, search_parallel = TRUE,
                            search_space = list(backtest_parallel = TRUE)), "Cannot tune")
})

test_that("tuning ranks current craft fits reproducibly and predicts", {
  args <- list(series = tuning_series(), n_models = 3, top_k = 2,
               n_testing = 2, search_space = list(transition_alpha = c(0.5, 1, 2)),
               fixed_args = tuning_args(), seed = 77, verbose = FALSE,
               keep_all_fits = TRUE)
  a <- do.call(craft_autotune, args)
  b <- do.call(craft_autotune, args)
  expect_s3_class(a, "craft_autotune")
  expect_equal(a$results, b$results)
  expect_true(all(a$results$ok))
  expect_equal(a$results$model_id, 1:3)
  expect_equal(nrow(a$failures), 0)
  expect_length(a$top_models, 2)
  expect_length(a$all_fits, 3)
  expect_equal(vapply(a$top_models, `[[`, numeric(1), "score"),
               sort(a$results$score)[1:2])
  fit <- a$top_models[[1]]$fit
  expect_s3_class(fit, "craft_fit")
  expect_s3_class(craft_predict(fit, n_draws = 10), "craft_prediction")
  expect_equal(fit$backtest$settings$actual_n_testing, 2)
})

test_that("fixed arguments and sampled arguments have documented precedence", {
  a <- craft_autotune(tuning_series(), n_models = 1, n_testing = 1,
                      search_space = list(window = 5),
                      fixed_args = tuning_args(), window = 8, verbose = FALSE)
  expect_equal(a$top_models[[1]]$fit$window, 5)
  b <- craft_autotune(tuning_series(), n_models = 1, n_testing = 1,
                      search_space = list(), fixed_args = tuning_args(), verbose = FALSE)
  expect_true(b$results$ok)
  expect_null(b$all_fits)
})

test_that("failed candidates are recorded and excluded", {
  count <- 0L
  window_draw <- function() { count <<- count + 1L; c(5, 500)[count] }
  a <- craft_autotune(tuning_series(), n_models = 2, n_testing = 1,
                      search_space = list(window = window_draw),
                      fixed_args = tuning_args(), keep_all_fits = TRUE, verbose = FALSE)
  expect_equal(a$results$ok, c(TRUE, FALSE))
  expect_equal(a$results$score[2], Inf)
  expect_length(a$top_models, 1)
  expect_equal(a$failures$model_id, 2)
  expect_null(a$all_fits[[2]])
  expect_warning(b <- craft_autotune(tuning_series(), n_models = 1, n_testing = 1,
                     search_space = list(window = 500), verbose = FALSE), "No candidate")
  expect_length(b$top_models, 0)
  expect_false(b$results$ok)
})

test_that("metric direction and composite ranks follow their definitions", {
  score <- CRAFT:::.craft_autotune_score_fit
  overall <- data.frame(return_rmse = 2, return_mae = 1, level_rmse = 4,
                        level_mae = 3, coverage_90 = 0.8, direction_accuracy = 0.75)
  fit <- list(backtest = list(summary = list(overall = overall)))
  expect_equal(score(fit, "direction_accuracy")$score, -0.75)
  expect_equal(score(fit, "coverage_90_gap")$score, 0.1)
  expect_equal(score(fit, "level_mae")$score, 3)
  expect_error(score(list()), "empty")
  rank_scores <- CRAFT:::.craft_autotune_composite_scores
  d <- data.frame(ok = c(TRUE, TRUE, FALSE), return_rmse = c(1, 3, NA),
                  return_mae = c(2, 1, NA))
  expect_equal(rank_scores(d), c(0.5, 0.5, Inf))
  expect_equal(rank_scores(d, c(return_rmse = 1, return_mae = 0)), c(0, 1, Inf))
})

test_that("parallel installed namespace search matches serial search", {
  skip_on_cran()
  args <- list(series = tuning_series(), n_models = 2, n_testing = 1,
               search_space = list(transition_alpha = c(0.5, 1)),
               fixed_args = tuning_args(), seed = 5, verbose = FALSE)
  a <- do.call(craft_autotune, args)
  b <- do.call(craft_autotune, c(args, list(search_parallel = TRUE, search_workers = 2)))
  expect_equal(a$results, b$results)
  expect_true(all(b$results$ok))
})
