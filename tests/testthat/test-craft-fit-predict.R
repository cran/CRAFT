test_that("craft_fit returns a forecast object and craft_predict reuses it", {
  skip_if_not_installed("changepoint")
  set.seed(2)
  n <- 80
  series <- data.frame(
    a = cumprod(1 + rnorm(n, 0.001, 0.01)),
    b = cumprod(1 + rnorm(n, 0.0005, 0.012))
  )

  fit <- craft_fit(
    series,
    window = 5,
    n_draws = 40,
    n_factors = 1,
    min_segment = 5,
    max_regimes_per_factor = 3,
    n_testing = 0,
    verbose = FALSE,
    return_train_probs = FALSE
  )

  expect_s3_class(fit, "craft_fit")
  expect_true(is.data.frame(fit$trajectory_draws))
  expect_named(fit$return_dists, names(series))
  expect_named(fit$level_dists, names(series))
  expect_named(fit$last_level, names(series))
  expect_true(is.list(fit$diagnostics))

  pred <- craft_predict(fit, n_draws = 20, seed = 3)
  expect_s3_class(pred, "craft_prediction")
  expect_equal(nrow(pred$trajectory_draws), 20)
  expect_named(pred$return_dists, names(series))
  expect_named(pred$level_dists, names(series))

  pred_new <- craft_predict(fit, newdata = tail(series, 10), n_draws = 20, seed = 4)
  expect_s3_class(pred_new, "craft_prediction")
  expect_equal(nrow(pred_new$trajectory_draws), 20)

  probs <- craft_predict(fit, type = "probabilities")
  expect_named(probs, fit$targets)

  trajectory_draws <- craft_predict(fit, n_draws = 10, seed = 5, type = "trajectory_draws")
  expect_equal(nrow(trajectory_draws), 10)
})
