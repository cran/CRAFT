test_that("svd_changepoint_regimes returns regimes and transformation helpers", {
  skip_if_not_installed("changepoint")
  set.seed(1)
  x <- data.frame(
    a = c(rnorm(35, 0, 0.2), rnorm(35, 2, 0.2)),
    b = c(rnorm(35, 0, 0.2), rnorm(35, -1, 0.2))
  )

  fit <- svd_changepoint_regimes(
    x,
    denoise = FALSE,
    reduction_criteria = "fixed",
    n_dim = 1,
    minseglen = 5,
    min_regime_size = 5,
    max_regimes_per_dim = 3,
    seed = 1
  )

  expect_s3_class(fit, "svd_changepoint_regimes_v2")
  expect_true(is.data.frame(fit$regimes))
  expect_true(is.data.frame(fit$diagnostics))
  expect_true(is.function(fit$transform_newdata))
  transformed <- fit$transform_newdata(tail(x, 2))
  expect_true(is.data.frame(transformed) || is.list(transformed))
})

test_that("trajectory_forecast fits one distribution per draw column", {
  set.seed(1)
  trajectory_draws <- data.frame(a = rnorm(50, 0.01, 0.02), b = rnorm(50, -0.01, 0.03))

  fit <- trajectory_forecast(
    trajectory_draws,
    probs = seq(0.05, 0.95, length.out = 40),
    min_unique = 5,
    verbose = FALSE
  )

  expect_s3_class(fit, "smooth_trajectory_forecast")
  expect_named(fit, names(trajectory_draws))
  expect_true(is.list(attr(fit, "diagnostics")))
})
