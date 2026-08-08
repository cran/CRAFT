test_that("time_volume_weights returns normalized non-negative weights", {
  vol <- data.frame(day = as.Date("2024-01-01") + 0:2, a = c(10, 20, 30), b = c(5, 5, 15))
  w <- time_volume_weights(vol, date_col = "day")

  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(w >= 0))
  expect_equal(names(w), as.character(vol$day))
})

test_that("smooth_quantile_distribution exposes density, distribution, quantile and random functions", {
  set.seed(1)
  fit <- smooth_quantile_distribution(rnorm(120), probs = seq(0.05, 0.95, length.out = 50))

  expect_named(fit, c("data", "dfun", "pfun", "qfun", "rfun", "q_fit", "d_fit", "original_x"))
  expect_true(is.data.frame(fit$data))
  expect_true(is.function(fit$dfun))
  expect_true(is.function(fit$pfun))
  expect_true(is.function(fit$qfun))
  expect_true(is.function(fit$rfun))
  expect_true(all(is.finite(fit$qfun(c(0.25, 0.5, 0.75)))))
  expect_length(fit$rfun(5), 5)
})
