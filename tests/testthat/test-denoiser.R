test_that("denoiser returns selected data with report metadata", {
  set.seed(1)
  x <- data.frame(
    trend = seq_len(40) + rnorm(40, sd = 0.1),
    noise = rnorm(40),
    constant = 1,
    label = rep(c("a", "b"), each = 20)
  )

  out <- denoiser(x, min_features = 1, verbose = FALSE, return_report = TRUE)

  expect_named(out, c("data", "features", "report", "settings"))
  expect_gte(length(out$features), 1)
  expect_true(all(out$features %in% names(out$data)))
  expect_true(all(c("feature", "selected", "score", "reason") %in% names(out$report)))
})
