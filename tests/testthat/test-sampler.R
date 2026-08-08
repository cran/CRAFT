test_that("regime_sampler samples rows and indices from named probabilities", {
  trafo <- data.frame(value = 1:4)
  clusters <- data.frame(regime = c("low", "low", "high", "high"))
  probs <- list(regime = c(low = 0.1, high = 0.9))

  rows <- regime_sampler(trafo, clusters, probs, size = 20, seed = 1)
  idx <- regime_sampler(trafo, clusters, probs, size = 20, return_type = "index", seed = 1)

  expect_s3_class(rows, "data.frame")
  expect_equal(nrow(rows), 20)
  expect_type(idx, "integer")
  expect_true(all(idx %in% seq_len(nrow(trafo))))
})
