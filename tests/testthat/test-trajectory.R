test_that("trajectory_reframer builds expected past and future trajectories", {
  x <- c(100, 110, 121, 100)

  past <- trajectory_reframer(x, k = 2, mode = "backward", ts_name = "x")
  expect_s3_class(past, "trajectory_reframer")
  expect_equal(colnames(past$trajectory), c("x_cum_lag_1", "x_cum_lag_2"))
  expect_equal(unname(past$trajectory[2, 1]), 0.1, tolerance = 1e-12)
  expect_equal(unname(past$trajectory[3, 2]), 0.21, tolerance = 1e-12)

  future <- trajectory_reframer(x, k = 1, mode = "forward", ts_name = "x")
  expect_equal(unname(future$trajectory[1, 1]), 0.1, tolerance = 1e-12)
})

test_that("trajectory_embedding aligns past and future trajectory rows", {
  x <- data.frame(a = c(100, 102, 101, 104, 105), b = c(50, 51, 53, 54, 55))
  trajectories <- trajectory_embedding(x, trajectory_window = 2)

  expect_named(
    trajectories,
    c("past_trajectories", "future_trajectories", "latest_trajectory", "time_index")
  )
  expect_equal(ncol(trajectories$past_trajectories), 4)
  expect_equal(ncol(trajectories$future_trajectories), 4)
  expect_equal(nrow(trajectories$past_trajectories), nrow(trajectories$future_trajectories))
  expect_equal(nrow(trajectories$latest_trajectory), 1)
  expect_type(trajectories$time_index, "logical")
})
