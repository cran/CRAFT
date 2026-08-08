## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(CRAFT)

set.seed(42)
n <- 120
market_pulse <- cumsum(rnorm(n, mean = 0, sd = 0.004))
rate_pulse <- cumsum(rnorm(n, mean = 0, sd = 0.003))

series <- data.frame(
  asset_a = cumprod(1 + 0.0010 + market_pulse + rnorm(n, 0, 0.008)),
  asset_b = cumprod(1 + 0.0005 + 0.6 * market_pulse - 0.2 * rate_pulse +
                      rnorm(n, 0, 0.010)),
  asset_c = cumprod(1 + 0.0008 - 0.3 * market_pulse + 0.5 * rate_pulse +
                      rnorm(n, 0, 0.007))
)

tail(series, 3)

## ----trajectories-------------------------------------------------------------
trajectories <- trajectory_embedding(series, trajectory_window = 5)

names(trajectories)
dim(trajectories$past_trajectories)
dim(trajectories$future_trajectories)

head(round(trajectories$past_trajectories[, 1:5], 4), 3)

## ----trajectory-forecast------------------------------------------------------
asset_a_cols <- grep(
  "^asset_a_cum_lead_",
  colnames(trajectories$future_trajectories),
  value = TRUE
)
asset_a_draws <- as.data.frame(trajectories$future_trajectories[, asset_a_cols, drop = FALSE])

asset_a_fc <- trajectory_forecast(
  asset_a_draws,
  probs = seq(0.05, 0.95, length.out = 40),
  min_unique = 5,
  verbose = FALSE
)

names(asset_a_fc)
asset_a_fc[["asset_a_cum_lead_5"]]$qfun(c(0.05, 0.50, 0.95))

## ----craft-fit----------------------------------------------------------------
fit <- craft_fit(
  series,
  window = 5,
  n_draws = 80,
  n_factors = 1,
  min_segment = 5,
  max_regimes_per_factor = 3,
  n_testing = 0,
  return_train_probs = FALSE,
  verbose = FALSE,
  seed = 123
)

class(fit)
fit$valid_joint_acc
names(fit$return_dists)

## ----inspect-fit--------------------------------------------------------------
dim(fit$trajectory_draws)
names(fit$diagnostics)

round(fit$sampler_weights, 3)
round(fit$diagnostics$posterior_entropy, 3)

## ----craft-predict------------------------------------------------------------
pred <- craft_predict(fit, n_draws = 50, seed = 7)

class(pred)
dim(pred$trajectory_draws)

horizon_name <- tail(names(pred$return_dists$asset_a), 1)
pred$return_dists$asset_a[[horizon_name]]$qfun(c(0.05, 0.50, 0.95))

## ----prediction-types---------------------------------------------------------
labels <- craft_predict(fit, type = "labels")
probs <- craft_predict(fit, type = "probabilities")
draws <- craft_predict(fit, type = "trajectory_draws", n_draws = 10, seed = 9)

labels
lapply(probs, function(x) round(x[1, ], 3))
dim(draws)

