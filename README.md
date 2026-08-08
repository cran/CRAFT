# CRAFT

CRAFT provides Conditional Regime Analog Forecasting with Trajectories tools for multivariate time
series. It builds lag and lead trajectory embeddings, detects SVD-based
changepoint regimes, estimates transition probabilities, samples plausible
future-trajectory analogues, and fits smooth empirical forecast distributions.

## Installation

After CRAN release:

```r
install.packages("CRAFT")
```

For local development from this folder:

```r
install.packages(".", repos = NULL, type = "source")
```

## Basic workflow

```r
library(CRAFT)

set.seed(1)
series <- data.frame(
  asset_a = cumprod(1 + rnorm(100, 0.001, 0.01)),
  asset_b = cumprod(1 + rnorm(100, 0.0005, 0.012))
)

fit <- craft_fit(
  series,
  window = 5,
  n_draws = 100,
  n_factors = 1,
  min_segment = 5,
  max_regimes_per_factor = 3,
  n_testing = 0,
  verbose = FALSE
)

fit$valid_joint_acc
pred <- craft_predict(fit, n_draws = 25, seed = 1)

names(pred$return_dists)
```

Lower-level functions are available for custom pipelines:

- `trajectory_embedding()` builds aligned past and future trajectory matrices.
- `svd_changepoint_regimes()` detects regimes in transformed time-series data.
- `regime_sampler()` samples historical future-trajectory rows from regime probabilities.
- `trajectory_forecast()` fits smooth empirical forecast distributions.
- `denoiser()` selects informative numeric features before regime detection.
- `time_volume_weights()` creates volume-derived observation weights.

## CRAN submission note

Before submitting to CRAN, replace the placeholder maintainer email in
`DESCRIPTION` with the real maintainer email address.

## License

MIT + file LICENSE
