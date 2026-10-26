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

## Hyperparameter tuning

`craft_autotune()` samples candidate settings and ranks their internal backtest
scores. Lower scores are better; direction accuracy is negated for ranking.

```r
tuned <- craft_autotune(
  series, n_models = 2, top_k = 1, n_testing = 2,
  search_space = list(transition_alpha = c(0.5, 1, 2)),
  fixed_args = list(window = 5, horizon = 5, n_draws = 40,
                   n_factors = 1, min_segment = 5,
                   max_regimes_per_factor = 3),
  seed = 123, verbose = FALSE
)
tuned$results
if (length(tuned$top_models)) {
  best_fit <- tuned$top_models[[1]]$fit
  prediction <- craft_predict(best_fit, n_draws = 25, seed = 7)
}
```

Search specifications may be candidate vectors, scalar constants, zero-argument
sampler functions, or lists such as `list(type = "numeric", min = 0.1,
max = 10, log = TRUE)`. Sampled settings override `fixed_args`, which override
arguments in `...`. An empty `search_space = list()` repeats the fixed settings
with different model seeds. Optional `search_parallel = TRUE` uses two workers
by default and disables nested backtest parallelism. Serial backtests are also
the default during serial search; pass `backtest_parallel = TRUE` explicitly
to enable them.

Keep `horizon` and `backtest_min_train` fixed when comparing different windows
or segment lengths so candidates are evaluated on the same outcomes. Inspect
`tuned$failures` before using the selected fit. These backtests select parameters;
use a separate, later holdout period to assess the selected model's performance.

## License

GPL-3
