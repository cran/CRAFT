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

`craft_autotune()` uses staged random search with scoring-only backtests,
bounded caching, dynamically scheduled parallel workers, and deferred full
fitting of selected candidates. The accelerated engine is now the only public
tuner; calls previously using `craft_autotune_fast()` should use this name.

```r
tuned <- craft_autotune(
  series, n_models = 12, top_k = 2, n_testing = 5,
  stage_budgets = c(2, 5), survivor_fraction = 0.3,
  search_space = list(transition_alpha = c(0.25, 0.5, 1, 2),
                      sampler_weight = c("uniform", "bal_acc")),
  fixed_args = list(window = 5, horizon = 5, n_draws = 250,
                   n_factors = 1, min_segment = 5,
                   max_regimes_per_factor = 3),
  search_parallel = TRUE, search_workers = 2,
  seed = 42, verbose = FALSE
)
tuned$results[, c("model_id", "score", "n_testing", "status", "selected")]
if (length(tuned$top_models)) {
  prediction <- craft_predict(tuned$top_models[[1]]$fit)
}
```

An explicit `search_space` replaces the default search space. In the example,
only transition smoothing and sampler weights are searched; `window` stays at
5. Sampled settings override `fixed_args`, which override `...`, so omit fixed
parameters from the search space. Specifications support scalar constants,
candidate vectors, zero-argument sampler functions, and typed choice, integer,
numeric, or log-scale ranges. An empty list repeats the fixed configuration
with different model seeds.

Screened candidates have lower-budget scores; compare finalists when selecting
a model. `top_models` does this automatically. Set `stage_budgets = n_testing`
to evaluate every candidate at the full budget. Screening and fewer draws can
change rankings. Reserve a later holdout period for performance evaluation.

The tuner uses common origins across feasible candidates. Fix `horizon` when
tuning windows so candidates also forecast the same outcomes. Caches are scoped
to each search, worker, and training prefix. They reuse trajectory embeddings
and identical deterministic regime fits, including SVD work. `cache_size`
limits entries per worker, not bytes. Bootstrap, explicitly seeded, or plotting
regime calls bypass regime caching.

Parallel workers receive the history once and return compact scoring data.
Nested parallelism and training-probability output are disabled. By default,
only selected finalists are fitted in full; `keep_all_fits = TRUE` additionally
refits screened candidates, increasing runtime and memory use.

Inspect `failures`, `refit_failures`, `stages`, and `performance` for diagnostics.
The installed benchmark at
`system.file("benchmarks", "compare_autotune_speed.R", package = "CRAFT")`
compares cached and uncached full-budget tuning with staged tuning.

## License

GPL-3
