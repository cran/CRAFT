# CRAFT 0.4.0

- Replaced the original `craft_autotune()` with the accelerated implementation.
  Staged screening, scoring-only backtests, bounded caches, dynamic parallel
  scheduling, and deferred full fitting are now available under this one name.
- Removed the separate `craft_autotune_fast()` function and export. Existing
  calls to that name should be changed to `craft_autotune()`.
- The returned class is `craft_autotune`. To evaluate every candidate at the
  final budget, set `stage_budgets = n_testing`.
- Consolidated source, reference documentation, README, and vignette; updated
  regression tests and the benchmark for the single public tuning interface.

# CRAFT 0.3.0

- Added `craft_autotune_fast()` with scoring-only backtests that skip forecast
  distributions and correlation diagnostics, then fit selected models in full.
- Added staged candidate screening with stable seeds and reuse of completed
  origins, plus a common backtest-origin grid across candidate configurations.
- Added bounded caches for trajectory embeddings and deterministic regime fits,
  including their SVD computations. Stochastic and plotting regime calls bypass
  regime caching. Cache entries never cross search calls or training prefixes.
- Added dynamically scheduled parallel scoring with one history transfer per
  worker, compact scoring results, and no nested parallelism.
- Added regression tests, usage documentation, and a reproducible tuner benchmark.
- Preserved `craft_fit()` and `craft_autotune()` public interfaces.

# CRAFT 0.2.0

- Added `craft_autotune()` random-search hyperparameter tuning for `craft_fit()`.
- Added scalar, choice, integer, numeric, log-scale, and function search spaces.
- Added backtest metric and composite ranking, retained top fits, candidate failure
  diagnostics, and optional parallel search with bounded worker counts.
- Added tuning examples, reference documentation, and regression tests.
- Corrected stale license and maintainer information in the README.

# CRAFT 0.1.0

- Initial CRAN-preparation release.
- Added trajectory embedding helpers for lag and lead cumulative changes.
- Added SVD changepoint regime detection and transition-based regime forecasting through `craft_fit()` and `craft_predict()`.
- Added smooth empirical forecast distribution fitting.
- Added tests, manual pages, README, and vignette.
