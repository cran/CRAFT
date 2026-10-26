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
