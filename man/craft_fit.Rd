\name{craft_fit}
\alias{craft_fit}
\title{Fit a Trajectory-Aware Regime Forecast}
\usage{
craft_fit(
  series,
  volumes = NULL,
  window = 20,
  n_draws = 1000,
  n_factors = 3,
  var_threshold = 0.95,
  horizon = NULL,
  seed = 42,
  transition_alpha = 1,
  transition_min_count = 1,
  transition_key_sep = "\\r",
  valid_frac = 0.2,
  valid_split = c("last", "random", "blocked"),
  n_testing = 0L,
  backtest_parallel = TRUE,
  backtest_workers = NULL,
  backtest_min_train = NULL,
  sampler_weight = c("uniform", "weighted_bal_acc", "bal_acc", "weighted_acc", "acc"),
  tail_penalty_lambda = 0,
  tail_penalty_threshold = 3,
  tail_penalty_power = 2,
  tail_penalty_aggregation = c("mean", "max", "sum"),
  tail_penalty_scope = c("forecast", "all"),
  volume_alpha = 1,
  max_regimes_per_factor = 6,
  factor_selection = "fixed",
  min_segment = 30L,
  past_regime_args = list(),
  future_regime_args = list(),
  return_train_probs = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{series}{A data frame or matrix of numeric time series in columns.}
\item{volumes}{Optional matrix or data frame of volumes with the same dimensions as \code{series}.}
\item{window}{Trajectory window length used for lag and lead embeddings.}
\item{n_draws}{Number of future-trajectory rows to sample.}
\item{n_factors}{Number of SVD factors when fixed factor selection is used.}
\item{var_threshold}{Cumulative variance threshold for variance-based factor selection.}
\item{horizon}{Optional forecast horizon. If \code{NULL}, all horizons from 1 to \code{window} are retained.}
\item{seed}{Optional random seed.}
\item{transition_alpha}{Additive smoothing parameter for transition probabilities.}
\item{transition_min_count}{Minimum count for memorized transition rows before fallback to global probabilities.}
\item{transition_key_sep}{Internal separator used for transition-state keys.}
\item{valid_frac}{Fraction of rows reserved for transition-model validation.}
\item{valid_split}{Validation split strategy.}
\item{n_testing}{Number of internal backtest origins. Use zero to skip backtesting.}
\item{backtest_parallel}{Logical. If \code{TRUE}, run internal backtests in parallel when possible.}
\item{backtest_workers}{Optional number of parallel workers.}
\item{backtest_min_train}{Optional minimum training length for backtest origins.}
\item{sampler_weight}{Metric used to weight target regimes during future-trajectory sampling.}
\item{tail_penalty_lambda}{Strength of the optional tail penalty applied to candidate rows.}
\item{tail_penalty_threshold}{Robust z-score threshold for tail penalties.}
\item{tail_penalty_power}{Power applied to excess robust z-scores.}
\item{tail_penalty_aggregation}{Aggregation used for row-level tail penalties.}
\item{tail_penalty_scope}{Whether the tail penalty uses forecast-horizon columns or all forward-trajectory columns.}
\item{volume_alpha}{Exponent for volume-derived time weights.}
\item{max_regimes_per_factor}{Maximum regimes per SVD factor.}
\item{factor_selection}{SVD factor selection strategy passed to \code{svd_changepoint_regimes()}.}
\item{min_segment}{Minimum segment length for changepoint regime detection.}
\item{past_regime_args}{Additional arguments passed to past-trajectory regime detection.}
\item{future_regime_args}{Additional arguments passed to future-trajectory regime detection.}
\item{return_train_probs}{Logical. If \code{TRUE}, include training transition probabilities.}
\item{verbose}{Logical. If \code{TRUE}, emit progress messages.}
}
\description{
Runs the CRAFT end-to-end workflow: trajectory embedding, past and future
regime detection, transition probability estimation, future-trajectory
sampling, correlation diagnostics, distribution fitting, and optional
internal backtesting.
}
\value{
An object of class \code{craft_fit}, represented as a
list containing trajectories, regime models, labels, transition probabilities,
sampled trajectory draws, fitted return and level distributions, diagnostics,
classifier details, optional backtest results, validation metrics, and the information needed by \code{craft_predict()}, including \code{trajectory_draws}.
}
\examples{
\donttest{
set.seed(1)
series <- data.frame(
  a = cumprod(1 + rnorm(80, 0.001, 0.01)),
  b = cumprod(1 + rnorm(80, 0.0005, 0.012))
)
fit <- craft_fit(
  series,
  window = 5,
  n_draws = 50,
  n_factors = 1,
  min_segment = 5,
  max_regimes_per_factor = 3,
  n_testing = 0,
  verbose = FALSE
)
names(fit$return_dists)
}
}
