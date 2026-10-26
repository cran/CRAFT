\name{craft_autotune}
\alias{craft_autotune}
\title{Tune CRAFT Hyperparameters by Random Search}
\description{
Samples hyperparameter configurations, fits each with \code{craft_fit()}, and
ranks candidates using internal rolling-origin backtest summaries.
}
\usage{
craft_autotune(
  series, volumes = NULL, n_models = 30L, top_k = 5L, n_testing = 5L,
  search_space = NULL,
  score_metric = c("return_rmse", "return_mae", "level_rmse", "level_mae",
    "coverage_90_gap", "direction_accuracy", "composite"),
  seed = 42, fixed_args = list(), composite_weights = NULL,
  keep_all_fits = FALSE, search_parallel = FALSE, search_workers = NULL,
  verbose = TRUE, ...
)
}
\arguments{
\item{series}{A numeric matrix or data frame of time series in columns.}
\item{volumes}{Optional volumes passed to \code{craft_fit()}.}
\item{n_models}{Positive integer number of candidate configurations.}
\item{top_k}{Positive integer maximum number of successful fits to retain.}
\item{n_testing}{Positive integer number of internal backtest origins.}
\item{search_space}{Named list of search specifications, or \code{NULL} for
the default space. An empty list repeats the fixed configuration with different
model seeds. See Details.}
\item{score_metric}{Backtest metric used for ranking. Lower scores are better.}
\item{seed}{Non-negative integer seed. Candidate \code{i} uses
\code{(seed + i) \%\% .Machine$integer.max}; this also controls its backtests.}
\item{fixed_args}{Named list of arguments to \code{craft_fit()}. Sampled
arguments override these, and these override arguments in \code{...}.}
\item{composite_weights}{Optional named numeric vector or named list of finite
non-negative component weights. Unspecified components have weight one.
At least one supplied weight must be positive.}
\item{keep_all_fits}{Whether to retain all successful fits in \code{all_fits}.}
\item{search_parallel}{Whether to evaluate candidates on a socket cluster.}
\item{search_workers}{Positive integer worker limit; \code{NULL} uses two
workers, capped by \code{n_models}.}
\item{verbose}{Whether to emit search progress messages.}
\item{...}{Additional named arguments to \code{craft_fit()}.}
}
\details{
Search names must match \code{craft_fit()} arguments. The arguments
\code{series}, \code{volumes}, \code{seed}, and \code{n_testing} are controlled
by the tuner and cannot appear in search specifications, fixed arguments or
\code{...}. Specifications may be scalar constants, atomic candidate vectors,
zero-argument sampler functions, or typed lists:
\describe{
\item{\code{list(type = "choice", values = c(...))}}{Uniform sampling from
non-empty atomic values.}
\item{\code{list(type = "integer", min = a, max = b)}}{Uniform sampling of
integers between inclusive integer bounds.}
\item{\code{list(type = "numeric", min = a, max = b, log = FALSE)}}{Uniform
numeric sampling, or uniform sampling on the log scale when \code{log = TRUE}.
Log-scale bounds must be positive.}
}
Other lists are passed as constant list-valued arguments. Sampler functions
should use R's random-number generator and avoid external state.

The default space searches window lengths 10, 15, 20, 30; 2--5 factors;
variance thresholds 0.85--0.98; fixed or variance factor selection; 3--6 regimes
per factor; minimum segments 10, 15, 20, 30, 40; log-uniform transition smoothing
0.25--5; transition minimum counts 1--5; uniform, weighted balanced accuracy,
or balanced accuracy sampler weights; and tail-penalty strength, threshold,
power, and aggregation. It requires a sufficiently long series; use a smaller
custom space for short histories.

Error metrics are minimized directly. Direction accuracy is negated, and
coverage scoring is the absolute gap from 0.90. Composite scoring takes a
weighted mean of normalized ranks for \code{return_rmse}, \code{return_mae},
\code{level_rmse}, \code{level_mae}, \code{coverage_90_gap},
\code{direction_error} (one minus accuracy), and \code{mean_p_actual_gap}
(absolute gap from 0.50). Ties use average ranks; constant components have rank
zero. Unavailable components are omitted and weights are renormalized over
available components. Composite scores are relative to the candidate pool.

Backtest fits use only observations through their origin. Retained candidate
fits use the full supplied history. Failed fits, failed backtest origins,
empty summaries, and non-finite selected scores are recorded with infinite
scores and excluded from selection. If all candidates fail, a warning is issued
and \code{top_models} is empty. Candidate drawing errors stop the search.

Search is serial by default, with serial internal backtests unless explicitly
enabled. Parallel search disables nested backtest parallelism, even if enabled
in fixed arguments; \code{backtest_parallel} cannot be tuned in parallel search.
The installed CRAFT namespace is loaded automatically on socket workers.

When tuning windows or minimum segment lengths, fix \code{horizon} and
\code{backtest_min_train} so scores compare the same origins and outcomes.
Tuning scores are selection statistics; use a separate later holdout period
to assess the selected model's predictive performance.
}
\value{
A list of class \code{craft_autotune} containing \code{call},
\code{score_metric}, \code{search_space}, \code{results} (candidate-order
scores, metrics, sampled parameters, and errors), \code{top_models} (up to
\code{top_k} entries ordered by increasing score, each containing
\code{model_id}, \code{params}, \code{score}, and a \code{craft_fit} object),
\code{failures}, \code{all_fits} (\code{NULL} unless requested, otherwise a
candidate-order list with \code{NULL} for failures), and \code{settings}.
Ties retain candidate order. The retained fits can be passed to
\code{craft_predict()}.
}
\seealso{\code{\link{craft_fit}}, \code{\link{craft_predict}}}
\examples{
set.seed(1)
series <- data.frame(a = cumprod(1 + rnorm(80, 0.001, 0.01)),
                     b = cumprod(1 + rnorm(80, 0.0005, 0.012)))
tuned <- craft_autotune(
  series, n_models = 2, top_k = 1, n_testing = 1,
  search_space = list(transition_alpha = c(0.5, 1, 2)),
  fixed_args = list(window = 5, n_draws = 30, n_factors = 1,
                   min_segment = 5, max_regimes_per_factor = 3),
  verbose = FALSE
)
tuned$results[, c("model_id", "ok", "score")]
if (length(tuned$top_models)) {
  prediction <- craft_predict(tuned$top_models[[1]]$fit, n_draws = 20)
}
}
