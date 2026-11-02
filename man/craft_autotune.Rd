\name{craft_autotune}
\alias{craft_autotune}
\title{Accelerated Staged Hyperparameter Tuning for CRAFT}
\description{
Random search with scoring-only backtests, staged screening, bounded shared
computation caches, dynamically scheduled parallel workers, and deferred full
fitting of selected configurations.
}
\usage{
craft_autotune(
  series, volumes = NULL, n_models = 30L, top_k = 5L, n_testing = 5L,
  search_space = NULL,
  score_metric = c("return_rmse", "return_mae", "level_rmse", "level_mae",
    "coverage_90_gap", "direction_accuracy", "composite"),
  seed = 42, fixed_args = list(), composite_weights = NULL,
  keep_all_fits = FALSE, search_parallel = FALSE, search_workers = NULL,
  stage_budgets = NULL, survivor_fraction = 0.3,
  cache = TRUE, cache_size = 64L, verbose = TRUE, ...
)
}
\arguments{
\item{series}{Numeric matrix or data frame with observations in rows.}
\item{volumes}{Optional volume matrix or data frame matching \code{series}.}
\item{n_models}{Positive integer number of candidate configurations.}
\item{top_k}{Positive integer maximum number of complete fitted models to select.}
\item{n_testing}{Positive integer requested final number of backtest origins.}
\item{search_space}{Named list of scalar constants, atomic candidate vectors,
zero-argument sampler functions, or typed specifications described in Details.
\code{NULL} uses the default search space. An empty list repeats the fixed configuration
with different model seeds.}
\item{score_metric}{Metric to minimize. Direction accuracy is negated; coverage
uses the absolute gap from 0.90. Composite scoring uses weighted normalized ranks.}
\item{seed}{Non-negative integer seed. Candidate seeds use
\code{(seed + model_id) \%\% .Machine$integer.max}; origin seeds add the stable
point ID modulo the same bound.}
\item{fixed_args}{Named arguments to \code{craft_fit()}. Sampled arguments override
these, which override \code{...}.}
\item{composite_weights}{Optional named non-negative finite weights for
\code{return_rmse}, \code{return_mae}, \code{level_rmse}, \code{level_mae},
\code{coverage_90_gap}, \code{direction_error}, and \code{mean_p_actual_gap}.
Unspecified weights are one; at least one supplied weight must be positive.}
\item{keep_all_fits}{If \code{FALSE}, fit only the selected finalists in full.
If \code{TRUE}, additionally refit every successfully scored candidate, including
screened candidates; this can substantially increase runtime and memory use.}
\item{search_parallel}{Whether to score candidates using socket workers.}
\item{search_workers}{Positive integer worker limit; \code{NULL} uses two,
capped by the number of feasible candidates.}
\item{stage_budgets}{Increasing positive integer vector of cumulative origin
budgets ending at \code{n_testing}. \code{NULL} uses
\code{unique(c(min(2L, n_testing), n_testing))}. Set \code{stage_budgets = n_testing}
to evaluate every candidate at the full budget without screening. Budgets are
capped by the number of distinct feasible origins.}
\item{survivor_fraction}{Fraction in \code{(0, 1]} promoted after each non-final
stage. At least \code{top_k} successful candidates are retained when available.}
\item{cache}{Whether to reuse embeddings and deterministic regime models.}
\item{cache_size}{Positive integer maximum number of least-recently-used cache
entries per process, combining embeddings and regime models. This limits entries,
not bytes; lower this value for large histories or constrained memory.}
\item{verbose}{Whether to report stage progress.}
\item{...}{Additional named arguments to \code{craft_fit()}.}
}
\details{

Search specifications may be scalar constants, atomic candidate vectors,
zero-argument sampler functions, or typed lists:
\describe{
\item{\code{list(type = "choice", values = c(...))}}{Uniform sampling from
non-empty atomic values without missing entries.}
\item{\code{list(type = "integer", min = a, max = b)}}{Uniform integer
sampling between inclusive integer bounds.}
\item{\code{list(type = "numeric", min = a, max = b, log = FALSE)}}{Uniform
numeric sampling, or sampling uniformly on the log scale when \code{log = TRUE}.
Log-scale bounds must be positive.}
}
Other lists are constant list-valued arguments. Sampled settings override
\code{fixed_args}, which override \code{...}. To hold a parameter fixed, omit it
from the search space. An explicit search space replaces the default space.
The default space searches windows 10, 15, 20, 30; 2--5 factors; variance
thresholds 0.85--0.98; fixed or variance factor selection; 3--6 regimes per
factor; minimum segments 10, 15, 20, 30, 40; transition smoothing and minimum
counts; sampler weights; and tail-penalty settings.

Composite scoring averages weighted normalized ranks for error metrics,
absolute coverage gap from 0.90, direction error (one minus accuracy), and
absolute mean probability-calibration gap from 0.50. Unavailable components
are omitted and weights renormalized. Constant components have rank zero;
ties use average ranks. Composite scores depend on the candidate pool.

Scoring fits stop after generating and aligning sampled trajectories. They skip
smooth return and level distributions, correlation diagnostics, and full output
construction. Backtest metrics use empirical sampled trajectories. Full-history fits are deferred until selection; no
full-history model is needed to rank a candidate.

All feasible configurations share one origin grid. Its earliest origin is the
largest minimum training length required by any feasible candidate, and its
latest origin leaves room for every candidate's requested forecast horizon.
Configurations with invalid windows or no individually feasible origins are
recorded as failures. If individually feasible configurations have no common
origin grid, the search stops with an error. Fix \code{horizon} when tuning
windows to ensure that candidates also forecast the same outcomes. With a
\code{NULL} horizon, each candidate scores all horizons through its own window.

The first two origin points span the common grid, with further points added in
grid order. Each origin retains its point ID and seed across stages. Completed
origins are reused rather than fitted again. Each stage ranks only candidates
that reached that stage, and final selection considers only candidates that
completed the final budget. Screening can discard a strong configuration because
its early score is noisy. Composite scores are relative to the candidate pool
within each stage and cannot be compared across stages.

Caches are local to one search and one process. Keys distinguish training-prefix
length, trajectory window, side, forecast horizon, and regime arguments. Reusing
a regime model also reuses its denoising, SVD, and changepoint computations.
Different structural settings require different regime fits. Regime calls with
bootstrap stability, an explicit regime seed, or plotting bypass regime caching.
No regime is fitted on observations after its training origin.

Parallel search transfers the immutable input history to each worker once and
uses load-balanced candidate tasks. Workers return compact prediction errors and
calibration columns, not fitted objects or sampled trajectory matrices. Workers
and their caches persist across stages and are stopped before full refitting.
Cache hit counts and elapsed times can vary with scheduling; scores do not for
the same seeds and deterministic specifications. Candidate sampler functions run
on the controller and should avoid external state.

Nested backtest parallelism and training-probability output are disabled.
\code{backtest_parallel}, \code{backtest_workers}, \code{verbose}, and
\code{return_train_probs} cannot be tuned; their fixed values are overridden by
the fast execution policy. \code{series}, \code{volumes}, \code{seed}, and
\code{n_testing} must be supplied through the tuner's own arguments.

Selected complete fits use the original candidate seed and number of draws.
Their \code{backtest} field contains cached summaries and compact scoring
predictions for the origins actually evaluated. It omits unneeded prediction
columns and per-origin fitted objects. Refitting failures are recorded separately
and remaining finalists are tried until \code{top_k} fits are selected or the
finalists are exhausted. No usable final fit produces a warning and an empty
\code{top_models} list.

Set \code{stage_budgets = n_testing} to disable screening. For the same
origin grid, candidate configurations, and seeds, enabling or disabling the
deterministic cache preserves empirical backtest scores.
Use a later holdout period to evaluate the selected model; tuning scores are
selection statistics.
}
\value{
A list of class \code{craft_autotune}, with:
\describe{
\item{\code{results}}{Candidate-order scores and metrics, sampled parameters,
\code{stage}, evaluated \code{n_testing}, \code{status}, \code{ok},
\code{error}, and \code{selected}. Status distinguishes \code{screened},
\code{finalist}, \code{failed}, and \code{refit_failed}. A screened candidate's
score uses fewer origins and should not be ranked against finalists.}
\item{\code{top_models}}{Selected finalists ordered by score, each containing
\code{model_id}, \code{params}, \code{score}, and a complete \code{craft_fit}
object usable with \code{craft_predict()}.}
\item{\code{failures}, \code{refit_failures}}{Scoring/feasibility failures and
full-refitting failures, respectively.}
\item{\code{all_fits}}{\code{NULL} unless requested; otherwise a candidate-order
list of complete fits with \code{NULL} for unsuccessful refits.}
\item{\code{stages}}{Each stage's budget, origins, result table, and promoted
candidate IDs. Ties retain candidate order.}
\item{\code{performance}}{Per-candidate, per-stage origin-fit counts, elapsed
seconds, cache hits, and cache misses. Excludes final full-history refits.}
\item{\code{settings}}{Requested and actual budgets, common origin grid, seeds,
cache configuration, and worker count.}
\item{\code{call}, \code{score_metric}, \code{search_space}}{The matched call,
selected metric, and search specifications.}
}
}
\seealso{\code{\link{craft_fit}}, \code{\link{craft_predict}}}
\examples{
set.seed(1)
series <- data.frame(a = cumprod(1 + rnorm(100, .001, .01)),
                     b = cumprod(1 + rnorm(100, .0005, .012)))
tuned <- craft_autotune(
  series, n_models = 4, top_k = 1, n_testing = 3,
  search_space = list(transition_alpha = c(.5, 1, 2)),
  fixed_args = list(window = 5, horizon = 5, n_draws = 40,
                   n_factors = 1, min_segment = 5,
                   max_regimes_per_factor = 3),
  verbose = FALSE
)
tuned$results[, c("model_id", "score", "n_testing", "status", "selected")]
if (length(tuned$top_models)) {
  prediction <- craft_predict(tuned$top_models[[1]]$fit, n_draws = 20)
}
}
