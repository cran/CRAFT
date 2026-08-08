\name{svd_changepoint_regimes}
\alias{svd_changepoint_regimes}
\title{Detect Regimes from SVD Scores and Changepoints}
\usage{
svd_changepoint_regimes(
  X,
  series_cols = NULL,
  denoise = TRUE,
  denoiser_args = NULL,
  reduction_criteria = c(
    "variance", "elbow", "kaiser", "bic", "fixed",
    "regime_strength", "hybrid"
  ),
  n_dim = NULL,
  var_threshold = 0.95,
  max_dim = NULL,
  center = TRUE,
  scale = TRUE,
  na_action = c("impute", "stop"),
  cpt_method = "PELT",
  cpt_penalty = "MBIC",
  minseglen = 5,
  filter_weak_dims = FALSE,
  min_regime_strength = 0.1,
  max_selected_dims = NULL,
  max_regimes_per_dim = 6,
  min_regime_size = 30,
  merge_similar_regimes = TRUE,
  merge_threshold = 0.35,
  stability_bootstrap = FALSE,
  n_bootstrap = 100,
  bootstrap_noise = 0.05,
  cpt_tolerance = NULL,
  min_cpt_stability = 0.6,
  newdata_temperature = 1,
  regime_prefix = "regime",
  return_scores = TRUE,
  return_probabilities = TRUE,
  plot = FALSE,
  plot_type = c("scores", "heatmap", "summary"),
  seed = NULL
)
}
\arguments{
\item{X}{A data frame or matrix containing numeric time-series columns.}
\item{series_cols}{Optional names of columns to use. If \code{NULL}, all numeric columns are used.}
\item{denoise}{Logical. If \code{TRUE}, select features with \code{denoiser()} before fitting.}
\item{denoiser_args}{Optional list of arguments passed to \code{denoiser()}.}
\item{reduction_criteria}{Criterion used to select SVD dimensions.}
\item{n_dim}{Number of dimensions when \code{reduction_criteria = "fixed"}.}
\item{var_threshold}{Cumulative explained variance threshold for variance-based selection.}
\item{max_dim}{Optional maximum SVD dimension considered.}
\item{center}{Logical. If \code{TRUE}, center selected columns.}
\item{scale}{Logical. If \code{TRUE}, scale selected columns.}
\item{na_action}{How to handle missing values: impute column means or stop.}
\item{cpt_method}{Changepoint method passed to \code{changepoint::cpt.meanvar()}.}
\item{cpt_penalty}{Penalty passed to \code{changepoint::cpt.meanvar()}.}
\item{minseglen}{Minimum segment length for changepoint fitting.}
\item{filter_weak_dims}{Logical. If \code{TRUE}, filter dimensions with weak regime strength.}
\item{min_regime_strength}{Minimum regime-strength threshold.}
\item{max_selected_dims}{Optional maximum number of selected dimensions.}
\item{max_regimes_per_dim}{Maximum regimes retained per SVD dimension after merging.}
\item{min_regime_size}{Minimum desired regime segment size.}
\item{merge_similar_regimes}{Logical. If \code{TRUE}, merge short or similar adjacent regimes.}
\item{merge_threshold}{Distance threshold used when merging similar regimes.}
\item{stability_bootstrap}{Logical. If \code{TRUE}, assess changepoint stability by bootstrapping.}
\item{n_bootstrap}{Number of bootstrap replicates.}
\item{bootstrap_noise}{Noise scale used for bootstrap perturbations.}
\item{cpt_tolerance}{Maximum distance for matching bootstrapped changepoints.}
\item{min_cpt_stability}{Minimum bootstrap stability for retaining changepoints.}
\item{newdata_temperature}{Temperature used when transforming new observations to regime probabilities.}
\item{regime_prefix}{Prefix used for generated regime-label columns.}
\item{return_scores}{Logical. If \code{TRUE}, include selected SVD scores.}
\item{return_probabilities}{Logical. If \code{TRUE}, include regime probabilities.}
\item{plot}{Logical. If \code{TRUE}, produce a plot.}
\item{plot_type}{Plot type for \code{plot = TRUE}.}
\item{seed}{Optional random seed.}
}
\description{
Projects multivariate time-series columns into SVD scores, detects
changepoints in selected score dimensions, merges or filters regimes, and
returns regime labels, profiles, diagnostics, and helper functions for new
data transformation and plotting.
}
\value{
An object of class \code{svd_changepoint_regimes_v2}, represented as a list
with regime labels, selected dimensions, diagnostics, profiles, SVD
metadata, changepoint results, and helper functions.
}
\examples{
\donttest{
set.seed(1)
x <- data.frame(a = c(rnorm(40, 0), rnorm(40, 2)), b = c(rnorm(40, 0), rnorm(40, -1)))
fit <- svd_changepoint_regimes(
  x,
  denoise = FALSE,
  reduction_criteria = "fixed",
  n_dim = 1,
  minseglen = 5,
  min_regime_size = 5,
  max_regimes_per_dim = 3,
  seed = 1
)
head(fit$regimes)
}
}
