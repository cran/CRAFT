\name{craft_predict}
\alias{craft_predict}
\title{Predict from a Fitted CRAFT Model}
\usage{
craft_predict(
  object,
  newdata = NULL,
  n_draws = NULL,
  seed = NULL,
  type = c("all", "trajectory_draws", "distributions", "probabilities", "labels")
)
}
\arguments{
\item{object}{A fitted CRAFT object returned by \code{craft_fit()}.}
\item{newdata}{Optional matrix or data frame containing the most recent series history. If supplied, it must contain the fitted asset columns and more rows than the fitted trajectory window. If \code{NULL}, the latest trajectory stored in \code{object} is reused.}
\item{n_draws}{Optional number of sampled future-trajectory draws. Defaults to the number of draws stored in \code{object}.}
\item{seed}{Optional random seed for sampling.}
\item{type}{Prediction output type. \code{"all"} returns the full prediction object; \code{"trajectory_draws"} returns sampled future-trajectory rows; \code{"distributions"} returns return and level forecast distributions; \code{"probabilities"} returns transition probabilities; \code{"labels"} returns latest regime labels.}
}
\description{
Uses a fitted CRAFT model to transform the latest trajectory, estimate future
regime probabilities, sample future-trajectory analogues, and fit smooth
forecast distributions.
}
\value{
For \code{type = "all"}, an object of class \code{craft_prediction}, a list
with \code{latest_regime_labels}, \code{latest_regime_probabilities}, \code{trajectory_draws},
\code{sampler_trace}, \code{corr_by_horizon}, \code{return_dists},
\code{level_dists}, and \code{diagnostics}. Other \code{type} values return
the requested subset.
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
pred <- craft_predict(fit, n_draws = 25, seed = 1)
names(pred$return_dists)
}
}
