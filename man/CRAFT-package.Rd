\docType{package}
\name{CRAFT-package}
\alias{CRAFT-package}
\alias{CRAFT}
\title{Conditional Regime Analog Forecasting with Trajectories}
\description{
CRAFT provides tools for building lag and lead trajectory embeddings from
multivariate time series, detecting regimes with singular value
decomposition and changepoint analysis, estimating regime-transition
probabilities, sampling future-trajectory analogues, and fitting smooth
empirical forecast distributions.
}
\details{
The main workflow is implemented by \code{craft_fit()} and \code{craft_predict()}.
Lower-level helpers are also exported for users who want to build custom
workflows: \code{trajectory_embedding()}, \code{svd_changepoint_regimes()},
\code{regime_sampler()}, \code{trajectory_forecast()},
\code{smooth_quantile_distribution()}, \code{denoiser()}, and
\code{time_volume_weights()}.
}
\author{Giancarlo C.}
\seealso{
\code{craft_fit()}, \code{craft_predict()}, \code{trajectory_embedding()},
\code{svd_changepoint_regimes()}
}
\keyword{package}
