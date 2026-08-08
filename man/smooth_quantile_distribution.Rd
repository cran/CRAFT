\name{smooth_quantile_distribution}
\alias{smooth_quantile_distribution}
\title{Fit a Smooth Empirical Distribution from Quantiles}
\usage{
smooth_quantile_distribution(
  x,
  probs = seq(0.001, 0.999, length.out = 1000),
  spar = NULL,
  density_spar = NULL,
  eps = 1e-06,
  enforce_monotone = TRUE,
  plot_result = FALSE
)
}
\arguments{
\item{x}{Numeric vector of observations. Non-finite values are removed.}
\item{probs}{Probability grid strictly inside \code{(0, 1)}.}
\item{spar}{Optional smoothing parameter passed to \code{stats::smooth.spline()} for the quantile function.}
\item{density_spar}{Optional smoothing parameter passed to \code{stats::smooth.spline()} for the implied density.}
\item{eps}{Small positive lower bound used in derivative and probability calculations.}
\item{enforce_monotone}{Logical. If \code{TRUE}, enforce monotonicity of the smoothed quantile curve.}
\item{plot_result}{Logical. If \code{TRUE}, plot the fitted quantile curve and implied density.}
}
\description{
Fits a smooth empirical quantile function and derives an interpolated
density, cumulative distribution, quantile function, and random sampler.
}
\value{
A list with \code{data}, \code{dfun}, \code{pfun}, \code{qfun},
\code{rfun}, \code{q_fit}, \code{d_fit}, and \code{original_x}.
}
\examples{
set.seed(1)
fit <- smooth_quantile_distribution(rnorm(100), probs = seq(0.05, 0.95, length.out = 50))
fit$qfun(c(0.25, 0.5, 0.75))
fit$dfun(0)
}
