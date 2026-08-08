\name{trajectory_forecast}
\alias{trajectory_forecast}
\title{Fit Smooth Forecast Distributions for Trajectory Draws}
\usage{
trajectory_forecast(
  trajectory_draws,
  last_level = NULL,
  recon_mode = NULL,
  probs = seq(0.001, 0.999, length.out = 1000),
  spar = NULL,
  density_spar = NULL,
  eps = 1e-06,
  enforce_monotone = TRUE,
  plot_result = FALSE,
  min_obs = 10L,
  min_unique = 10L,
  jitter_if_needed = TRUE,
  jitter_scale = 1e-08,
  zero_density_outside = TRUE,
  seed = 1L,
  verbose = FALSE,
  return_diagnostics = TRUE
)
}
\arguments{
\item{trajectory_draws}{Data frame of sampled future-trajectory draws, with one variable per column.}
\item{last_level}{Optional scalar or vector of latest observed levels used to reconstruct forecast levels.}
\item{recon_mode}{Optional reconstruction mode: \code{"differences"}, \code{"scaled_growth"}, \code{"log_growth"}, or \code{NULL}.}
\item{probs}{Probability grid passed to \code{smooth_quantile_distribution()}.}
\item{spar}{Optional quantile smoothing parameter.}
\item{density_spar}{Optional density smoothing parameter.}
\item{eps}{Small positive numerical tolerance.}
\item{enforce_monotone}{Logical. If \code{TRUE}, enforce monotone quantile fits.}
\item{plot_result}{Logical. If \code{TRUE}, plot fitted distributions.}
\item{min_obs}{Minimum observations used for each fitted distribution. Shorter samples are bootstrapped.}
\item{min_unique}{Minimum unique values before jittering is considered.}
\item{jitter_if_needed}{Logical. If \code{TRUE}, add small jitter to discrete or constant samples.}
\item{jitter_scale}{Scale of jitter relative to the sample scale.}
\item{zero_density_outside}{Logical. If \code{TRUE}, return zero density outside the fitted support.}
\item{seed}{Optional random seed used for augmentation and jitter.}
\item{verbose}{Logical. If \code{TRUE}, emit fitting progress messages.}
\item{return_diagnostics}{Logical. If \code{TRUE}, attach diagnostics as an attribute.}
}
\description{
Fits one smoothed empirical forecast distribution per column in a data
frame of sampled future-trajectory draws. Optionally reconstructs returns or
differences into forecast levels.
}
\value{
A named list of distribution objects. Diagnostics and settings are attached
as attributes when requested.
}
\examples{
set.seed(1)
trajectory_draws <- data.frame(a = rnorm(40, 0.01, 0.02), b = rnorm(40, -0.01, 0.03))
fc <- trajectory_forecast(
  trajectory_draws,
  probs = seq(0.05, 0.95, length.out = 30),
  verbose = FALSE
)
fc$a$qfun(0.5)
}
