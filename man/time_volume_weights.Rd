\name{time_volume_weights}
\alias{time_volume_weights}
\title{Compute Time Weights from Trading Volumes}
\usage{
time_volume_weights(
  volume_df,
  date_col = NULL,
  alpha = 1,
  normalize = TRUE,
  eps = 1e-08,
  winsorize = FALSE,
  winsor_q = 0.99
)
}
\arguments{
\item{volume_df}{A data frame containing one or more numeric volume columns.}
\item{date_col}{Optional name of a date column to exclude from the numeric volume aggregation and use as output names.}
\item{alpha}{Non-negative exponent applied to aggregate volume scores. Values below 1 compress high-volume observations; values above 1 emphasize them.}
\item{normalize}{Logical. If \code{TRUE}, return weights summing to one.}
\item{eps}{Small positive offset added before exponentiation.}
\item{winsorize}{Logical. If \code{TRUE}, cap aggregate volumes at \code{winsor_q}.}
\item{winsor_q}{Quantile used as the cap when \code{winsorize = TRUE}.}
}
\description{
Aggregates numeric volume columns row-wise and converts the aggregate
volume profile into non-negative time weights.
}
\value{A numeric vector of weights, named by \code{date_col} when supplied.}
\examples{
vol <- data.frame(day = as.Date("2024-01-01") + 0:2, a = c(10, 20, 30), b = c(5, 5, 10))
time_volume_weights(vol, date_col = "day")
}
