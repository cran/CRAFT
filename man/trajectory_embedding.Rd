\name{trajectory_embedding}
\alias{trajectory_embedding}
\title{Build Past and Future Trajectory Matrices}
\usage{
trajectory_embedding(ts_set, trajectory_window)
}
\arguments{
\item{ts_set}{A data frame or matrix of numeric time series in columns.}
\item{trajectory_window}{Positive integer number of lag and lead horizons.}
}
\description{
Builds past and future cumulative-trajectory matrices for all columns in a
multivariate time-series data set. Rows containing unavailable lag or lead
values are removed from the aligned training matrices, while
\code{latest_trajectory} retains the latest past trajectory row for
forecasting.
}
\value{
A list with \code{past_trajectories}, \code{future_trajectories},
\code{latest_trajectory}, and \code{time_index}. The \code{time_index} vector
indicates rows retained in the aligned past and future trajectory matrices.
}
\examples{
prices <- data.frame(a = c(100, 101, 103, 104, 106), b = c(50, 52, 51, 53, 54))
trajectories <- trajectory_embedding(prices, trajectory_window = 2)
dim(trajectories$past_trajectories)
dim(trajectories$future_trajectories)
}
