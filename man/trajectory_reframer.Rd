\name{trajectory_reframer}
\alias{trajectory_reframer}
\title{Build Lag or Lead Trajectory Features for One Series}
\usage{
trajectory_reframer(ts, k, mode = c("backward", "forward"), ts_name = NULL)
}
\arguments{
\item{ts}{A numeric vector or time-series-like object.}
\item{k}{Positive integer number of lag or lead horizons to create.}
\item{mode}{Direction of the trajectory. \code{"backward"} creates cumulative lag changes; \code{"forward"} creates cumulative lead changes.}
\item{ts_name}{Optional name used as the prefix for output columns.}
}
\description{
Creates cumulative scaled lag or lead trajectory features for a single time
series. For a backward trajectory, horizon \code{h} at time \code{t} is
\code{x[t] / x[t - h] - 1}. For a forward trajectory, horizon \code{h} is
\code{x[t + h] / x[t] - 1}.
}
\value{
An object of class \code{trajectory_reframer}, represented as a list with
\code{trajectory}, \code{ts}, \code{k}, and \code{mode}.
}
\examples{
x <- c(100, 102, 105, 103, 108)
traj <- trajectory_reframer(x, k = 2, mode = "backward", ts_name = "price")
head(traj$trajectory)
}
