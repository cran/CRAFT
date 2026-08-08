\name{regime_sampler}
\alias{regime_sampler}
\title{Sample Candidate Rows from Regime Probabilities}
\usage{
regime_sampler(
  trafo_df,
  cluster_df,
  prob_vectors,
  weights = NULL,
  size = 1000L,
  return_type = c("stack", "index"),
  seed = NULL,
  eps = 1e-12,
  ...
)
}
\arguments{
\item{trafo_df}{Data frame of candidate transformed observations to sample.}
\item{cluster_df}{Data frame of regime labels for each candidate row. Must have the same number of rows as \code{trafo_df}.}
\item{prob_vectors}{Named list of probability vectors or one-row matrices/data frames, one per target regime column. Names should correspond to columns of \code{cluster_df}.}
\item{weights}{Optional numeric weights for target regime columns.}
\item{size}{Number of rows or indices to sample with replacement.}
\item{return_type}{Return sampled rows as \code{"stack"} or sampled row indices as \code{"index"}.}
\item{seed}{Optional random seed.}
\item{eps}{Small positive probability floor.}
\item{...}{Additional arguments. The legacy argument name \code{return} is accepted as an alias for \code{return_type}.}
}
\description{
Combines per-target regime probability vectors into row-level sampling
probabilities and samples candidate future-trajectory rows accordingly.
}
\value{
If \code{return_type = "stack"}, a data frame sampled from \code{trafo_df}.
If \code{return_type = "index"}, an integer vector of sampled row indices.
}
\examples{
trafo <- data.frame(x = 1:4)
clusters <- data.frame(r = c("low", "low", "high", "high"))
probs <- list(r = c(low = 0.2, high = 0.8))
regime_sampler(trafo, clusters, probs, size = 5, seed = 1)
}
