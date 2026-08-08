\name{denoiser}
\alias{denoiser}
\title{Select Informative Numeric Features}
\usage{
denoiser(
  x,
  target = NULL,
  min_non_missing = 0.7,
  min_unique = 3L,
  min_mad = 1e-08,
  max_cor = 0.95,
  keep_quantile = 0.7,
  min_score = NULL,
  min_features = 1L,
  max_features = NULL,
  top_n = NULL,
  supervised_weight = if (!is.null(target)) 0.6 else 0,
  pca_weight = 0.45,
  temporal_weight = 0.25,
  snr_weight = 0.15,
  info_weight = 0.15,
  lag_max = 5L,
  impute = c("median", "zero"),
  return_report = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{x}{A matrix or data frame of candidate features.}
\item{target}{Optional target vector, or the name of a column in \code{x}, used for supervised scoring.}
\item{min_non_missing}{Minimum non-missing fraction required before a feature can be selected.}
\item{min_unique}{Minimum number of unique finite values required.}
\item{min_mad}{Minimum median absolute deviation required.}
\item{max_cor}{Maximum absolute pairwise correlation allowed among selected features.}
\item{keep_quantile}{Score quantile used as an adaptive selection threshold.}
\item{min_score}{Optional absolute minimum score threshold.}
\item{min_features}{Minimum number of features to retain when possible.}
\item{max_features}{Optional maximum number of retained features.}
\item{top_n}{Optional number of top-scoring features to retain.}
\item{supervised_weight}{Weight assigned to supervised target association.}
\item{pca_weight}{Weight assigned to PCA loading contribution.}
\item{temporal_weight}{Weight assigned to temporal autocorrelation structure.}
\item{snr_weight}{Weight assigned to signal-to-noise proxy scoring.}
\item{info_weight}{Weight assigned to information and completeness scoring.}
\item{lag_max}{Maximum lag used for temporal-structure scoring.}
\item{impute}{Imputation method used during scoring.}
\item{return_report}{Logical. If \code{TRUE}, return selected data, feature names, report, and settings.}
\item{verbose}{Logical. If \code{TRUE}, emit a summary message.}
}
\description{
Ranks numeric features using completeness, variation, temporal structure,
PCA contribution, signal-to-noise, and optional supervised association,
then returns a reduced feature set.
}
\value{
By default, a matrix or data frame with selected features and attributes
\code{denoiser_report}, \code{selected_features}, and \code{removed_features}.
If \code{return_report = TRUE}, a list with \code{data}, \code{features},
\code{report}, and \code{settings}.
}
\examples{
set.seed(1)
x <- data.frame(signal = cumsum(rnorm(30)), noise = rnorm(30), const = 1)
out <- denoiser(x, min_features = 1, verbose = FALSE, return_report = TRUE)
out$features
}
