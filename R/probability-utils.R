.mc_normalize_prob <- function(p, eps = 1e-12) {
  p <- as.numeric(p)

  if (!length(p)) {
    return(p)
  }

  p[!is.finite(p) | p < 0] <- 0

  if (sum(p) <= eps) {
    p[] <- 1 / length(p)
  } else {
    p <- p / sum(p)
  }

  p
}

.mc_softmax <- function(x) {
  x <- as.numeric(x)

  if (!length(x)) {
    return(x)
  }

  x[!is.finite(x)] <- -Inf

  if (all(!is.finite(x))) {
    return(rep(1 / length(x), length(x)))
  }

  z <- exp(x - max(x, na.rm = TRUE))

  if (sum(z) <= 0 || any(!is.finite(z))) {
    return(rep(1 / length(x), length(x)))
  }

  z / sum(z)
}
