smooth_quantile_distribution <- function(
  x,
  probs = seq(0.001, 0.999, length.out = 1000),
  spar = NULL,
  density_spar = NULL,
  eps = 1e-6,
  enforce_monotone = TRUE,
  plot_result = FALSE
) {
  ####################################################
  #######SUPPORT FUN##################################
  gradient_quantile <- function(x, p) {
  n <- length(x)
  
  if (n < 3) {
    stop("Need at least 3 points to compute gradient.")
  }
  
  g <- numeric(n)
  
  # Q'(p), where x = Q(p)
  g[1] <- (x[2] - x[1]) / (p[2] - p[1])
  g[n] <- (x[n] - x[n - 1]) / (p[n] - p[n - 1])
  
  for (i in 2:(n - 1)) {
    g[i] <- (x[i + 1] - x[i - 1]) / (p[i + 1] - p[i - 1])
  }
  
  pmax(g, .Machine$double.eps)
}


trapz_simple <- function(x, y) {
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}


cumulative_trapz <- function(x, y) {
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  
  out <- numeric(length(x))
  
  if (length(x) > 1) {
    out[-1] <- cumsum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }
  
  out
}  

  ####################################################
  x <- x[is.finite(x)]
  
  if (length(x) < 10) {
    stop("Need at least 10 finite observations.")
  }
  
  if (any(probs <= 0 | probs >= 1)) {
    stop("probs must be strictly inside (0, 1).")
  }
  
  probs <- sort(unique(probs))
  
  # Empirical quantiles
  q_emp <- as.numeric(stats::quantile(x, probs = probs, type = 8, names = FALSE))
  
  # Smooth quantile function Q(p)
  q_fit <- stats::smooth.spline(
    x = probs,
    y = q_emp,
    spar = spar
  )
  
  q_smooth <- as.numeric(stats::predict(q_fit, probs)$y)
  
  # Defensive monotonicity correction
  if (enforce_monotone) {
    q_smooth <- cummax(q_smooth)
  }
  
  # Remove duplicated quantile values, otherwise inverse interpolation can fail
  keep <- !duplicated(q_smooth)
  
  probs_u <- probs[keep]
  q_u <- q_smooth[keep]
  
  if (length(q_u) < 10) {
    stop("Smoothed quantile function has too few unique values. Try lower smoothing.")
  }
  
  # Estimate density from derivative of quantile function:
  # f(Q(p)) = 1 / Q'(p)
  dq_dp <- gradient_quantile(q_u, probs_u)
  density_u <- 1 / pmax(dq_dp, eps)
  
  # Smooth density over x-grid
  d_fit <- stats::smooth.spline(
    x = q_u,
    y = density_u,
    spar = density_spar
  )
  
  density_smooth <- as.numeric(stats::predict(d_fit, q_u)$y)
  density_smooth <- pmax(density_smooth, 0)
  
  # Normalize density numerically
  area <- trapz_simple(q_u, density_smooth)
  
  if (!is.finite(area) || area <= 0) {
    stop("Density normalization failed. Try different smoothing parameters.")
  }
  
  density_smooth <- density_smooth / area
  
  # Reconstruct smooth CDF from normalized density
  cdf_u <- cumulative_trapz(q_u, density_smooth)
  cdf_u <- cdf_u / max(cdf_u)
  
  # Ensure valid boundaries
  cdf_u <- pmin(pmax(cdf_u, 0), 1)
  
  # Main output table
  out <- data.frame(
    p = probs_u,
    x = q_u,
    q_empirical = q_emp[keep],
    q_smooth = q_u,
    density = density_smooth,
    cdf = cdf_u
  )
  
  # Density function
  dfun <- function(x_new) {
    stats::approx(
      x = out$x,
      y = out$density,
      xout = x_new,
      rule = 2,
      ties = "ordered"
    )$y
  }
  
  # CDF function
  pfun <- function(q_new) {
    stats::approx(
      x = out$x,
      y = out$cdf,
      xout = q_new,
      rule = 2,
      ties = "ordered"
    )$y
  }
  
  # Quantile function
  qfun <- function(p_new) {
    p_new <- pmin(pmax(p_new, min(out$cdf)), max(out$cdf))
    
    stats::approx(
      x = out$cdf,
      y = out$x,
      xout = p_new,
      rule = 2,
      ties = "ordered"
    )$y
  }
  
  # Random generator
  rfun <- function(n) {
    if (length(n) != 1 || n <= 0) {
      stop("n must be a positive integer.")
    }
    
    u <- stats::runif(n)
    qfun(u)
  }
  
  if (plot_result) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(mfrow = c(1, 2))
    
    plot(
      probs,
      q_emp,
      pch = 16,
      cex = 0.45,
      main = "Smoothed quantile function",
      xlab = "Probability",
      ylab = "Quantile"
    )
    
    lines(
      out$p,
      out$q_smooth,
      lwd = 2
    )
    
    plot(
      out$x,
      out$density,
      type = "l",
      lwd = 2,
      main = "Implied smoothed density",
      xlab = "x",
      ylab = "Density"
    )
    
    rug(x)
  }
  
  return(list(
    data = out,
    dfun = dfun,
    pfun = pfun,
    qfun = qfun,
    rfun = rfun,
    q_fit = q_fit,
    d_fit = d_fit,
    original_x = x
  ))
}


