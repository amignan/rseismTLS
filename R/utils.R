#' Rejection sampling
#'
#' Estimates the credible interval of a specific probability distribution
#'
#' @param distr.x a vector of x increments
#' @param distr.dens a vector of density values for each x
#' @param interval the interval size (optional, default =.9 for 90% interval)
#' @param nsamples the number of samples used for rejection sampling (optional, default =1e6)
#' @return the credible interval
rejection_sampling <- function(distr.x, distr.dens, interval = .9, nsamples = 1e6){
  x.range <- range(distr.x)
  dx <- unique(diff(distr.x))[1]
  x <- x.range[1] + runif(nsamples) * diff(x.range)
  y <- runif(nsamples)
  require(signal)   #interp1()
  pr_x <- interp1(distr.x, distr.dens * dx, x)
  accepted_samples <- x[y < pr_x]
  credible_interval <- quantile(x = accepted_samples, probs = c(.5 - interval / 2, .5 + interval / 2))
  return(credible_interval)
}
