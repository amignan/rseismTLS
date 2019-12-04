#' Complete Negative Log Likelihood Function
#'
#' Estimate the negative log likelihood of the statistical model of Mignan et al. (2017) for
#' a specific set of input parameters.
#'
#' See Eq. A2 of Broccardo et al. (2017). The function is here structured to be used by optim().
#'
#' @param data a list containing all the necessary data:
#' * `seism` an earthquake catalogue data frame of parameters `t` (time in dec. days) and `m` (magnitude)
#' * `inj` matching injection profile data frame of parameters `t` (time in dec. days), `dV` (flow rate in cubic
#' metre/day) and `V` (cumulative injected volume in cubic metres)
#' * `ts` shut-in time (in decimal days)
#' * `Tmax` upper range of the time window (in decimal days)
#' @param par a vector of the input parameters `a_fb`, `tau` and `b`
#' @return the negative log likelihood estimate for specific `par` values
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{model_par.val}
negloglik.val <- function(data, par){
  theta <- list(a_fb = par[1], tau = par[2], b = par[3])

  N <- nrow(data$seism)
  seism.inj <- subset(data$seism, data$seism$t <= data$ts)
  seism.post <- subset(data$seism, data$seism$t > data$ts)
  Ninj <- nrow(seism.inj)
  Npost <- nrow(seism.post)
  dV <- interp1(data$inj$t, data$inj$dV, seism.inj$t)
  dV.ts <- tail(dV, 1)
  V.ts <- tail(data$inj$V, 1)

  C1 <- N * (theta$a_fb - theta$b * data$m0) / log10(exp(1))
  C2 <- sum(log(dV), na.rm = T)
  C3 <- Npost * log(dV.ts)
  C4 <- -1 / theta$tau * (sum(seism.post$t - data$ts))
  C5 <- -10^(theta$a_fb - theta$b * data$m0) * (V.ts + dV.ts * theta$tau * (1-exp(-(data$Tmax - data$ts)/theta$tau)))
  C6 <- N * log(theta$b)
  C7 <- N * log(log(10))
  C8 <- -theta$b * log(10) * sum(data$seism$m)
  C9 <- N  * theta$b * log(10) * (data$m0 - mbin/2)   # Mmax -> +Inf
  nLL <- -(C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9)
  return(nLL)
}

#' Model Parameter Maximum Likelihood Estimation
#'
#' Provides the maximum likelihood estimates (MLE) of the parameters of the statistical model of Mignan et al. (2017)
#' using the complete negative log likelihood function \code{negloglik.val}.
#'
#' Optimization done using the optim() function of the stats package.
#'
#' @param data a list containing all the necessary data:
#' * `seism` an earthquake catalogue data frame of parameters `t` (time in dec. days) and `m` (magnitude)
#' * `inj` matching injection profile data frame of parameters `t` (time in dec. days), `dV` (flow rate in cubic
#' metre/day) and `V` (cumulative injected volume in cubic metres)
#' * `ts` shut-in time (in decimal days)
#' * `Tmax` upper range of the time window (in decimal days)
#' @param par a vector of the initial values for the parameters `a_fb`, `tau` and `b` to be optimized over
#' @return a list of the parameters' MLEs:
#' * `a_fb` the underground feedback activation (in /cubic metre)
#' * `tau` the mean relaxation time (in days)
#' * `b` the slope of the Gutenberg-Richter law
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{negloglik.val}
model_par.val <- function(data, par = c(0, 1, 1)) {
  require(signal)       # interp1()
  res <- optim(par = par, fn = rseismTLS::negloglik.val, data = data)
  return(list(a_fb = res$par[1], tau = res$par[2], b = res$par[3], nLL = res$value))
}

#' Data Binning
#'
#' Bucketize the two main data sets (injection `inj` and seismicity `seism`) in time intervals `tint`.
#'
#' Time `t` represents the centre of each bin of the time intervals `tint`.
#'
#' @param seism an earthquake catalogue data frame of parameters:
#' * `t` the occurrence time (in decimal days)
#' * `m` the magnitude
#' @param inj matching injection profile data frame of parameters:
#' * `t` the occurrence time (in decimal days)
#' * `dV` the injected volume (in cubic metres)
#' * `V` the cumulative injected volume (in cubic metres)
#' @param tint time intervals of constant length tbin (in decimal days)
#' @return A list of 2 data frames binned in time:
#' * `seism.binned` with parameters `t` and `rate` (number of events per `tbin`)
#' * `inj.binned` with parameters `t`, `dV` (in cubic metres per `tbin`) and `V`
data.bin <- function(seism, inj, tint) {
  seism.binned <- hist(seism$t, breaks = tint, plot = F)

  dt <- unique(diff(tint))[1]
  require(signal)       # interp1()
  inj.binned <- data.frame(t = seism.binned$mids,
                           dV = interp1(c(0, inj$t), c(0, inj$dV), seism.binned$mids) * dt,
                           V = interp1(c(0, inj$t), c(0, inj$V), seism.binned$mids))

  return(list(seism.rate = data.frame(t = seism.binned$mids, rate = seism.binned$counts),
              inj.binned = subset(inj.binned, !is.na(dV))))
}

#' Underground feedback activation
#'
#' Estimate the underground feedback activation \out{<i>a<sub>fb</sub></i>}, which is the
#' \out{<i>a</i>}-value of the Gutenberg-Richter law (Gutenberg and Richter, 1944)
#' normalized by the injected volume (e.g., Dinske and Shapiro, 2013).
#'
#' This is equivalent to Shapiro's Seismogenic Index (e.g., Dinske and Shapiro, 2013)
#' \out{<i>SI = log<sub>10</sub>(N<sub>tot</sub> / V<sub>tot</sub>) + b * m<sub>c</sub></i>}
#' but here theoretically agnostic, following the notation of Mignan et al. (2017).
#'
#' @param Ntot total number of events above `mc` for `Vtot`
#' @param mc the completeness magnitude
#' @param b the slope of the Gutenberg-Richter law
#' @param Vtot total volume of fluids injected
#' @return The numeric value of the underground feedback activation
#' @references Dinske C., Shapiro S.A. (2013), Seismotectonic state of reservoirs inferred
#' from magnitude distributions of fluid-induced seismicity. J. Seismol., 17, 13-25
#' \href{https://link.springer.com/article/10.1007/s10950-012-9292-9}{doi: 10.1007/s10950-012-9292-9}
#' @references Gutenberg B., Richter C.F. (1944), Frequency of earthquakes in California.
#' Bull. Seismol. Soc. Am.,
#' \href{https://authors.library.caltech.edu/47734/1/185.full.pdf}{34, 184-188}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#'\href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
a_fb.val <- function(Ntot, b, mc, Vtot) {
  log10(Ntot / Vtot) + b * mc
}

#' Statistical model of induced seismicity
#'
#' Statistical model proposed by Mignan et al. (2017) that predicts the
#' co-injection induced seismity rate via the linear relationship with flow rate
#' and the post-injection rate via exponential decay.
#'
#' The linear relationship is in agreement with the literature (Dinske and Shapiro,
#' 2013; Mignan, 2016; van der Elst et al., 2016). The exponential model was
#' verified to perform best when tested on 6 stimulations (Mignan et al., 2017).
#'
#' The flow rate is here defined as `dV`, the volume injected in the period (t-dt, dt) with `dt`
#' defined in the `data.bin()` function. The seismicity rate is modelled for the same time vector
#' as `inj` for co-injection and uses `t.postinj` for post-injection.
#'
#' @param method the method to be used: "`co-injection`", "`post-injection`", or "`full sequence`"
#' @param theta the list of model parameters:
#' * `a_fb` the underground feedback activation (`NULL` for "`post-injection`")
#' * `tau` the mean relaxation time (`NULL` for "`co-injection`")
#' * `b` the slope of the Gutenberg-Richter law
#' * `mc` the completeness magnitude
#' @param inj the binned injection profile data frame with parameters
#' * `t` the occurrence time (in decimal days)
#' * `dV` the flow rate (in cubic metres per bin)
#' @param shutin the list of shut-in parameters (required by "`post-injection`")
#' * `t` the shut-in time (in decimal days)
#' * `rate` the rate of induced seismicity at shut-in
#' @param t.postinj the time vector for which a seismicity rate is predicted
#' @return A data frame of the modelled seismicity rate with parameters:
#' * `t` the time (in decimal days)
#' * `rate` the seismicity rate per bin
#' @references Dinske C., Shapiro S.A. (2013), Seismotectonic state of reservoirs inferred
#' from magnitude distributions of fluid-induced seismicity. J. Seismol., 17, 13-25
#' \href{https://link.springer.com/article/10.1007/s10950-012-9292-9}{doi: 10.1007/s10950-012-9292-9}
#' @references Mignan A. (2016), Static behaviour of induced seismicity. Nonlin. Processes Geophys.,
#' 23, 107-113, \href{https://www.nonlin-processes-geophys.net/23/107/2016/}{doi: 10.5194/npg-23-107-2016}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @references van der Elst N.J., Page M.T., Weiser D.A., Goebel T.H.W., Hosseini S.M. (2016),
#' Induced earthquake magnitudes are as large as (statistically) expected. J. Geophys. Res., 121 (6),
#' 4575-4590, \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JB012818}{doi: 10.1002/2016JB012818}
model_rate.val <- function(method, theta, inj = NULL, shutin = NULL, t.postinj = NULL) {
  if(method != 'co-injection' & method != 'post-injection' & method != 'full sequence')
    stop('method not found. Use co-injection, post-injection, or full sequence')
  if (method == 'co-injection') {
    if(is.null(inj)) stop('injection profile missing')
    if(is.null(theta$a_fb) | is.null(theta$b) | is.null(theta$mc))
      stop('parameters a_fb, b and/or mc missing in theta')
    rate <- 10 ^ theta$a_fb * 10 ^ (-theta$b * theta$mc) * inj$dV
    t <- inj$t
  }
  if (method == 'post-injection') {
    if(is.null(shutin)) stop('shutin data missing')
    if(is.null(t.postinj)) stop('t.postinj time vector missing')
    if(is.null(theta$tau)) stop('parameter tau missing in theta')
    rate <- shutin$rate * exp(-(t.postinj - shutin$t) / theta$tau)
    t <- t.postinj
  }
  if (method == 'full sequence') {
    if(is.null(inj)) stop('injection profile required')
    if(is.null(theta$a_fb) | is.null(theta$b) | is.null(theta$mc) | is.null(theta$tau))
      stop('parameters a_fb, tau, b and mc required in theta')
    if(!is.null(shutin)) print('warning: shutin data overwritten by model')
    if(is.null(t.postinj)) stop('t.postinj time vector missing')
    rate.stimul <-  10 ^ theta$a_fb * 10 ^ (-theta$b * theta$mc) * inj$dV
    rate.relax <- rate.stimul[length(rate.stimul)] * exp(-(t.postinj - t.postinj[1]) / theta$tau)
    rate <- c(rate.stimul, rate.relax)
    t <- c(inj$t, t.postinj)
  }
  return(data.frame(t = t, rate = rate))
}

#' Poisson Log-Likelihood Function
#'
#' Computes the log-likelihood of the Poisson distribution
#'
#' @param obs a vector of observed rates
#' @param pred a vector of predicted rates
#' @return The Poisson log-likelihood estimate
poi.loglik <- function(obs = obs, pred = pred) {
  sum(obs * log(pred) - pred - log(factorial(obs)), na.rm = T)
}

#' Induced Seismicity Model Log-Likelihood Function
#'
#' Estimates the log-likelihood of the induced seismicity statistical
#' model `ratemodel.val()` using the Poisson formulation `poi.loglik()`.
#'
#' The function automatically determines the method ("`co-injection`",
#' "`post-injection`", or "`full sequence`") from `data` and `t.shutin`.
#' `t.shutin` can be set to `Inf` if unknown in the future.
#'
#' @param data The observed seismicity rate data frame
#' * `t` the time (in decimal days)
#' * `rate` the rate, or number of events at time `t`
#' @param theta.mle The list of model parameters
#' * `a_fb` the underground feedback activation (`NULL` for "`post-injection`")
#' * `tau` the mean relaxation time (`NULL` for "`co-injection`")
#' @param theta.GR The list of Gutenberg-Richter law parameters
#' * `mc` the completeness magnitude
#' * `b` the slope of the Gutenberg-Richter law
#' @param t.shutin The shut-in time (in decimal days)
#' @param inj the binned injection profile data frame with parameters (not required for "`post-injection`")
#' * `t` the occurrence time (in decimal days)
#' * `dV` the injected volume (in cubic metres) at `t`
#' @return the model's log-likelihood estimate
ratemodel.loglik <- function(data, theta.mle, theta.GR, t.shutin, inj = NULL) {
  if (is.null(t.shutin)) stop('t.shutin missing. Define t.shutin = Inf if unknown in the future')
  tmin <- min(data$t)
  tmax <- max(data$t)
  if (tmax <= t.shutin) method <- 'co-injection'
  if (tmin < t.shutin & tmax > t.shutin) method <- 'full sequence'
  if (tmin >= t.shutin) method <- 'post-injection'

  theta <- list(a_fb = theta.mle$a_fb, tau = theta.mle$tau, b = theta.GR$b, mc = theta.GR$mc)

  if (method == 'co-injection') rate.pred <- ratemodel.val('co-injection', theta, inj = inj)
  if (method == 'post-injection') {
    indpost <- data$t >= t.shutin
    shutin = list(rate = data$rate[indpost][1], t = data$t[indpost][1])
    rate.pred <- ratemodel.val('post-injection', theta, shutin = shutin, t.postinj = data$t[indpost])
  }
  if (method == 'full sequence') {
    indpost <- data$t >= t.shutin
    rate.pred <- ratemodel.val('full sequence', theta, inj = inj, t.postinj = data$t[indpost])
  }
  loglik <- poi.loglik(obs = data$rate, pred = rate.pred$rate)
  return(loglik)
}

#' Induced Seismicity model MLE
#'
#' Maximum Likelihood Estimation (MLE) of the induced seismicity model
#' parameters `a_fb` and `tau`.
#'
#' The maximization is done on the `ratemodel.loglik()` result.
#'
#' @param data The observed seismicity rate data frame
#' * `t` the time (in decimal days)
#' * `rate` the rate, or number of events at time `t`
#' @param theta.GR The list of Gutenberg-Richter law parameters
#' * `mc` the completeness magnitude
#' * `b` the slope of the Gutenberg-Richter law
#' @param t.shutin The shut-in time (in decimal days)
#' @param inj the binned injection profile data frame with parameters (not required for "`post-injection`")
#' * `t` the occurrence time (in decimal days)
#' * `dV` the injected volume (in cubic metres) at `t`
#' @return The list of estimated parameters
#' * `a_fb` the underground feedback activation (`NULL` for "`post-injection`")
#' * `tau` the mean relaxation time (`NULL` for "`co-injection`")
ratemodel.mle <- function(data, theta.GR, t.shutin, inj = NULL) {
  if (is.null(t.shutin)) stop('t.shutin missing. Define t.shutin = Inf if unknown in the future')
  tmin <- min(data$t)
  tmax <- max(data$t)
  if (tmax <= t.shutin) method <- 'co-injection'
  if (tmin < t.shutin & tmax > t.shutin) method <- 'full sequence'
  if (tmin >= t.shutin) method <- 'post-injection'

  if (method == 'co-injection') {
    a_fb.i <- seq(-8, 2, .1)

    loglik <- sapply(1:length(a_fb.i), function(i) ratemodel.loglik(data, list(a_fb = a_fb.i[i]),
                                                                    theta.GR, t.shutin, inj = inj))
    indmax <- which(loglik == max(loglik, na.rm = T))
    theta.mle <- list(a_fb = a_fb.i[indmax])
  }
  if (method == 'post-injection') {
    tau.i <- c(seq(.1,.9,.1), seq(1, 30))

    loglik <- sapply(1:length(tau.i), function(i) ratemodel.loglik(data, list(tau = tau.i[i]),
                                                                   theta.GR, t.shutin))
    indmax <- which(loglik == max(loglik, na.rm = T))
    theta.mle <- list(tau = tau.i[indmax])
  }
  if (method == 'full sequence') {
    a_fb.i <- seq(-8, 2, .1)
    tau.j <- c(seq(.1,.9,.1), seq(1, 30))

    loglik <- sapply(1:length(a_fb.i),
                     function(i) sapply(1:length(tau.j),
                                        function(j) ratemodel.loglik(data, list(a_fb = a_fb.i[i], tau = tau.j[j]),
                                                                     theta.GR, t.shutin, inj = inj)))

    indmax <- which(loglik == max(loglik, na.rm = T), arr.ind = T)
    theta.mle <- list(a_fb = a_fb.i[indmax[2]], tau = tau.j[indmax[1]])
  }
  if(theta.mle$a_fb == -8 | theta.mle$a_fb == 2) print('WARNING: parameter space not fully searched for a_fb')
  if(theta.mle$tau == .1 | theta.mle$tau == 30) print('WARNING: parameter space not fully searched for tau')

  return(theta.mle)
}

