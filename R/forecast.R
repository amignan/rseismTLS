#' Negative Log Likelihood Function (Point Data)
#'
#' Estimate the negative log likelihood of the statistical model of Mignan et al. (2017) for
#' a specific set of input parameters, for three possible time windows (full sequence, injection phase or
#' post-injection phase).
#'
#' See Eq. A2 of Broccardo et al. (2017) for `'full sequence'` and Eq. A3 for `'injection'`.
#' Note that the function is here structured for the parameters to be optimized by `optim()``, which is done
#' automatically in `model_par.mle_point()`. The `par` vector depends on the selected `window`, which is dealt with
#' in `model_par.mle_point()` (see Details section there).
#'
#' @param data a list containing all the necessary point data:
#' * `seism` an earthquake catalogue data frame of parameters `t` (time in days) and `m` (magnitude)
#' * `inj` matching injection profile data frame of parameters `t` (time in days), `dV` (flow rate in cubic
#' metre/day) and `V` (cumulative injected volume in cubic metres)
#' * `m0` the minimum magnitude cutoff of the earthquake catalogue
#' * `ts` shut-in time (in days) (not required for "`injection`")
#' * `Tmax` upper range of the time window (in days) (not required for "`injection`")
#' * `lambda0` seismicity rate (per day) at shut-in (only required for "`post-injection`")
#' @param par a vector of the input parameters
#' @param window the window to be used: "`injection`", "`post-injection`", or "`full sequence`"
#' @return the negative log likelihood estimate for the specified `par` values
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{model_par.mle_point}, \code{model_par.mle_hist}, \code{negloglik_hist.val}
negloglik_point.val <- function(data, par, window = 'full sequence'){
  if(window == 'full sequence'){
    theta <- list(a_fb = par[1], tau = par[2], b = par[3])

    N <- nrow(data$seism)
    seism.inj <- subset(data$seism, data$seism$t <= data$ts)
    seism.post <- subset(data$seism, data$seism$t > data$ts)
    Ninj <- nrow(seism.inj)
    Npost <- nrow(seism.post)
    require(signal)       # interp1()
    dV <- interp1(data$inj$t, data$inj$dV, seism.inj$t)
    dV.ts <- tail(dV, 1)
    V.ts <- tail(data$inj$V, 1)

    C1 <- N * (theta$a_fb - theta$b * data$m0) / log10(exp(1))
    log_dV <- log(dV); log_dV[!is.finite(log_dV)] <- NA  # case of dV = 0 during injection
    C2 <- sum(log_dV, na.rm = T)
    C3 <- Npost * log(dV.ts)
    C4 <- -1 / theta$tau * (sum(seism.post$t - data$ts))
    C5 <- -10^(theta$a_fb - theta$b * data$m0) * (V.ts + dV.ts * theta$tau * (1-exp(-(data$Tmax - data$ts)/theta$tau)))
    C6 <- N * log(theta$b)
    C7 <- N * log(log(10))
    C8 <- -theta$b * log(10) * sum(data$seism$m)
    C9 <- N  * theta$b * log(10) * (data$m0 - mbin/2)   # Mmax -> +Inf
    nLL <- -(C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9)
  }
  if(window == 'injection'){
    if(max(data$seism$t) > max(data$inj$t))
      stop('Seismicity temporal range must be within the injection profile range')
    theta <- list(a_fb = par[1], b = par[2])

    N <- nrow(data$seism)
    require(signal)       # interp1()
    dV <- interp1(data$inj$t, data$inj$dV, data$seism$t)
    V <- tail(data$inj$V, 1)

    C1 <- N * (theta$a_fb - theta$b * data$m0) / log10(exp(1))
    log_dV <- log(dV); log_dV[!is.finite(log_dV)] <- NA  # case of dV = 0 during injection
    C2 <- sum(log_dV, na.rm = T)
    C5b <- -10^(theta$a_fb - theta$b * data$m0) * V
    C6 <- N * log(theta$b)
    C7 <- N * log(log(10))
    C8 <- -theta$b * log(10) * sum(data$seism$m)
    C9 <- N  * theta$b * log(10) * (data$m0 - mbin/2)   # Mmax -> +Inf
    nLL <- -(C1 + C2 + C5b + C6 + C7 + C8 + C9)
  }
  if(window == 'post-injection'){
    theta <- list(tau = par[1], b = par[2])

    N <- nrow(data$seism)

    C3b <- N * log(data$lambda0)
    C4 <- -1 / theta$tau * (sum(data$seism$t - data$ts))
    C5b <- -data$lambda0 * theta$tau * (1 - exp(-(data$Tmax - data$ts)/theta$tau))
    C6 <- N * log(theta$b)
    C7 <- N * log(log(10))
    C8 <- -theta$b * log(10) * sum(data$seism$m)
    C9 <- N * theta$b * log(10) * (data$m0 - mbin/2)   # Mmax -> +Inf
    nLL <- -(C3b + C4 + C5b + C6 + C7 + C8 + C9)
  }
  return(nLL)
}

#' Model Parameter Maximum Likelihood Estimation (Point Data)
#'
#' Provides the maximum likelihood estimates (MLE) of the parameters of the statistical model of Mignan et al. (2017)
#' using the negative log likelihood function `negloglik_point.val()` for one of the three possible time windows (full sequence, injection phase or
#' post-injection phase).
#'
#' Optimization done using the `optim()`` function of the stats package.
#'
#' @param data a list containing all the necessary point data:
#' * `seism` an earthquake catalogue data frame of parameters `t` (time in days) and `m` (magnitude)
#' * `inj` matching injection profile data frame of parameters `t` (time in days), `dV` (flow rate in cubic
#' metre/day) and `V` (cumulative injected volume in cubic metres)
#' * `m0` the minimum magnitude cutoff of the earthquake catalogue
#' * `ts` shut-in time (in days) (not required for "`injection`")
#' * `Tmax` upper range of the time window (in days) (not required for "`injection`")
#' * `lambda0` seismicity rate (per day) at shut-in (only required for "`post-injection`")
#' @param theta.init an optional list of the initial values for the parameters to be optimized over
#' * `a_fb` the underground feedback activation (in /cubic metre)
#' * `tau` the mean relaxation time (in days)
#' * `b` the slope of the Gutenberg-Richter law
#' @param window the window to be used: "`injection`", "`post-injection`", or "`full sequence`"
#' @return a list of the parameters' MLEs:
#' * `a_fb` the underground feedback activation (in /cubic metre)
#' * `tau` the mean relaxation time (in days)
#' * `b` the slope of the Gutenberg-Richter law
#' * `nLL` the negative log likelihood value for the optimized parameters
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{negloglik_point.val}, \code{negloglik_hist.val}, \code{model_par.mle_hist}
model_par.mle_point <- function(data, theta.init = list(a_fb = -1, tau = 1, b = 1), window = 'full sequence') {
  if(window == 'full sequence') {
    par <- numeric(3)
    par[1] <- theta.init$a_fb; par[2] <- theta.init$tau; par[3] <- theta.init$b
  }
  if(window == 'injection') {
    par <- numeric(2)
    par[1] <- theta.init$a_fb; par[2] <- theta.init$b
  }
  if(window == 'post-injection') {
    par <- numeric(2)
    par[1] <- theta.init$tau; par[2] <- theta.init$b
  }
  res <- optim(par = par, fn = rseismTLS::negloglik_point.val, data = data, window = window)
  if(window == 'full sequence') {a_fb <- res$par[1]; tau <- res$par[2]; b <- res$par[3]}
  if(window == 'injection') {a_fb <- res$par[1]; tau <- NA; b <- res$par[2]}
  if(window == 'post-injection') {a_fb <- NA; tau <- res$par[1]; b <- res$par[2]}

  return(list(a_fb = a_fb, tau = tau, b = b, nLL = res$value))
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

  tbin <- unique(diff(tint))[1]
  require(signal)       # interp1()
  inj.binned <- data.frame(t = seism.binned$mids,
                           dV = interp1(c(0, inj$t), c(0, inj$dV), seism.binned$mids) * tbin,
                           V = interp1(c(0, inj$t), c(0, inj$V), seism.binned$mids))

  return(list(seism.binned = data.frame(t = seism.binned$mids, rate = seism.binned$counts),
              inj.binned = subset(inj.binned, !is.na(dV))))
}

#' Underground feedback activation
#'
#' Estimate the underground feedback activation \out{<i>a<sub>fb</sub></i>}, which is the
#' \out{<i>a</i>}-value of the Gutenberg-Richter law (Gutenberg and Richter, 1944)
#' normalized by the injected volume (e.g., Dinske and Shapiro, 2013).
#'
#' This is equivalent to Shapiro's Seismogenic Index (e.g., Dinske and Shapiro, 2013)
#' \out{<i>SI = log<sub>10</sub>(N<sub>tot</sub> / V<sub>tot</sub>) + b * m<sub>0</sub></i>}
#' but here theoretically agnostic, following the notation of Mignan et al. (2017).
#'
#' @param Ntot total number of events above `m0` for `Vtot`
#' @param m0 minimum magnitude threshold
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
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
a_fb.val <- function(Ntot, b, m0, Vtot) {
  log10(Ntot / Vtot) + b * m0
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
#' The binned injection profile is defined in the `data.bin()` function. The seismicity rate
#' is modelled for the same time vector as `inj` for the injection and uses `t.postinj` for post-injection.
#'
#' @param window the window to be used: "`injection`", "`post-injection`", or "`full sequence`"
#' @param theta the list of model parameters:
#' * `a_fb` the underground feedback activation (`NULL` for "`post-injection`")
#' * `tau` the mean relaxation time (`NULL` for "`injection`")
#' * `b` the slope of the Gutenberg-Richter law
#' * `m0` the minimum magnitude threshold
#' @param inj the binned injection profile data frame with parameters
#' * `t` the occurrence time (in decimal days)
#' * `dV` the flow rate (in cubic metres per bin)
#' @param shutin the list of shut-in parameters
#' * `t` the shut-in time (in decimal days) (required for "`full sequence`" and "`post-injection`")
#' * `rate` the rate of induced seismicity at shut-in (required for "`post-injection`")
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
model_rate.val <- function(window, theta, inj = NULL, shutin = NULL, t.postinj = NULL) {
  if(window != 'injection' & window != 'post-injection' & window != 'full sequence')
    stop("window not found. Use 'injection', 'post-injection', or 'full sequence'")
  if (window == 'injection') {
    if(is.null(inj)) stop('injection profile missing')
    if(is.null(theta$a_fb) | is.null(theta$b) | is.null(theta$m0))
      stop('parameters a_fb, b and/or m0 missing in theta')
    rate <- 10 ^ theta$a_fb * 10 ^ (-theta$b * theta$m0) * inj$dV
    t <- inj$t
  }
  if (window == 'post-injection') {
    if(is.null(shutin)) stop('shutin data missing')
    if(is.null(t.postinj)) stop('t.postinj time vector missing')
    if(is.null(theta$tau)) stop('parameter tau missing in theta')
    rate <- shutin$rate * exp(-(t.postinj - shutin$t) / theta$tau)
    t <- t.postinj
  }
  if (window == 'full sequence') {
    if(is.null(inj)) stop('injection profile required')
    if(is.null(theta$a_fb) | is.null(theta$b) | is.null(theta$m0) | is.null(theta$tau))
      stop('parameters a_fb, tau, b and m0 required in theta')
    if(is.null(shutin)) stop('shutin data missing')
    if(is.null(t.postinj)) stop('t.postinj time vector missing')
    rate.stimul <-  10 ^ theta$a_fb * 10 ^ (-theta$b * theta$m0) * inj$dV
    rate.relax <- rate.stimul[length(rate.stimul)] * exp(-(t.postinj - shutin$t) / theta$tau)
    rate <- c(rate.stimul, rate.relax)
    t <- c(inj$t, t.postinj)
  }
  return(data.frame(t = t, rate = rate))
}

#' Negative Log Likelihood Function (Histogram Data)
#'
#' Estimate the negative log likelihood of the Non-Homogeneous Poisson Process (NHPP) derived from the statistical
#' model of Mignan et al. (2017) for histogram data and a specific set of input parameters, for three possible time
#' windows (full sequence, injection phase or post-injection phase).
#'
#' Note that the function is here structured for the parameters to be optimized by `optim()``, which is done
#' automatically in `model_par.mle_hist()`. The `par` vector depends on the selected `window`, which is dealt with
#' in `model_par.mle_hist()` (see Details section there).
#'
#' @param data a list containing all the necessary histogram data:
#' * `seism.binned` an earthquake-rate histogram data frame of parameters `t` (time in days) and `rate` (per time bin)
#' * `inj.binned` matching binned injection profile data frame of parameters `t` (time in days), `dV` (flow rate in cubic
#' metre/time bin) and `V` (cumulative injected volume in cubic metres)
#' * `b` the b-value of the earthquake catalogue (not required for "`post-injection`")
#' * `m0` the minimum magnitude cutoff of the earthquake catalogue (not required for "`post-injection`")
#' * `ts` shut-in time (in days) (not required for "`injection`")
#' * `lambda0` seismicity rate (per day) at shut-in (only required for "`post-injection`")
#' @param par a vector of the input parameters
#' @param window the window to be used: "`injection`", "`post-injection`", or "`full sequence`"
#' @return the negative log likelihood estimate for the specified `par` values
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{model_par.mle_hist}, \code{model_par.mle_point}, \code{negloglik_point.val}, \code{data.bin}
negloglik_hist.val <- function(data, par, window = 'full sequence'){
  obs <- data$seism.binned$rate

  if(window == 'full sequence'){
    theta <- list(a_fb = par[1], tau = par[2])

    ind.post <- which(data$seism.binned$t > tail(data$inj.binned$t, 1))
    res <- rseismTLS::model_rate.val('full sequence',
                                     list(a_fb = theta$a_fb, tau = theta$tau, b = data$b, m0 = data$m0),
                                     inj = data$inj.binned, shutin = list(t = data$ts),
                                     t.postinj = data$seism.binned$t[ind.post])
  }
  if(window == 'injection'){
    theta <- list(a_fb = par)

    res <- rseismTLS::model_rate.val('injection',
                                     list(a_fb = theta$a_fb, b = data$b, m0 = data$m0), inj = data$inj.binned)
  }
  if(window == 'post-injection'){
    theta <- list(tau = par)

    ind.post <- data$seism.binned$t > data$ts
    res <- rseismTLS::model_rate.val('post-injection',
                                     list(tau = theta$tau), shutin = list(t = data$ts, rate = data$lambda0),
                                     t.postinj = data$seism.binned$t[ind.post])

  }
  pred <- res$rate
  nLL <- -sum(obs * log(pred) - pred - log(factorial(obs)), na.rm = T)
  return(nLL)
}

#' Model Parameter Maximum Likelihood Estimation (Histogram Data)
#'
#' Provides the maximum likelihood estimates (MLE) of the parameters of the Non-Homogeneous Poisson Process (NHPP) derived from
#' the statistical model of Mignan et al. (2017) using the negative log likelihood function `negloglik_hist.val()` for
#' one of the three possible time windows (full sequence, injection phase or post-injection phase).
#'
#' Optimization done using the `optim()`` function of the stats package. In contrast to point-data MLE parameter
#' optimization, the b-value is here an input instead of an output.
#'
#' @param data a list containing all the necessary histogram data:
#' * `seism.binned` an earthquake-rate histogram data frame of parameters `t` (time in days) and `rate` (per time bin)
#' * `inj.binned` matching binned injection profile data frame of parameters `t` (time in days), `dV` (flow rate in cubic
#' metre/time bin) and `V` (cumulative injected volume in cubic metres)
#' * `b` the b-value of the earthquake catalogue (not required for "`post-injection`")
#' * `m0` the minimum magnitude cutoff of the earthquake catalogue (not required for "`post-injection`")
#' * `ts` shut-in time (in days) (not required for "`injection`")
#' * `lambda0` seismicity rate (per day) at shut-in (only required for "`post-injection`")
#' @param theta.init an optional list of the initial values for the parameters to be optimized over
#' * `a_fb` the underground feedback activation (in /cubic metre)
#' * `tau` the mean relaxation time (in days)
#' @param window the window to be used: "`injection`", "`post-injection`", or "`full sequence`"
#' @return a list of the parameters' MLEs:
#' * `a_fb` the underground feedback activation (in /cubic metre)
#' * `tau` the mean relaxation time (in days)
#' * `nLL` the negative log likelihood value for the optimized parameters
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{negloglik_hist.val}, \code{negloglik_point.val}, \code{model_par.mle_point}, \code{data.bin}
model_par.mle_hist <- function(data, theta.init = list(a_fb = -1, tau = 1), window = 'full sequence') {
  if(window == 'full sequence') {
    par <- numeric(2)
    par[1] <- theta.init$a_fb; par[2] <- theta.init$tau
    res <- optim(par = par, fn = rseismTLS::negloglik_hist.val, data = data, window = window)
  }
  if(window == 'injection') {
    par <- theta.init$a_fb
    res <- optim(par = par, fn = rseismTLS::negloglik_hist.val, data = data, window = window,
                 method = 'Brent', lower = -10, upper = 5)   #wide a_fb range [-10, 5]
  }
  if(window == 'post-injection') {
    par <- theta.init$tau
    res <- optim(par = par, fn = rseismTLS::negloglik_hist.val, data = data, window = window,
                 method = 'Brent', lower = 0, upper = 50)    #wide tau range [0, 50]
  }
  if(window == 'full sequence') {a_fb <- res$par[1]; tau <- res$par[2]}
  if(window == 'injection') {a_fb <- res$par[1]; tau <- NA}
  if(window == 'post-injection') {a_fb <- NA; tau <- res$par[1]}

  return(list(a_fb = a_fb, tau = tau, nLL = res$value))
}
