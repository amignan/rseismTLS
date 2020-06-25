#' Negative Log Likelihood Function (Point Data)
#'
#' Estimates the negative log likelihood of the statistical model of Mignan et al. (2017) for
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
#' @param par a vector of the input parameters `a_fb`, `tau`, and `b` (in this order)
#' @param window the window to be used: "`injection`", "`post-injection`", or "`full sequence`"
#' @return the negative log likelihood estimate for the specified `par` values
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{model_par.mle_point}, \code{model_par.mle_hist}, \code{negloglik_hist.val}, \code{loglik_point.array}
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
#' Bucketizes the two main data sets (injection `inj` and seismicity `seism`) in time intervals `tint`.
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
#' Estimates the underground feedback activation \out{<i>a<sub>fb</sub></i>}, which is the
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
#' Estimates the negative log likelihood of the Non-Homogeneous Poisson Process (NHPP) derived from the statistical
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

#' Time transformation
#'
#' Transforms earthquake occurrence times following the Ogata (1988) method based on the integration of the
#' induced seismicity model of `model_rate.val()`, to then be used as input in `stat.KS_uniform()`.
#'
#' See examples of induced seismicity applications in Mignan et al. (2017, 2019) and Broccardo et al. (2017).
#'
#' @param data a list containing all the necessary point data:
#' * `seism` an earthquake catalogue data frame of parameters `t` (time in days) and `m` (magnitude)
#' * `inj` matching injection profile data frame of parameters `t` (time in days), `dV` (flow rate in cubic
#' metre/day) and `V` (cumulative injected volume in cubic metres)
#' * `m0` the minimum magnitude cutoff of the earthquake catalogue
#' * `ts` shut-in time (in days)
#' @param theta the list of fitted model parameters:
#' * `a_fb` the underground feedback activation
#' * `tau` the mean relaxation time
#' * `b` the slope of the Gutenberg-Richter law
#' @return the transformed time vector
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2019), Autonomous Decision-Making Against Induced
#' Seismicity in Deep Fluid Injections. In: Ferrari A., Laloui L. (eds), Energy Geotechnics, SEG 2018, Springer
#' Series in Geomechanics and Geoengineering, 369-376,
#' \href{https://link.springer.com/chapter/10.1007/978-3-319-99670-7_46}{doi: 10.1007/978-3-319-99670-7_46}
#' @references Ogata Y. (1988), Statistical Models for Earthquake Occurrences and Residual Analysis for Point
#' Processes. J. Am. Stat. Assoc.,
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.1988.10478560#.XfJUpZNKjOQ}{83 (401), 9-27}
#' @seealso \code{stat.KS_uniform}, \code{model_rate.val}
t.transform <- function(data, theta) {
  require(caTools)	 # trapz()
  require(signal)	 # interp1()

  N <- nrow(data$seism)
  t.transf <- numeric(N)

  indinj <- which(data$seism$t <= data$ts)
  inj_atEvent <- data.frame(t = data$seism$t[indinj], dV = interp1(data$inj$t, data$inj$dV, data$seism$t[indinj]))
  rate.pred_atEvent <- rseismTLS::model_rate.val('full sequence',
                                                 list(a_fb = theta$a_fb, tau = theta$tau, b = theta$b, m0 = data$m0),
                                                 inj = inj_atEvent, shutin = list(t = data$ts),
                                                 t.postinj = data$seism$t[-indinj])

  rate.pred_atEvent <- rbind(data.frame(t = 0, rate = 0), rate.pred_atEvent)
  for(i in 2:(N+1)) t.transf[i-1] <- trapz(rate.pred_atEvent$t[1:i], rate.pred_atEvent$rate[1:i])
  return(t.transf)
}

#' Kolmogorov-Smirnov Goodness of Fit Test (uniform)
#'
#' Kolmogorov-Smirnov Goodness of Fit Test (K-S test) to quantify the level of performance of the model
#' `model_rate.val()` and to potentially make a residual analysis (e.g., Ogata, 1988). The test is made
#' assuming a uniform distribution (i.e. distribution of a stationary Poisson process of intensity 1) for the
#' time transformed in `t.transform()` for the induced seismicity model.
#'
#' See examples of induced seismicity applications in Mignan et al. (2017, 2019) and Broccardo et al. (2017).
#'
#' @param t.transf a vector of transformed times
#' @return the results of the K-S test as a data frame:
#' * `t` the sequence of N events, from 1 to N
#' * `lim95up` the upper 95% K-S confidence limit
#' * `lim95down` the lower 95% K-S confidence limit
#' * `lim99up` the upper 99% K-S confidence limit
#' * `lim99down` the lower 99% K-S confidence limit
#' * `boolean95` logical vector, true if events are with the 95% K-S confidence interval
#' * `boolean99` logical vector, true if events are with the 99% K-S confidence interval
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2019), Autonomous Decision-Making Against Induced
#' Seismicity in Deep Fluid Injections. In: Ferrari A., Laloui L. (eds), Energy Geotechnics, SEG 2018, Springer
#' Series in Geomechanics and Geoengineering, 369-376,
#' \href{https://link.springer.com/chapter/10.1007/978-3-319-99670-7_46}{doi: 10.1007/978-3-319-99670-7_46}
#' @references Ogata Y. (1988), Statistical Models for Earthquake Occurrences and Residual Analysis for Point
#' Processes. J. Am. Stat. Assoc.,
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.1988.10478560#.XfJUpZNKjOQ}{83 (401), 9-27}
#' @seealso \code{t.transform}, \code{model_rate.val}
stat.KS_uniform <- function(t.transf) {
  Yn <- diff(c(0, t.transf))
  Un <- 1 - exp(-Yn)
  Un_sort <- sort(Un)

  N <- length(t.transf)
  ti <- seq(N)
  Fx_u <- ti / N
  #95% for cum number of events
  bound_left_95 <- Fx_u - 2 / sqrt(N) * sqrt(Fx_u * (1 - Fx_u))
  bound_right_95 <- Fx_u + 2 / sqrt(N) * sqrt(Fx_u * (1 - Fx_u))
  d_KS_95 <- 1.358 / sqrt(N)
  #99% for cum number of events
  bound_left_99 <- Fx_u - 3 / sqrt(N) * sqrt(Fx_u * (1 - Fx_u))
  bound_right_99 <- Fx_u + 3 / sqrt(N) * sqrt(Fx_u * (1 - Fx_u))
  d_KS_99 <- 1.628 / sqrt(N)

  boolean95 <- seq(N) >= N * (t.transf / N - d_KS_95) & seq(N) <= N * (t.transf / N + d_KS_95)
  boolean99 <- seq(N) >= N * (t.transf / N - d_KS_99) & seq(N) <= N * (t.transf / N + d_KS_99)
  return(data.frame(t = ti,
                    lim95up = N * (Fx_u + d_KS_95),
                    lim95down = N * (Fx_u - d_KS_95),
                    lim99up = N * (Fx_u + d_KS_99),
                    lim99down = N * (Fx_u - d_KS_99),
                    boolean95 = boolean95,
                    boolean99 = boolean99))
}

#' Prior distribution estimation
#'
#' Estimates the prior distributions of the 3 parameters of the induced seismicity model of
#' Mignan et al. (2017), `model_rate.val()`, including the 3 joint prior distributions.
#'
#' The priors of `b` and `a_fb` are defined as Beta distributions while the prior of `tau` is defined as a
#' Gamma distribution. Read Broccardo et al. (2017) for details. The `par` file of Mignan et al. (2017)
#' (`par_Mignan_etal_SciRep2017.dat`) is provided as standard input. The proposed `bi`, `ai` and `taui` increments
#' are based on the ranges given in that file.
#'
#' @param par a data frame of model parameters fitted from past stimulations
#' * `b` the slope of the Gutenberg-Richter law
#' * `a_fb` the underground feedback activation (in m^-3)
#' * `tau` the mean relaxation time (in days)
#' @return the list of the prior distributions for each of the 3 model parameters:
#' * `bi` the vector of `b` increments
#' * `b.prior` the prior Beta distribution of `b`
#' * `ai` the vector of `a_fb` increments
#' * `a.prior` the prior Beta distribution of `a_fb`
#' * `taui` the vector of `tau` increments
#' * `tau.prior` the prior Gamma distribution of `tau`
#' * `b_a.prior` the array of the (`b`, `a_fb`) joint distribution
#' * `b_tau.prior` the array of the (`b`, `tau`) joint distribution
#' * `a_tau.prior` the array of the (`a_fb`, `tau`) joint distribution
#' * `joint.prior_norm` the 3-dimensional array of the normalized joint prior distribution
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{model_rate.val}, \code{model_joint_prior.distr}, \code{par_Mignan_etal_SciRep2017.dat}
model_prior.distr <- function(par, ai = seq(-5,1,.01), bi = seq(.5,2.,.01), taui = seq(.1,15,.05)) {
  require(MASS)    # fitdistr()

  n.a <- length(ai); n.b <- length(bi); n.tau <- length(taui)

  # b prior
  b0.norm <- (par$b - min(bi)) / (max(bi) - min(bi))
  b0.mu <- mean(b0.norm)
  b0.var <- var(b0.norm)
  b.alpha <- ((1 - b0.mu) / b0.var - 1 / b0.mu) * b0.mu ^ 2
  b.beta <- b.alpha * (1 / b0.mu - 1)
  b.prior <- dbeta( (bi - min(bi)) / (max(bi) - min(bi)), b.alpha, b.beta ) / (max(bi) - min(bi))

  # a prior
  a0.norm <- (par$a_fb - min(ai)) / (max(ai) - min(ai))
  a0.mu <- mean(a0.norm)
  a0.var <- var(a0.norm)
  a.alpha <- ((1 - a0.mu) / a0.var - 1 / a0.mu) * a0.mu ^ 2
  a.beta <- a.alpha * (1 / a0.mu - 1)
  a.prior <- dbeta( (ai - min(ai)) / (max(ai) - min(ai)), a.alpha, a.beta ) / (max(ai) - min(ai))

  # tau prior
  tau.par <- fitdistr(par$tau, 'Gamma')
  tau.shape <- tau.par$estimate[1]
  tau.rate <- tau.par$estimate[2]
  tau.prior <- dgamma(taui, tau.shape, tau.rate)

  abin <- unique(diff(ai))[1]
  bbin <- unique(diff(bi))[1]
  taubin <- unique(diff(taui))[1]

  # prior normalization for integration
  a.prior_norm <- a.prior / (sum(a.prior) * abin)
  b.prior_norm <- b.prior / (sum(b.prior) * bbin)
  tau.prior_norm <- tau.prior / (sum(tau.prior) * taubin)

  # bimarginal priors & joint prior distributions
  a.prior_norm2D <- matrix(rep(a.prior_norm, n.tau), nrow = n.a, ncol = n.tau)
  tau.prior_norm2D <- matrix(rep(tau.prior_norm, each = n.a), nrow = n.a, ncol = n.tau)
  a_tau.prior_norm2D <- a.prior_norm2D * tau.prior_norm2D
  joint.prior_norm <- array(NA, dim = c(n.a, n.tau, n.b))
  for(i in 1:n.b) joint.prior_norm[,, i] <- a_tau.prior_norm2D * b.prior_norm[i]
  b_a.prior <- sapply(1:length(ai), function(j) sapply(1:length(bi), function(i) sum(joint.prior_norm[j, , i]))) * taubin
  b_tau.prior <- sapply(1:length(taui), function(j) sapply(1:length(bi), function(i) sum(joint.prior_norm[, j, i]))) * abin
  a_tau.prior <- sapply(1:length(taui), function(j) sapply(1:length(ai), function(i) sum(joint.prior_norm[i, j, ]))) * bbin

  # partial prior (no tau)
  a.prior_norm2Dbis <- matrix(rep(a.prior_norm, n.b), nrow = n.a, ncol = n.b)
  b.prior_norm2D <- matrix(rep(b.prior_norm, each = n.a), nrow = n.a, ncol = n.b)
  joint.prior.partial_norm <- a.prior_norm2Dbis * b.prior_norm2D

  prior <- list(ai = ai, a.prior = a.prior, bi = bi, b.prior = b.prior, taui = taui, tau.prior = tau.prior,
                b_a.prior = b_a.prior, b_tau.prior = b_tau.prior, a_tau.prior = a_tau.prior,
                joint.prior_norm = joint.prior_norm, joint.prior.partial_norm = joint.prior.partial_norm)
  return(prior)
}

#' Log Likelihood Distribution (Point Data)
#'
#' Estimates the log likelihood distribution of the statistical model of Mignan et al. (2017) for
#' a specific set of input parameters. Function similar to `negloglik_point.val()` but for 3 vectors of input
#' parameters instead of 3 parameter values.
#'
#' See Eq. A2 of Broccardo et al. (2017) for details. Used for the hierarchical Bayesian modelling.
#'
#' @param data a list containing all the necessary point data:
#' * `seism` an earthquake catalogue data frame of parameters `t` (time in days) and `m` (magnitude)
#' * `inj` matching injection profile data frame of parameters `t` (time in days), `dV` (flow rate in cubic
#' metre/day) and `V` (cumulative injected volume in cubic metres)
#' * `m0` the minimum magnitude cutoff of the earthquake catalogue
#' * `ts` shut-in time (in days) (not required for "`injection`")
#' * `Tmax` upper range of the time window (in days) (not required for "`injection`")
#' @param par.space an array of the input parameters
#' * `bi` the vector of `b` increments
#' * `ai` the vector of `a_fb` increments
#' * `tau` the vector of `tau` increments (optional for `type = partial`)
#' @param type by default `complete` (full model), otherwise `partial` (injection phase only)
#' @return the log likelihood distribution for the specified `par.space` array
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{negloglik_point.val}, \code{model_posterior.distr}
loglik_point.array <- function(data, par.space, type = 'complete') {
  if(type == 'complete') {
    n.a <- length(par.space$ai); n.b <- length(par.space$bi); n.tau <- length(par.space$taui)
    LL <- array(NA, c(n.a, n.tau, n.b))

    N <- nrow(data$seism)
    seism.inj <- subset(data$seism, data$seism$t <= data$ts)
    seism.post <- subset(data$seism, data$seism$t > data$ts)
    Ninj <- nrow(seism.inj)
    Npost <- nrow(seism.post)
    dV <- interp1(data$inj$t, data$inj$dV, seism.inj$t)
    dV[dV <= 0] <- NA       #avoid -Inf for log(dV) - model applies to positive dV only
    dV.ts <- tail(dV, 1)
    V.ts <- tail(data$inj$V, 1)

    # components of length 1 or n.b
    C2 <- sum(log(dV), na.rm = T)
    C3 <- Npost * log(dV.ts)
    C6 <- N * log(par.space$bi)
    C7 <- N * log(log(10))
    C8 <- -par.space$bi * log(10) * sum(data$seism$m)
    C9 <- N  * par.space$bi * log(10) * (data$m0 - mbin/2)   # Mmax -> +Inf

    # other components as matrix [n.a, n.tau] for efficient computation
    Ai <- matrix(rep(par.space$ai, n.tau), nrow = n.a, ncol = n.tau)
    Taui <- matrix(rep(par.space$taui, each = n.a), nrow = n.a, ncol = n.tau)
    C4 <- -1 / Taui * (sum(seism.post$t - data$ts))
    for(i in 1:n.b){
      C1 <- N * (Ai - par.space$bi[i] * data$m0) / log10(exp(1))
      C5 <- -10^(Ai - par.space$bi[i] * data$m0) * (V.ts + dV.ts * Taui * (1-exp(-(data$Tmax - data$ts) / Taui)))

      LL[,,i] <- C1 + C2 + C3 + C4 + C5 + C6[i] + C7 + C8[i] + C9[i]
    }
  } else {
    n.a <- length(par.space$ai); n.b <- length(par.space$bi)
    LL <- array(NA, c(n.a, n.b))

    N <- nrow(data$seism)
    dV <- interp1(data$inj$t, data$inj$dV, data$seism$t)
    dV[dV <= 0] <- NA       #avoid -Inf for log(dV) - model applies to positive dV only
    V <- tail(data$inj$V, 1)

    Ai <- matrix(rep(par.space$ai, n.b), nrow = n.a, ncol = n.b)
    Bi <- matrix(rep(par.space$bi, each = n.a), nrow = n.a, ncol = n.b)

    C1 <- N * (Ai - Bi * data$m0) / log10(exp(1))
    C2 <- sum(log(dV), na.rm = T)
    C5b <- -10^(Ai - Bi * data$m0) * V
    C6 <- N * log(Bi)
    C7 <- N * log(log(10))
    C8 <- -Bi * log(10) * sum(data$seism$m)
    C9 <- N  * Bi * log(10) * (data$m0 - mbin/2)   # Mmax -> +Inf
    LL <- C1 + C2 + C5b + C6 + C7 + C8 + C9
  }
  return(LL)
}

#' Posterior distribution estimation
#'
#' Estimates the posterior distributions of the 3 parameters of the induced seismicity model of
#' Mignan et al. (2017), `model_rate.val()`, including joint posterior distributions.
#'
#' Read Broccardo et al. (2017) for details.
#'
#' @param prior the list of model parameter increments as defined in `model_prior.distr`:
#' * `bi` the vector of `b` increments
#' * `ai` the vector of `a_fb` increments
#' * `taui` the vector of `tau` increments (optional for `type = partial`)
#' * `joint.prior_norm` the 3-dimensional array of the normalized joint prior distribution  (for `type = complete`)
#' * `joint.prior.partial_norm` the 3-dimensional array of the normalized joint prior distribution (for `type = partial`)
#' @param LL the log likelihood distribution computed from `loglik_point.array`
#' @param type by default `complete` (full model), otherwise `partial` (injection phase only)
#' @return a list of the posterior distributions:
#' * `bi` the vector of `b` increments
#' * `ai` the vector of `a_fb` increments
#' * `taui` the vector of `tau` increments
#' * `a.post` the vector of the `a_fb` posterior distribution
#' * `b.post` the vector of the `b` posterior distribution
#' * `tau.post` the vector of the `tau` posterior distribution
#' * `b_a.post` the array of the (`b`, `a_fb`) joint distribution
#' * `b_tau.post` the array of the (`b`, `tau`) joint distribution
#' * `a_tau.post` the array of the (`a_fb`, `tau`) joint distribution
#' * `joint.post_norm` the 3-dimensional array of the normalized joint posterior distribution
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{model_prior.distr}, \code{model_joint_prior.distr}, \code{loglik_point.array}
model_posterior.distr <- function(prior, LL, type = 'complete'){
  if(type == 'complete') {
    abin <- unique(diff(prior$ai))[1]
    bbin <- unique(diff(prior$bi))[1]
    taubin <- unique(diff(prior$taui))[1]

    # joint posterior distribution
    # with log and minus max(LL) to avoid overflow
    joint.post_norm <- exp((LL - max(LL)) + log(prior$joint.prior_norm)) /
      (sum(exp((LL - max(LL)) + log(prior$joint.prior_norm))) * abin * taubin * bbin)

    # marginal posteriors
    a.post <- sapply(1:length(prior$ai), function(i) sum(joint.post_norm[i,,])) * taubin * bbin
    b.post <- sapply(1:length(prior$bi), function(i) sum(joint.post_norm[,, i])) * taubin * abin
    tau.post <- sapply(1:length(prior$taui), function(i) sum(joint.post_norm[, i,])) * abin * bbin

    # bimarginal posteriors
    b_a.post <- sapply(1:length(prior$ai), function(j) sapply(1:length(prior$bi), function(i) sum(joint.post_norm[j,, i]))) * taubin
    b_tau.post <- sapply(1:length(prior$taui), function(j) sapply(1:length(prior$bi), function(i) sum(joint.post_norm[, j, i]))) * abin
    a_tau.post <- sapply(1:length(prior$taui), function(j) sapply(1:length(prior$ai), function(i) sum(joint.post_norm[i, j, ]))) * bbin

    return(list(ai = prior$ai, bi = prior$bi, taui = prior$taui,
                a.post = a.post, b.post = b.post, tau.post = tau.post,
                b_a.post = b_a.post, b_tau.post = b_tau.post, a_tau.post = a_tau.post,
                joint.post_norm = joint.post_norm))
  } else {
    abin <- unique(diff(prior$ai))[1]
    bbin <- unique(diff(prior$bi))[1]

    # joint posterior distribution
    # with log and minus max(LL) to avoid overflow
    joint.post.partial_norm <- exp((LL - max(LL)) + log(prior$joint.prior.partial_norm)) /
      (sum(exp((LL - max(LL)) + log(prior$joint.prior.partial_norm))) * abin * bbin)

    # marginal posteriors
    a.post <- sapply(1:length(prior$ai), function(i) sum(prior$joint.post_norm_partial[i, ])) * bbin
    b.post <- sapply(1:length(prior$bi), function(i) sum(prior$joint.post_norm_partial[, i])) * abin

    return(list(ai = prior$ai, bi = prior$bi, a.post = a.post, b.post = b.post,
                joint.post.partial_norm = joint.post.partial_norm))
  }
}

#' Model Parameter Bayesian Estimation (Point Data)
#'
#' Provides the posterior mean, posterior mode (i.e. Maximum A Posteriori (MAP) estimate), and the
#' MLE (the later only if the log likelihood distribution `LL` is specified).
#'
#' Read Broccardo et al. (2017) for details.
#'
#' @param posterior the posterior distribution parameters estimated by `model_posterior.distr()`
#' @param type by default `complete` (full model), otherwise `partial` (injection phase only)
#' @param LL the log likelihood distribution estimated by `loglik_point.array()`
#' @return a data frame of the 3 parameter estimates based on 3 approaches:
#' * `a_fb.MAP` the MAP estimate of the underground feedback activation (in /cubic metre)
#' * `b.MAP` the MAP estimate of the slope of the Gutenberg-Richter law
#' * `tau.MAP` the MAP estimate of the mean relaxation time (in days) (if `type = complete`)
#' * `a_fb.mean` the posterior mean estimate of the underground feedback activation (in /cubic metre)
#' * `b.mean` the posterior mean estimate of the slope of the Gutenberg-Richter law
#' * `tau.mean` the posterior mean estimate of the mean relaxation time (in days) (if `type = complete`)
#' * `a_fb.MLE` the MLE of the underground feedback activation (in /cubic metre) - optional result
#' * `b.MLE` the MLE of the slope of the Gutenberg-Richter law - optional result
#' * `tau.MLE` the MLE of the mean relaxation time (in days) - optional result (if `type = complete`)
#' @references Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
#' Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
#' @seealso \code{negloglik_point.val}, \code{model_par.mle_point}, \code{model_par.mle_hist}
model_par.bayesian <- function(posterior, type = 'complete', LL = NULL){
  abin <- unique(diff(posterior$ai))[1]
  bbin <- unique(diff(posterior$bi))[1]
  if(type == 'complete') taubin <- unique(diff(posterior$taui))[1]

  #MAP
  a.MAP <- posterior$ai[posterior$a.post == max(posterior$a.post)]
  b.MAP <- posterior$bi[posterior$b.post == max(posterior$b.post)]
  if(type == 'complete') tau.MAP <- posterior$taui[posterior$tau.post == max(posterior$tau.post)]
  else tau.MAP <- NA

  #mean
  a.mean <- sum(posterior$ai * posterior$a.post) * abin
  b.mean <- sum(posterior$bi * posterior$b.post) * bbin
  if(type == 'complete') tau.mean <- sum(posterior$taui * posterior$tau.post) * taubin
  else tau.mean <- NA

  par <- data.frame(a_fb.MAP = a.MAP, b.MAP = b.MAP, tau.MAP = tau.MAP,
                       a_fb.mean = a.mean, b.mean = b.mean, tau.mean = tau.mean)

  if(!is.null(LL)){
    #MLE
    indmax <- which(LL == max(LL), arr.ind = T)
    a.MLE <- posterior$ai[indmax[1]]
    b.MLE <- posterior$bi[indmax[3]]
    if(type == 'complete') tau.MLE <- posterior$taui[indmax[2]] else tau.MLE < -NA

    par <- data.frame(par, data.frame(a_fb.MLE = a.MLE, b.MLE = b.MLE, tau.MLE = tau.MLE))
  }
  return(par)
}


