#' Data Binning
#'
#' Bucketize the two main raw data sets (injection `inj` and seismicity `seism`) in `dt` bins
#' and compute the seismicity rate `seism.rate` and injected volume `inj.binned` at times t
#'
#' Time `t` represents the past `dt` increment and is defined in the interval for which
#' data are available.
#'
#' @param seism an earthquake catalogue data frame of parameters:
#' * `t` the occurrence time (in decimal days)
#' * `m` the magnitude
#' @param inj matching injection profile data frame of parameters:
#' * `t` the occurrence time (in decimal days)
#' * `dV` the injected volume (in cubic metres)
#' * `V` the cumulative injected volume (in cubic metres)
#' @param dt time increment in same unit as `t`
#' @return A list of 2 data frames binned in time:
#' * `seism.rate` with parameters `t` and `rate` (number of events per `dt`)
#' * `inj.binned` with parameters `t`, `dV` and `V`
data.bin <- function(seism, inj, dt) {
  t.bin.seism <- seq(0, max(seism$t) + dt, dt)
  t.bin.inj <- seq(dt,  max(inj$t), dt)

  seism.bin <- hist(seism$t, breaks = t.bin.seism, plot = F)

  require(signal)       #interp1()
  V.bin <- interp1(inj$t, inj$V, t.bin.inj, method = "linear")
  dV.bin <- diff(c(0, V.bin))

  return(list(seism.rate = data.frame(t = seism.bin$mids + dt / 2, rate = seism.bin$counts),
              inj.binned = data.frame(t = t.bin.inj, V = V.bin, dV = dV.bin)))
}

#' Underground feedback activation
#'
#' Estimate the underground feedback activation \out{<i>a<sub>fb</sub></i>}, which is the
#' \out{<i>a</i>}-value of the Gutenberg-Richter law (Gutenberg and Richter, 1944)
#' normalized by the injected volume (e.g., Dinske and Shapiro, 2013).
#'
#' This is equivalent to Shapiro's Seismogenic Index (e.g., Dinske and Shapiro, 2013)
#' \out{<i>SI = log<sub>10</sub>(N<sub>tot</sub> / V<sub>tot</sub>) + b * m<sub>c</sub></i>}
#' but here theoretically agnostic, following the notation of Mignan et al. (2017)
#'
#' @param Ntot total number of events above `mc` for `Vtot`
#' @param b slope of the Gutenberg-Richter law
#' @param mc completeness magnitude
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
#' Description - COMING SOON
#'
#' Details - COMING SOON
#'
#' @param method COMING SOON
#' @param theta COMING SOON
#' @param inj COMING SOON
#' @param shutin COMING SOON
#' @param t.postinj COMING SOON
#' @return COMING SOON
ratemodel.val <- function(method, theta, inj = NULL, shutin = NULL, t.postinj = NULL) {
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


