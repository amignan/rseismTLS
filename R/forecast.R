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


