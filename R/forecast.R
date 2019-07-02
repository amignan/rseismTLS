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
