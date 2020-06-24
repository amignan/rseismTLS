#' 2006 Basel EGS injection profile
#'
#' A dataset containing the injected fluid volume during the 2006 stimulation
#' at the Basel, Switzerland, Enhanced Geothermal System. Digitzed (with 1-minute
#' resolution) from Fig. 5a of Haering M.O., Schanz U., Ladner F., Dyer B.C.
#' (2008), Characterisation of the Basel 1 enhanced geothermal system. Geothermics,
#' 37, 469-495, doi: 10.1016/j.geothermics.2008.06.002 - Post-injection bleed-off
#' not included.
#'
#' @format A data frame with 9310 rows and 3 variables:
#' \describe{
#'   \item{t}{time (in decimal days) since 2006-12-02 00:00:00.00}
#'   \item{dV}{volume (in cubic metres) injected in the minute preceding t}
#'   \item{V}{cumulative volume (in cubic metres) injected up to time t}
#' }
#' @source \url{https://www.sciencedirect.com/science/article/pii/S0375650508000382}
"Basel2006_inj"

#' Simulated version of the 2006 Basel EGS earthquake catalogue
#'
#' A dataset containing a stochastic version of the earthquakes induced during the 2006 stimulation
#' at the Basel, Switzerland, Enhanced Geothermal System. Simulated using the thinning method (Lewis & Shedler, 1979) with the
#' model parameters `a_fb`, `b`, and `tau` estimated from the earthquake catalogue of Kraft & Deichmann (2014)
#' and the injection profile `Basel2006_inj`. Includes events above the completeness magnitude and for 12 days
#' from the start of the stimulation.
#'
#' @format A data frame with 809 rows and 2 variables:
#' \describe{
#'   \item{t}{earthquake occurrence time (in decimal days) since injection start at 2006-12-02 18:00:00.00}
#'   \item{m}{earthquake magnitude}
#' }
#' @references Kraft T., Deichmann N. (2014), High-precision relocation and focal mechanism
#' of the injection-induced seismicity at the Basel EGS. Geothermics, 52, 59-73,
#' \href{https://www.sciencedirect.com/science/article/pii/S0375650514000686}{doi: 10.1016/j.geothermics.2014.05.014}
#' @references Lewis P.A.W., Shedler G.S. (1979), Simulation of nonhomogeneous Poisson processes by thinning.
#' Naval Res. Logistics Q., 26, 403–413, \href{https://apps.dtic.mil/dtic/tr/fulltext/u2/a059904.pdf}{doi:10.1002/nav.3800260304}
#' @seealso \code{Basel2006_inj}
"Basel2006_seism_simulated"

#' Induced seismicity model parameter list
#'
#' List of parameters `b`, `a_fb`, `tau` and `mc` estimaed for different stimulation `sites`, as given in
#' Mignan et al. (2017).
#' @format A data frame with the following columns:
#' \describe{
#'   \item{site}{stimulate site}
#'   \item{mc}{completeness magnitude}
#'   \item{b}{slope of the Gutenberg-Richter law}
#'   \item{a_fb}{underground feedback activation (in m^-3)}
#'   \item{tau}{mean relaxation time (in days)}
#' }
#' @source \url{https://www.nature.com/articles/s41598-017-13585-9}
#' @references Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
#' traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
#' \href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
"par_Mignan_etal_SciRep2017"
