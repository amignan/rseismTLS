% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast.R
\name{t.transform}
\alias{t.transform}
\title{Time transformation}
\usage{
\method{t}{transform}(data, theta)
}
\arguments{
\item{data}{a list containing all the necessary point data:
\itemize{
\item \code{seism} an earthquake catalogue data frame of parameters \code{t} (time in days) and \code{m} (magnitude)
\item \code{inj} matching injection profile data frame of parameters \code{t} (time in days), \code{dV} (flow rate in cubic
metre/day) and \code{V} (cumulative injected volume in cubic metres)
\item \code{m0} the minimum magnitude cutoff of the earthquake catalogue
\item \code{ts} shut-in time (in days)
}}

\item{theta}{the list of fitted model parameters:
\itemize{
\item \code{a_fb} the underground feedback activation
\item \code{tau} the mean relaxation time
\item \code{b} the slope of the Gutenberg-Richter law
}}
}
\value{
the transformed time vector
}
\description{
Transforms earthquake occurrence times following the Ogata (1988) method based on the integration of the
induced seismicity model of \code{model_rate.val()}, to then be used as input in \code{stat.KS_uniform()}.
}
\details{
See examples of induced seismicity applications in Mignan et al. (2017, 2019) and Broccardo et al. (2017).
}
\references{
Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
\href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}

Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
\href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}

Mignan A., Broccardo M., Wiemer S., Giardini D. (2019), Autonomous Decision-Making Against Induced
Seismicity in Deep Fluid Injections. In: Ferrari A., Laloui L. (eds), Energy Geotechnics, SEG 2018, Springer
Series in Geomechanics and Geoengineering, 369-376,
\href{https://link.springer.com/chapter/10.1007/978-3-319-99670-7_46}{doi: 10.1007/978-3-319-99670-7_46}

Ogata Y. (1988), Statistical Models for Earthquake Occurrences and Residual Analysis for Point
Processes. J. Am. Stat. Assoc.,
\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.1988.10478560#.XfJUpZNKjOQ}{83 (401), 9-27}
}
\seealso{
\code{stat.KS_uniform}, \code{model_rate.val}
}
