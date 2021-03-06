% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast.R
\name{model_joint_prior.distr}
\alias{model_joint_prior.distr}
\title{Joint prior distribution estimation}
\usage{
model_joint_prior.distr(prior)
}
\arguments{
\item{prior}{a list of the prior distributions computed with \code{model_prior.distr()}
\itemize{
\item \code{bi} the vector of \code{b} increments
\item \code{b.prior} the prior Beta distribution of \code{b}
\item \code{ai} the vector of \code{a_fb} increments
\item \code{a.prior} the prior Beta distribution of \code{a_fb}
\item \code{taui} the vector of \code{tau} increments
\item \code{tau.prior} the prior Gamma distribution of \code{tau}
}}
}
\value{
the list of the joint prior distributions for each of the 3 model parameter combinations:
\itemize{
\item \code{b_a.prior} the array of the (\code{b}, \code{a_fb}) joint distribution
\item \code{b_tau.prior} the array of the (\code{b}, \code{tau}) joint distribution
\item \code{a_tau.prior} the array of the (\code{a_fb}, \code{tau}) joint distribution
\item \code{joint.prior_norm} the 3-dimensional array of the normalized joint prior distribution
}
}
\description{
Estimates the joint prior distributions related to the 3 parameters of the induced seismicity model of
Mignan et al. (2017), \code{model_rate.val()}.
}
\details{
Read Broccardo et al. (2017) for details.
}
\references{
Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian
Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367,
\href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL075251}{doi: 10.1002/2017GL075251}

Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form
traffic light system for actuarial decision-making during deep fluid injections. Sci. Rep., 7, 13607,
\href{https://www.nature.com/articles/s41598-017-13585-9}{doi: 10.1038/s41598-017-13585-9}
}
\seealso{
\code{model_rate.val}, \code{model_prior.distr}, \code{par_Mignan_etal_SciRep2017.dat}
}
