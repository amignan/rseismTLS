#' rseismTLS: Induced Seismicity Traffic Light System
#'
#' Provides the tools to build a dynamic traffic light system (TLS)
#' to mitigate the risk of induced seismicity during underground stimulation
#' by fluid injection.
#'
#' @section data
#' "Basel2006_inj", "Basel2006_seism_simulated", "par_Mignan_etal_SciRep2017"
#' @section rseismTLS forecast functions:
#' negloglik_point.val, model_par.mle_point, a_fb.val,
#' data.bin, model_rate.val, negloglik_hist.val, model_par.mle_hist,
#' t.transform, stat.KS_uniform,
#' model_prior.distr, loglik_point.array, model_posterior.distr, model_par.bayesian,
#' forecast.seism
#' @section rseismNet risk functions: N/A
#' @section utils
#' rejection_sampling
#'
#' @docType package
#' @name rseismTLS

NULL
