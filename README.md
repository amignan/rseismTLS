---
title: "rseismTLS"
author: "Arnaud Mignan"
date: "2020-06-30"
output:
  html_document:
    toc: true
    toc_depth: 2
    keep_md: true
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



<!-- badges: start -->
<!-- badges: end -->

The rise in the frequency of anthropogenic earthquakes due to deep fluid injections is posing serious economic, societal, and legal challenges to many geo-energy and waste-disposal projects. The rseismTLS package provides actuarial tools to analyse, forecast and mitigate induced seismicity during underground stimulation, based on a TLS (i.e. Traffic Light System) procedure. It reflects the need to quantify the dynamic nature of the industrial operations and underground feedback, and complements the standard TLS approach that is inherently heuristic and mostly based on expert elicitation (so-called clinical judgment).

What rseismTLS does:

- Models the seismicity induced by fluid injection in the underground, based on the linear relationship between flow rate and seismicity rate;
- Models the seismicity trailing effect (post-injection) based on an exponential relaxation function;
- Estimates the model parameters with two alternative approaches: frequentist and Bayesian;
- Forecasts future seismicity based on planned underground stimulation (COMING SOON);
- Estimates the magnitude threshold above which stimulation must be stopped to respect a risk-based safety criterion (COMING SOON);
- Computes the seismic risk for a point source based on the RISK-UE method and EMS98 scale (COMING LATER);
- Estimates the probability of reaching the magnitude threshold based on the seismicity forecast (COMING LATER).

What rseismTLS needs:

- An earthquake catalogue (one induced seismicity sequence so far included in `/data` for testing);
- The matching fluid injection profile (one injection profile so far included in `/data` for testing).

## Disclaimer

The rseismTLS package is provided "as is", without warranty of any kind. In no event shall the author or copyright holder be liable for any claim, damages or other liability (see MIT `LICENSE` file).

## Installation

You can install rseismTLS from github with:


```r
# install.packages("devtools")
devtools::install_github("amignan/rseismTLS")
```

## Input example

One dataset is so far provided (see `data/`), which is a simulated version of the 2006 Basel, Switzerland, EGS experiment (Häring et al., 2008; Kraft and Deichmann, 2014). The injection profile `Basel2006_inj` was digitized from Häring et al. (2008) and the catalogue `Basel2006_seism_simulated` derived from the dataset provided by Kraft and Deichmann (2014), here limited to the initial 12 days (c. 6 days of stimulation and 6 days of post-injection decay) and to magnitude information only (spatial component not included). To avoid any copyright infringement, we only provide a stochastic iteration of the Basel data (based on the thinning method (Lewis and Shedler, 1979) and the intensity function of Mignan et al. (2017) - see below). We will use this dataset to illustrate the various functionalities of rseismTLS. We will consistently use the following naming convention: `inj` for the injection profile and `seism` for the earthquake catalogue. User-created data sets must be in the following units: decimal days for time and cubic metres for injected volumes. Moreover, since the model so far does not include negative flow rates (i.e., bleed-off in post-injection phase), the injection flow rate must stop at the shut-in time and any negative flow rate during the injection phase be fixed to zero (this case is fortunately very rare). Following those rules will conform with the rseismTLS computational framework.


```r
## Load the simulated 2006 Basel EGS data ##
inj <- rseismTLS::Basel2006_inj
seism <- rseismTLS::Basel2006_seism_simulated
# time t [days], flow rate dV [m^3/day], cumulative volume V [m^3]
tail(inj, 5)       
#>          t       dV        V
#> 36 5.91488 4569.869  9652.23
#> 37 6.15537 4680.778 10777.91
#> 38 6.17162 2626.776 10820.59
#> 39 6.46357 2603.563 11580.71
#> 40 6.48125 2603.563 11626.74
#
# same time t reference, event magnitudes m
head(seism, 5)
#>           t        m
#> 1 0.9142506 2.030908
#> 2 0.9387515 1.005578
#> 3 1.2395748 1.018013
#> 4 1.4741568 1.210736
#> 5 1.5500578 1.191206
```


```r
# illustrate how decimal days are computed
initdate <- strptime("2006-12-02 00:00:00.00", "%Y-%m-%d %H:%M:%OS")
# start injection = 2 December, 6:00 pm
startdate <- strptime("2006-12-02 18:00:00.00", "%Y-%m-%d %H:%M:%OS")
t.start <- as.double(startdate)/86400-as.double(initdate)/86400
# end injection = 8 December, 11:33 am
enddate <- strptime("2006-12-08 11:33:00.00", "%Y-%m-%d %H:%M:%OS")
( t.shutin <- as.double(enddate)/86400-as.double(initdate)/86400 )
#> [1] 6.48125
# note that the time of the last row in inj (see above) must equal t.shutin

# plot data
par(mfrow = c(2,2), mar = c(4, 4, 1.5, 0.5))
plot(inj$t, inj$dV, type = 'l', xlim = c(0,12), main = 'Flow rate dV [m^3/day]')
abline(v = c(t.start, t.shutin), lty = 'dotted', col = 'red')

plot(inj$t, inj$V, type = 'l', xlim = c(0,12), main = 'cum. volume V [m^3]')
abline(v = c(t.start, t.shutin), lty = 'dotted', col = 'red')

plot(seism$t, seism$m, xlim = c(0,12), pch = 20, cex = seism$m, main = 'magnitude m')
for(i in 1:nrow(seism)) lines(c(seism$t[i], seism$t[i]), c(0, seism$m[i]))
abline(v = c(t.start, t.shutin), lty = 'dotted', col = 'red')

hist(seism$t, breaks = seq(0, 12, 1/24), xlim = c(0,12), main = 'hourly rate')
abline(v = c(t.start, t.shutin), lty = 'dotted', col = 'red')
```

<img src="figs/unnamed-chunk-3-1.png" width="100%" />

## Forecast functions

The following examples make use of all the forecast functions of rseismTLS, which are listed in `R/forecast.R`. They are based on the statistical model of Mignan et al. (2017):

$$\lambda (t, m \ge m_c; \theta) = \begin{cases} 10^{a_{fb}-b m_c} \dot V (t) & ; t \le t_{shut-in} \\ 
  10^{a_{fb}-b m_c} \dot V (t_{shut-in}) \exp (- \frac{t-t_{shut-in}}{\tau})  & ; t > t_{shut-in} \end{cases}$$

where $\lambda$ is the predicted seismicity rate above completeness magnitude $m_c$, $\dot V$ the flow rate, $a_{fb}$ the underground activation feedback (equivalent to the seismogenic index; e.g. Dinske and Shapiro, 2013) and $\tau$ the mean relaxation time (Mignan et al., 2017). The linear relationship between seismicity rate $\lambda$ and flow rate $\dot V$ is well recognised (Dinske and Shapiro, 2013; Mignan, 2016; van der Elst et al., 2016; Mignan et al., 2017; Broccardo et al., 2017). Use of the exponential decay during the post-injection phase was proposed by Mignan et al. (2017).

Only the frequentist approach developed by Mignan et al. (2017) is so far available. The bayesian approach developed by Broccardo et al. (2017) will be implemented at a later date. Physics-based models are so far not considered. 

### Frequentist approach (likelihood)

#### Data preprocessing
Note that we first need to use functions from the rseismNet package for data preprocessing, here to compute the completeness magnitude $m_c$ (see [rseismNet README](https://amignan.github.io/rseismNet/README.html) for details).


```r
## mandatory preprocessing ##
#devtools::install_github("amignan/rseismNet")   #in case rseismNet not yet installed
mbin <- .1
mc <- rseismNet::mc.val(seism$m, method = 'mode', mbin)      # see other methods in rseismNet
seism <- subset(seism, m > mc - mbin / 2)                    # complete earthquake catalogue
(b.Aki <- rseismNet::beta.mle(seism$m, mc, mbin) / log(10))  # Aki (1965) MLE method
#> [1] 1.606944
```

#### Maximum likelihood estimation method
The model parameters are estimated with `model_par.mle_point()` by minimizing the negative log-likelihood function `negloglik_point.val()`, as given in Broccardo et al. (2017:A2):


```r
data <- list(seism = seism, inj = inj, m0 = mc, ts = t.shutin, Tmax = 12)
( par.MLE <- rseismTLS::model_par.mle_point(data) )
#> Loading required package: signal
#> 
#> Attaching package: 'signal'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, poly
#> $a_fb
#> [1] 0.100735
#> 
#> $tau
#> [1] 1.151804
#> 
#> $b
#> [1] 1.606875
#> 
#> $nLL
#> [1] -2595.884
```

Note that the $b$ estimate is close to the MLE value obtained directly from the frequency-magnitude distribution (see `b.Aki`). One could also calculate $a_{fb}$ via the Shapiro Seismogenic Index (SI) method (e.g., Dinske and Shapiro, 2013) although that approach is less robust.


```r
( Vtot <- tail(inj$V, 1) )        # cubic meters
#> [1] 11626.74
( a_fb.SI <- rseismTLS::a_fb.val(nrow(seism), b.Aki, mc, Vtot) )
#> [1] 0.1996769
```

#### Model validation
The Kolmogorov-Smirnov Goodness of Fit Test (K-S test) can then be applied to quantify the level of performance of the model and to potentially make a residual analysis (e.g., Ogata, 1988). The test is made with `stat.KS_uniform()` assuming a uniform distribution (i.e. distribution of a stationary Poisson process of intensity 1) for the transformed time $\tilde{t}_i$ (Ogata, 1988; Broccardo et al., 2017; Mignan et al., 2017), which is defined as

$$\tilde{t}_i = \int^{t_i}_0 \lambda(t; \theta)dt$$
by the function `t.transform()`. Deviations from the 95% or 99% K-S confidence intervals can then be considered failures of the model.


```r
data <- list(seism = seism, inj = inj, m0 = mc, ts = t.shutin)
t.transf <- rseismTLS::t.transform(data, par.MLE)
#> Loading required package: caTools
res.KS <- rseismTLS::stat.KS_uniform(t.transf)

par(mfrow = c(2,1), mar = c(4, 4, 1.5, 0.5))
plot(t.transf, seism$m, pch = 20, cex = seism$m,
     xlab = 'Transformed time', ylab = 'Magnitude', main = 'Transformed time series')
for(i in 1:nrow(seism)) lines(c(t.transf[i], t.transf[i]), c(mc, seism$m[i]))

col.KS <- rep('darkgreen', length(t.transf))
col.KS[!res.KS$boolean95] <- 'orange'
col.KS[!res.KS$boolean99] <- 'red'

plot(t.transf, res.KS$t, type = 'b', cex = 0.5, pch = 20, col = col.KS, main = 'K-S test')
lines(res.KS$t, res.KS$lim95down, lty = 'dashed')
lines(res.KS$t, res.KS$lim95up, lty = 'dashed')
lines(res.KS$t, res.KS$lim99down, lty = 'dotted')
lines(res.KS$t, res.KS$lim99up, lty = 'dotted')
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

#### Model prediction visualisation
The predicted seismicity rate $\lambda$ per time bin `tbin` is calculated by `model_rate.val()`, which input is binned by `data.bin()`. This can be used to visualize the rate model against the binned data:


```r
# bucketize injection & seismicity data into tbin time bins
tbin <- 1/6                # [days]
tint <- seq(0, 12, tbin)    # [days]
data.binned <- rseismTLS::data.bin(seism, inj, tint)
seism.binned <- data.binned$seism.binned
inj.binned <- data.binned$inj.binned

ind.post <- which(seism.binned$t > tail(inj.binned$t, 1))
rate.pred <- rseismTLS::model_rate.val('full sequence', 
                                       list(a_fb = par.MLE$a_fb, tau = par.MLE$tau, b = par.MLE$b, m0 = mc), 
                                       inj = inj.binned, shutin = list(t = t.shutin), 
                                       t.postinj = seism.binned$t[ind.post])

par(mfrow = c(1, 1), mar = c(4, 4, 1.5, 0.5))
hist(seism$t, breaks = tint, main = 'full sequence', col = 'grey', border = 'white')
lines(rate.pred, col = 'darkred')
abline(v = c(t.start, t.shutin), col = 'red', lty = 'dotted')
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

We used the option `window = 'full sequence'`, which means that the rate sequence is continuous. Considering the `injection` and `post-injection` separately can lead to a discontinuity (play with different values of `tbin` to test the impact of binning). Those two options should only be used when data is limited to one of the two phases.


```r
rate_inj.pred <- rseismTLS::model_rate.val('injection', 
                                             list(a_fb = par.MLE$a_fb, b = par.MLE$b, m0 = mc),
                                             inj = inj.binned)

rate_postinj.pred <- rseismTLS::model_rate.val('post-injection',
                                               list(tau = par.MLE$tau),
                                               shutin = list(t = t.shutin, 
                                                             rate = seism.binned$rate[ind.post[1] - 1]),
                                               t.postinj = seism.binned$t[ind.post])

par(mar = c(4, 4, 1.5, 0.5))
hist(seism$t, breaks = tint, main = 'injection + post-injection', col = 'grey', border = 'white')
lines(rate_inj.pred$t, rate_inj.pred$rate, col = 'darkred')
lines(rate_postinj.pred$t, rate_postinj.pred$rate, col = 'darkblue')
abline(v = c(t.start, t.shutin), col = 'red', lty = 'dotted')
```

<img src="figs/unnamed-chunk-9-1.png" width="100%" />

In the case of independent time windows on each side of the shut-in time, the model parameters should be estimated for those specific periods. This is done as previously, but now adding the `window` type in `model_par.mle_point()`:


```r
data_inj <- list(seism = subset(seism, t <= t.shutin), inj = inj, m0 = mc, ts = t.shutin, Tmax = 12)
( par_inj.MLE <- rseismTLS::model_par.mle_point(data_inj, window = 'injection') )
#> $a_fb
#> [1] 0.1328267
#> 
#> $tau
#> [1] NA
#> 
#> $b
#> [1] 1.646938
#> 
#> $nLL
#> [1] -2166.312

data_post <- list(seism = subset(seism, t > t.shutin), m0 = mc, ts = t.shutin, Tmax = 12, 
                  lambda0 = seism.binned$rate[ind.post[1] - 1] / tbin)
( par_post.MLE <- rseismTLS::model_par.mle_point(data_post, list(a_fb = 0, tau = 1, b = 1), 
                                                 window = 'post-injection') )
#> $a_fb
#> [1] NA
#> 
#> $tau
#> [1] 1.138642
#> 
#> $b
#> [1] 1.473208
#> 
#> $nLL
#> [1] -430.3742
```

The parameter estimates here differ from the ones obtained previously due to the fact that there is no constraint to have a continuous function over the full sequence. For the Basel case, it does not lead to significant changes in $a_{fb}$ and $\tau$. However, we verify that the $b$-value decreases over time, as already known for this sequence (e.g., Mignan et al., 2017 - compare `par_inj.MLE$b` to `par_post.MLE$b`).

#### NB: when the earthquake count data is available instead of an earthquake catalogue
In the case in which only the seismicity rate `seism.binned` is available (e.g, catalogue not available but histogram data digitizable from a published figure), the model parameters can still be estimated. Then the Poisson log-likelihood function `negloglik_hist.val()` is used instead of `negloglik_point.val()`, assuming a Non-Homogeneous Poisson Process (NHPP). It is computed, and the model parameters optimized, in `model_par.mle_hist()`:


```r
# MLE values based on binned data
data <- list(seism.binned = seism.binned, inj.binned = inj.binned, b = b.Aki, m0 = mc, ts = t.shutin)
( par.MLE <- rseismTLS::model_par.mle_hist(data) )
#> $a_fb
#> [1] 0.09704507
#> 
#> $tau
#> [1] 1.15937
#> 
#> $nLL
#> [1] 149.1395

data_inj <- list(seism.binned = subset(seism.binned, t <= t.shutin), inj.binned = inj.binned, b = b.Aki, m0 = mc)
( par_inj.MLE <- rseismTLS::model_par.mle_hist(data_inj, window = 'injection') )
#> $a_fb
#> [1] 0.09193446
#> 
#> $tau
#> [1] NA
#> 
#> $nLL
#> [1] 97.43646

data_post <- list(seism.binned = subset(seism.binned, t > t.shutin), ts = t.shutin,
                  lambda0 = tail(seism.binned$rate[seism.binned$t <= t.shutin], 1))
( par_post.MLE <- rseismTLS::model_par.mle_hist(data_post, window = 'post-injection') )
#> $a_fb
#> [1] NA
#> 
#> $tau
#> [1] 1.141612
#> 
#> $nLL
#> [1] 51.51799
```


### Bayesian approach

The following Bayesian hierarchical framework was proposed by Broccardo et al. (2017). It is based on a nonhomogeneous Poisson process that follows the induced seismicity model of Mignan et al. (2017).

#### Prior distribution estimation

The prior distributions are selected as follows: Beta for $b$ and $a_{fb}$, and Gamma for $\tau$. The hyperparameters of the prior distributions are estimated in `model_prior.distr()` from the model parameters `par` fitted for various deep fluid stimulations by Mignan et al. (2017) (see dataset `par_Mignan_etal_SciRep2017.dat`). 


```r
par <- rseismTLS::par_Mignan_etal_SciRep2017
ndat <- nrow(par)
prior <- rseismTLS::model_prior.distr(par)
#> Loading required package: MASS
#> Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

par(mfrow = c(1,3))
plot(prior$ai, prior$a.prior, type = 'l')
points(par$a_fb, rep(0, ndat))
plot(prior$bi, prior$b.prior, type = 'l')
points(par$b_fb, rep(0, ndat))
plot(prior$taui, prior$tau.prior, type = 'l')
points(par$tau, rep(0, ndat))
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

In addition to the marginal prior distributions, `model_prior.distr()` also provides the 3 joint prior distributions:


```r
par(mfrow = c(1,3))
image(prior$bi, prior$ai, prior$b_a.prior)
image(prior$bi, prior$taui, prior$b_tau.prior)
image(prior$ai, prior$taui, prior$a_tau.prior)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

#### Posterior distribution estimation

The Bayesian framework enables the computation of posterior distribution for the model’s parameters, It first requies the computation of the log likelihood distribution by using `loglik_point.array()`. The posterior distributions are then obtained via `model_posterior.distr()`.


```r
data <- list(seism = seism, inj = inj, m0 = mc, ts = t.shutin, Tmax = 12)

# compute log likelihood distribution 
LL <- rseismTLS::loglik_point.array(data, prior)

# compute posterior distribution
posterior <- rseismTLS::model_posterior.distr(prior, LL)

par(mfrow = c(1,3))
plot(posterior$ai, posterior$a.post, type = 'l')
plot(posterior$bi, posterior$b.post, type = 'l')
plot(posterior$taui, posterior$tau.post, type = 'l')
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

Finally, the model parameters are estimated with `model_par.bayesian()`. Five estimates are given, the posterior mean, posterior mode (i.e. Maximum A Posteriori (MAP) estimate), lower and upper bounds of the posterior 90% credible interval, and the MLE (the later only if the log likelihood distribution `LL` is specified). 


```r
( res <- rseismTLS::model_par.bayesian(posterior, LL) )
#> $a_fb.MAP
#> [1] 0.07
#> 
#> $b.MAP
#> [1] 1.57
#> 
#> $tau.MAP
#> [1] 1.15
#> 
#> $a_fb.mean
#> [1] 0.06739245
#> 
#> $b.mean
#> [1] 1.571804
#> 
#> $tau.mean
#> [1] 1.160416
#> 
#> $a_fb.CI
#>          5%         95% 
#> -0.02311527  0.16074375 
#> 
#> $b.CI
#>       5%      95% 
#> 1.475901 1.669262 
#> 
#> $tau.CI
#>       5%      95% 
#> 1.025076 1.308376 
#> 
#> $a_fb.MLE
#> [1] 0.1
#> 
#> $b.MLE
#> [1] 1.61
#> 
#> $tau.MLE
#> [1] 1.15
```

We here plot the MAP estimates on the joint posterior distributions for comparison:


```r
par(mfrow = c(1,3))
image(posterior$bi, posterior$ai, posterior$b_a.post, xlim = c(1.4, 1.8), ylim = c(-.1, .3))
abline(v = res$b.MAP, h = res$a_fb.MAP)
image(posterior$bi, posterior$taui, posterior$b_tau.post, xlim = c(1.4, 1.8), ylim = c(.8, 1.4))
abline(v = res$b.MAP, h = res$tau.MAP)
image(posterior$ai, posterior$taui, posterior$a_tau.post, xlim = c(-.1, .3), ylim = c(.8, 1.4))
abline(v = res$a_fb.MAP, h = res$tau.MAP)
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />

#### Seismicity rate forecasting

Now that the functions to estimate the model parameters have been presented, we are ready to forecast future seismicity. The function `forecast.seism()` takes all the data and slices it in regular time windows of length `forecast.twin` (fixed below to $\Delta t = 1/6$ of a day) for pseudo-prospective forecasting. Note that we here use the default fast MAP method (i.e. ). The full posterior distribution can be used instead by selecting `method = 'bayesFull'`, which provides a better evaluation of the forecast distribution tail. Note however that this method becomes very slow in the case of a high-resolution parameter space (consider changing then the input parameter ranges `ai`, `bi` and `taui` in the function `model_prior.distr()`). Let us first plot the MAP parameters estimated from the observations up to time $t$.


```r
forecast.twin <- 1/6
forecast <- rseismTLS::forecast.seism(data, prior, forecast.twin)

plot(forecast$ti, forecast$a_fb, type = 'l', ylim = c(-2,2), col = 'red', ylab = 'MAP estimates')
lines(forecast$ti, forecast$a_fb.CI[, 1], lty = 'dotted', col = 'red')
lines(forecast$ti, forecast$a_fb.CI[, 2], lty = 'dotted', col = 'red')
lines(forecast$ti, forecast$b, col = 'blue')
lines(forecast$ti, forecast$b.CI[, 1], lty = 'dotted', col = 'blue')
lines(forecast$ti, forecast$b.CI[, 2], lty = 'dotted', col = 'blue')
lines(forecast$ti, forecast$tau, col = 'darkgreen')
lines(forecast$ti, forecast$tau.CI[, 1], lty = 'dotted', col = 'darkgreen')
lines(forecast$ti, forecast$tau.CI[, 2], lty = 'dotted', col = 'darkgreen')
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />

Note the high uncertainties when data is scarce. In that case, the parameter values are close to the prior estimates. The three parameters then converve to their site-specific values when more data comes in.
Now, let us plot the forecasted count of earthquakes in the period $(t, t+\Delta t]$:


```r
forecast.bins <- seq(forecast.twin, data$Tmax, forecast.twin)
res <- hist(seism$t, breaks = forecast.bins, main = 'Pseudo-prospective forecast (bayesMAP)', col = 'grey', border = 'white', ylim = c(0, 70))
polygon(x = c(forecast$ti, rev(forecast$ti)), y = c(forecast$N.CI[, 1], rev(forecast$N.CI[, 2])), col = rgb(1, 0, 0,0.5), border = NA)
boolean90 <- res$counts >= forecast$N.CI[, 1] & res$counts <= forecast$N.CI[, 2]
points(res$mids[boolean90], rep(-1, length(which(boolean90))), col = 'darkgreen', pch = 15, cex = .5)
points(res$mids[!boolean90], rep(-1, length(which(!boolean90))), col = 'darkred', pch = 15, cex = .5)
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" />

```r

( length(which(boolean90)) / length(res$mids) )
#> [1] 0.8309859
```

Can we reduce the number of bins falling outside the 90% credible interval? The fast `bayesMAP` method does not model the parameter tails so we might improve the forecast by modelling the entire posterior distribution. First we will reduce the resolution of the parameter space and then rerun `forecast.seism()`:


```r
prior2 <- rseismTLS::model_prior.distr(par, ai = seq(-5,1,.1), bi = seq(.5,2.,.1), taui = seq(.1,15,.1))
#> Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
forecast2 <- rseismTLS::forecast.seism(data, prior2, forecast.twin, method = 'bayesFull')

res <- hist(seism$t, breaks = forecast.bins, main = 'Pseudo-prospective forecast (bayesFull)', col = 'grey', border = 'white', ylim = c(0, 70))
polygon(x = c(forecast2$ti, rev(forecast2$ti)), y = c(forecast2$N.CI[, 1], rev(forecast2$N.CI[, 2])), col = rgb(1, 0, 0,0.5), border = NA)
boolean90 <- res$counts >= forecast2$N.CI[, 1] & res$counts <= forecast2$N.CI[, 2]
points(res$mids[boolean90], rep(-1, length(which(boolean90))), col = 'darkgreen', pch = 15, cex = .5)
points(res$mids[!boolean90], rep(-1, length(which(!boolean90))), col = 'darkred', pch = 15, cex = .5)
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" />

```r

( length(which(boolean90)) / length(res$mids) )
#> [1] 0.8732394
```



#### Maximum magnitude forecasting

TO BE COMPLETED.


## Risk functions

TO BE COMPLETED.

## References

Broccardo M., Mignan A., Wiemer S., Stojadinovic B., Giardini D. (2017), Hierarchical Bayesian Modeling of Fluid‐Induced Seismicity. Geophysical Research Letters, 44 (22), 11,357-11,367, doi: 10.1002/2017GL075251 

Dinske C., Shapiro S.A. (2013), Seismotectonic state of reservoirs inferred from magnitude distributions of fluid-induced seismicity. J. Seismol., 17, 13-25, doi: 10.1007/s10950-012-9292-9

Häring M.O., Schanz U., Ladner F., Dyer B.C. (2008), Characterisation of the Basel 1 enhanced geothermal system, Geothermics, 37 (5), 469-495, doi: 10.1016/j.geothermics.2008.06.002

Kraft T., Deichmann N. (2014), High-precision relocation and focal mechanism of the injection-induced seismicity at the Basel EGS.
Geothermics, 52, 59–73, doi: 10.1016/j.geothermics.2014.05.014

Lewis P.A.W., Shedler G.S. (1979), Simulation of nonhomogeneous Poisson processes by thinning. Naval Res. Logistics Q., 26, 403–413, doi:10.1002/nav.3800260304

Mignan A., Karnouvis D., Broccardo M., Wiemer S., Giardini D. (2019), Including seismic risk mitigation measures into the Levelized Cost Of Electricity in enhanced geothermal systems for optimal siting. Applied Energy, 238, 831-850, doi: 10.1016/j.apenergy.2019.01.109

Mignan A., Broccardo M., Wiemer S., Giardini D. (2019), Autonomous Decision-Making Against Induced Seismicity in Deep Fluid Injections. In: Ferrari A., Laloui L. (eds), Energy Geotechnics, SEG 2018, Springer Series in Geomechanics and Geoengineering, , 369-376, doi: 10.1007/978-3-319-99670-7_46

Mignan A., Broccardo M., Wiemer S., Giardini D. (2017), Induced seismicity closed-form traffic light system for actuarial decision-making during deep fluid injections. Scientific Reportsvolume, 7, 13607, doi: 10.1038/s41598-017-13585-9

Mignan A., Landtwing D., Kaestli P., Mena B., Wiemer S. (2015), Induced seismicity risk analysis of the 2006 Basel, Switzerland, Enhanced Geothermal System project: Influence of uncertainties on risk mitigation. Geothermics, 53, 133-146, doi: 10.1016/j.geothermics.2014.05.007

Ogata Y. (1988), Statistical Models for Earthquake Occurrences and Residual Analysis for Point Processes. J. Am. Stat. Assoc., 83 (401), 9-27

van der Elst N.J., Page M.T., Weiser D.A., Goebel T.H.W., Hosseini S.M. (2016), Induced earthquake magnitudes are as large as (statistically) expected. J. Geophys. Res., 121 (6), 4575-4590, doi: 10.1002/2016JB012818

