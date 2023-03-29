# ctsDoseResponse

### Overview

ctsDoseResponse is a tool for nonparametric estimation and inference of a continuous dose-response curve, as well as for null hypothesis testing of whether a dose-response is flat or not. It has two main functions ``ctsCausal`` and ``ctsCausalTest``.

The function for nonparametric estimation and inference, ``ctsCausal``, is implemented with three different methods. The method ``dr.isoreg`` is for [causal isotonic regression](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12372), the method ``dr.loclin`` is for [doubly robust estimation of continuous treatment effects](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12212) and the method ``dr.debiased`` is for debiased doubly robust estimation of continuous treatment effects.

``ctsCausalTest`` has two different methods for null hypothesis testing. ``dr.nd`` is for [nonparametric tests of causal null with nondiscrete exposures](https://www.tandfonline.com/doi/abs/10.1080/01621459.2020.1865168?journalCode=uasa20) and ``dr.cts`` is for [nonparametric doubly robust test for a continuous treatment effect](https://arxiv.org/abs/2202.03369).


### Installation

Please use the R command ``devtools::install_github("https://github.com/Yujian-Wu/ctsDoseResponse/tree/main", dependencies = T)`` to install the package and dependent packages. You might be prompted to updated some existing packages. If you wish to update the package `RcppArmadillo` and the version of your R is not most updated, please be advised to install it from sources that the package does **NOT** need compilation.


### Example
```
library(ctsDoseResponse)

set.seed(12345)

n <- 200
cols <- 3

W <- matrix(runif(n*cols), ncol = cols) # a 200 * 3 matrix of covariates
A <- rnorm(n, mean=W%*%rnorm(cols)) # a 200 * 1 vector of treatment variable
Y <- rnorm(n, mean = a^2 + rnorm(n)) # a 200 * 1 vector of response variable

### causal isotonic regression
iso.reg <- ctsCausal(Y, A, W, method='dr.isoreg', cross.fit = F)

### continuous dose-response curve with smoothing treatment
loc.lin <- ctsCausal(Y, A, W, method='dr.loclin')

### debiased continuous dose-response curve with smoothing treatment
bc.loc.lin <- ctsCausal(Y, A, W, method='dr.debiased')

### null hypothesis test for non-diescrete treatment
test.nd <- ctsCausalTest(Y, A, W, method='dr.nd')

### null hypothesis test for continuous treatment
test.cts <- ctsCausalTest(Y, A, W, method='dr.cts')
```
