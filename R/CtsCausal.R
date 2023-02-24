#' Nonparametric estimation and inference for continuous causal dose response curves
#'
#' This function performs estimation and inference of the causal dose-response curve E[Y | A] for a expose A. The exposure may be discrete, continuous, or an arbitrary mixture of discrete and continuous components.  See the accompanying paper for details.
#'
#' @import Rsolnp
#' @import SuperLearner
#' @import sets
#' @import fdrtool
#' @import KernSmooth
#' @import reshape
#' @import plyr
#' @import nprobust
#'
#'
#' @param Y \code{n x 1} numeric vector of observed outcome values.
#' @param A \code{n x 1} numeric vector of exposure values.
#' @param W \code{n x p} data.frame of potential confounders.
#' @param method controls which method to be used for the estimation and inference. Current supported methods are: \code{'loclin'}, the method by Kennedy et al. (2017) for a continuous treatment \code{A}; \code{'isoreg'}, the method by (Westling et al., 2020) for a isotonic regression typed dose response curve; \code{'debiased'}, the method of debiased local linear estimation
#' @param SL.libraries the libraries used for Superlearner. \code{c("SL.mean", "SL.glm", "SL.gam", "SL.earth")} will be used by default.
#' @param cross.fit boolean variable that controls whether cross fitting will be used for the estimation of nuisance parameters
#' @param num.folds number that controls the folds of cross fitting. Default is 10.
#' @param sigma.sq number that will be used for plug-in estimation for inference. Default is NULL.
#' @param binary.outcome boolean variable. If set to TRUE, the output will be binary, otherwise the output will be real numbers.
#' @param verbose if set to TRUE, the function will output the progression of estimation and inferece
#' @param infer.grid numeric vector on which the confidence interval will be constructed. Default is \code{seq(min(A), max(A), length.out=50)}
#' @param bw.seq numeric vector on which the bandwidth will be selected for method \code{'loclin'} and \code{'debiased'}. Default is \code{seq((max(A) - min(A)) / length(A), max(A), length.out=50)}
#' @param conf.level number that controls the confidence level for inference. Default set to 0.95.
#' @param rho real number that controls the ratio of the bandwidth of smoothing kernel and the bandwidth of the kernel for the estimation of bias. Default is 1.
#' @export ctsCausal
#' @examples
#' # Sample data
#' n <- 1000
#' W <- data.frame(W1 = runif(n))
#' Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W$W1)))
#' A <- (1-Z) * rnorm(n, mean = W$W1, sd = abs(1 + W$W1))
#' Y <- rexp(n, rate = 1+abs(W$W1 * A))
#' ctsCausal(Y, A, W, method='isoreg')

library(Rsolnp)
library(SuperLearner)
library(sets)
library(fdrtool)
library(KernSmooth)
library(reshape)
library(plyr)
library(nprobust)

ctsCausal <- function(Y, A, W, method,
                      SL.libraries=c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                      cross.fit=TRUE,
                      num.folds = 10,
                      sigma.sq=NULL,
                      binary.outcome=FALSE,
                      verbose=TRUE,
                      infer.grid=seq(min(A), max(A), length.out=50),
                      bw.seq=seq((max(A) - min(A)) / length(A), max(A), length.out=50),
                      conf.level=0.95,
                      rho=1,
                      kernel='epa'
){
  if(!is.data.frame(W)) {
    W <- data.frame(W)
    colnames(W) <- sapply(1:ncol(W), function(u) paste('W',u, sep = ''))
    }
  #### Isotonic regression and inference
  if(tolower(method) == 'isoreg'){
    ##### Causal Isotonic Regression
      nuisance <- causalDoseResponse(Y, A, W, control = list(cross.fit=cross.fit,
                                                  verbose=verbose,
                                                  V=num.folds,
                                                  mu.SL.library=SL.libraries,
                                                  g.SL.library=SL.libraries))
    #### Isotonic inference
    if (cross.fit){
      mono.curve <- causal.isoreg(Y, A, W, g.hats = nuisance$g.hat, mu.hats = nuisance$mu.predicted, mu.means = nuisance$m.means.predicted,
                                  binary.outcome, verbose = verbose, folds = nuisance$folds)
    }else{
      mono.curve <- causal.isoreg(Y, A, W, g.hats = nuisance$g.hat, mu.hats = nuisance$mu.predicted, mu.means = nuisance$m.means.predicted,
                                  binary.outcome, verbose = verbose)
    }

    mono.infer <- confint.causal.isoreg(Y, A, W, fit = mono.curve,
                                        scale.type=c('plug.in', 'DR'),
                                        mu.hats=nuisance$mu.predicted,
                                        g.hats=nuisance$g.hat,
                                        mu=nuisance$mu.hat,
                                        g=nuisance$g.hat.fun,
                                        mu.means=nuisance$m.means.predicted,
                                        x0.vals=infer.grid,
                                        conf=conf.level,
                                        sigma.sq=sigma.sq, verbose=verbose,binary.outcome)

    ### return the estimates
    if(is.null(mono.infer$binary.plug.in)){
      return(list(dose=mono.curve$x.vals, response=mono.curve$theta.hat, DRCI=mono.infer$dr.CI,
                  PlugCI=mono.infer$plug.in.CI))
    }else{
      return(list(dose=mono.curve$x.vals, response=mono.curve$theta.hat, DRCI=mono.infer$dr.CI,
                  PlugCI=mono.infer$plug.in.CI, PlugBiCI=mono.infer$binary.plug.in))
    }

  }else{
      ##### Causal Regression for nuisance parameters
    nuisance <- causalDoseResponse(Y, A, W, control = list(cross.fit=FALSE,
                                                    verbose=verbose,
                                                    V=num.folds,
                                                    mu.SL.library=SL.libraries,
                                                    g.SL.library=SL.libraries))

    if(tolower(method) == 'loclin'){
    #### Kennedy's continuous method

    loclin <- dr.ctseff(Y, A, W, bw.seq = bw.seq, a.vals = infer.grid,
                         mu = nuisance$mu.hat, g = nuisance$g.hat.fun, limited.mem = F, se = T)

    return(list(dose=loclin$res$a.vals, response=loclin$res$est, DRCI=data.frame(ci.ll=loclin$res$ci.ll,
                                                                              ci.ul=loclin$res$ci.ul),
                h=loclin$h.opt))
    }else if(tolower(method) == 'debiased'){
      #### Kenta's debiased method
      bc.reg <- debiased.ctseff(y=Y, a=A, x=W, bw.seq = bw.seq, eval.pts = infer.grid, mu = nuisance$mu.hat,
                                g = nuisance$g.hat.fun, tau = rho, kernel.type=kernel, verbose = verbose)

      bc.est <- bc.reg$mu - bc.reg$b

      return(list(dose=bc.reg$x, response=bc.est, IFCI=data.frame(ci.ll=bc.est - qnorm(1-(1-conf.level)/2)*bc.reg$se.infl.robust,
                                                                ci.ul=bc.est + qnorm(1-(1-conf.level)/2)*bc.reg$se.infl.robust),
             RBCI=data.frame(ci.ll=bc.est - qnorm(1-(1-conf.level)/2)*bc.reg$se.rb,
                             ci.ul=bc.est + qnorm(1-(1-conf.level)/2)*bc.reg$se.rb)))
    }
  }
}


#' Nonparametric test of flatness of a causal dose-response curve under exogeneity
#'
#' This function performs a hypothesis test that the causal dose-response curve theta(a) is flat on the support of the observed exposure A. It includes two different method for the test. The method by (Westling et al., 2020) allows the exposure may be discrete, continuous, or an arbitrary mixture of discrete and continuous components, while the method by (Weng et al., 2022) only considers the continuous treatment.  See the accompanying paper for details.
#'
#' @param Y \code{n x 1} numeric vector of observed outcome values.
#' @param A \code{n x 1} numeric vector of exposure values.
#' @param W \code{n x p} data.frame of potential confounders.
#' @param method decides which method to be used for the test. Current supported methods are: \code{'isoreg'}, the method by (Westling et al., 2020); \code{'continuous'}, the method by (Weng et al., 2022) for a continuous treatment curve.
#' @param SL.libraries the libraries used for Superlearner. If set to NULL, then the default libraries \code{c("SL.mean", "SL.glm", "SL.gam", "SL.earth")} will be used.
#' @param cross.fit boolean variable that controls whether cross fitting will be used for the estimation of nuisance parameters
#' @param conf.level the confidence level of the test Default set to 0.95.
#' @param dist the distribution used by (Weng et al., 2022) for the test. It has to be either \code{"TwoPoint"} or \code{"Rademachar"}.
#' @param verbose if set to TRUE, the function will output the progression of estimation and inferece
#' @return \code{causalNullTest} returns a list with the following elements:
#' \item{p.value}{The p value of the test.}
#' \item{test.stat}{Test statistic.}
#' If \code{method == "isoreg"}, the following elements are also included in the output:
#' \item{stat.ci.ll}{The lower limit of the confidence interval for the test statistic.}
#' \item{stat.ci.ul}{The upper limit of the confidence interval for the test statistic.}
#' @export ctsCausalTest
#' @examples
#' # Sample data
#' n <- 1000
#' W <- data.frame(W1 = runif(n))
#' Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W$W1)))
#' A <- (1-Z) * rnorm(n, mean = W$W1, sd = abs(1 + W$W1))
#' Y <- rexp(n, rate = 1+abs(W$W1 * A))
#' ctsCausalTest(Y, A, W, method='isoreg', cross.fit=T, verbose=TRUE)


ctsCausalTest <- function(Y, A, W, method, conf.level=0.95, dist="TwoPoint", cross.fit = TRUE, SL.library=NULL, verbose=TRUE){

  if(is.null(SL.library)) SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth") else SL.library <- SL.library

  if(tolower(method) == 'isoreg'){
    test <- causalNullTest(Y, A, W, control = list(mu.SL.library=SL.library, g.SL.library=SL.library, cross.fit=cross.fit, verbose=verbose, conf.level=conf.level))

    return(list(p.value=test$test$p.val, test.stat=test$test$obs.stat, stat.ci.ll=test$test$ci.ll, stat.ci.ul=test$test$ci.ul))
  }else if(tolower(method) == 'continuous'){
    test.fit <- causalDoseResponse(Y, A, W, control = list(method = 'loclin', cross.fit=F))
    test <- drdrtest(Y, A, W, arange = c(min(A), max(A)), pifunc = test.fit$g.hat.fun, mufunc = test.fit$mu.hat, cross.fit=cross.fit, dist=dist)
    return(list(p.value=test$p.value, test.stat=test$test.stat, stat.ci.ll=NULL, stat.ci.ul=NULL))
  }
}

