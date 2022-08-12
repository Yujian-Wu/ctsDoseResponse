causal.isoreg <- function(Y, X, W, g.hats, mu.hats, mu.means, binary.outcome=FALSE, verbose=FALSE, folds=NULL) {
  
  n <- length(X)
  if(class(W) == "data.frame" | class(W) == "numeric") W <- as.matrix(W)
  
  if(verbose) cat("\nEstimating theta...\n")
  
  x.ecdf <- ecdf(X)
  U <- x.ecdf(X)
  x.vals <- unique(sort(X))
  u.vals <- x.ecdf(x.vals)
  
  if(is.null(folds)) {
    term1 <- (Y - mu.hats) / g.hats
    step.fun.vals <- sapply(u.vals, function(u0) sum(term1[U <= u0])/n)
    
    int.vals <- sapply(u.vals, function(u0) mean(mu.means * as.numeric(U <= u0)))
    
    Gamma.hat <- step.fun.vals + int.vals
  } else {
    Gamma.hats <- sapply(unique(folds), function(v) {
      inds <- which(folds == v)
      n.v <- length(inds)
      term1 <- (Y[inds] - mu.hats[inds]) / g.hats[inds]
      step.fun.vals <- sapply(u.vals, function(u0) sum(term1[U[inds] <= u0])/n.v)
      
      int.vals <- sapply(u.vals, function(u0) mean(mu.means[inds] * as.numeric(U[inds] <= u0)))
      
      step.fun.vals + int.vals
    })
    Gamma.hat <- rowMeans(Gamma.hats)
  }
  
  gcm <- gcmlcm(c(0,u.vals), c(0,Gamma.hat), type='gcm')
  u.knots <- gcm$x.knots
  Psi.bar <- gcm$y.knots
  slope.knots <- gcm$slope.knots
  if(binary.outcome) slope.knots <- pmin(pmax(slope.knots, 0), 1)
  
  # values of thetahat at x
  theta.hat.uvals <- sapply(u.vals, function(u0) slope.knots[max(which(u.knots < u0))])
  
  ret <- list(x.vals=x.vals, theta.hat=theta.hat.uvals, Gamma.hat=Gamma.hat, u.knots=u.knots, slope.knots=slope.knots)
  
  if(verbose) cat('\nEstimation done!\n')
  return(ret)
}



confint.causal.isoreg <- function(Y, X, W, fit, mu.hats, g.hats, mu, g, mu.means, scale.type=c('plug.in', 'DR'), x0.vals=seq(min(X), max(X), length.out=50), conf=.95, sigma.sq=NULL, verbose=FALSE,  binary.outcome=FALSE) {
  cat("\nEstimating confidence intervals...\n")
  n <- length(X)
  ret <- NULL
  ret$x.vals <- x0.vals
  x.ecdf <- ecdf(X)
  ret$theta.hat <- sapply(x.ecdf(x0.vals), function(u0) fit$slope.knots[max(which(fit$u.knots < u0))])
  ## INFERENCE
  if(conf > .9999) stop("Confidence greater than 0.9999 not allowed.")
  # Get chernoff quantile
  pp <- 1 - (1 - conf)/2
  if(pp %in% chernoff_quantiles$p) {
    quantile <- chernoff_quantiles$Finv[which(chernoff_quantiles$p == pp)]
  }else {
    up_p <- min(chernoff_quantiles$p[chernoff_quantiles$p > pp])
    down_p <- max(chernoff_quantiles$p[chernoff_quantiles$p < pp])
    up_F <- chernoff_quantiles$Finv[chernoff_quantiles$p == up_p]
    down_F <- chernoff_quantiles$Finv[chernoff_quantiles$p == down_p]
    quantile <- down_F + ((up_F - down_F)/(up_p - down_p)) * (pp2 - down_p)
  }
  
  # Estimate derivative
  if(verbose) cat("\nEstimating derivative...\n")
  deriv.hat <- estimate.deriv(u.knots=fit$u.knots, slope.knots=fit$slope.knots, n=length(X))
  derivs <- approxfun(deriv.hat$u, deriv.hat$psi.prime.hat)(x.ecdf(x0.vals))
  if(any(derivs <= 0, na.rm=TRUE)) {
    warning("Difficulty estimating derivstive: some derivative estimates <= 0. Replacing with small value.")
    derivs[derivs <= 0] <- min(derivs[derivs > 0])
  }
  ret$psi.prime.hat <- derivs
  
  # Estimate kappa
  if(verbose) cat("\nEstimating scale...\n")
  kappas <- estimate.kappa(Y,X, W, mu.hats=mu.hats, g.hats=g.hats, mu=mu, g=g, mu.means=mu.means, fit=fit, scale.type=scale.type, x0.vals=x0.vals, sigma.sq=sigma.sq, verbose=verbose, binary.outcome=binary.outcome)
  if("plug.in" %in% tolower(scale.type)) {
    cat("\nEstimating plug in confidence interval...\n")
    ret$plug.in.kappa <- kappas$plug.in.kappa
    plug.in.scale <-  (4 * abs(derivs) * kappas$plug.in.kappa)^{1/3}
    ret$plug.in.CI <- cbind(ret$theta.hat - plug.in.scale * quantile / n^{1/3}, ret$theta.hat + plug.in.scale * quantile / n^{1/3})
    if(binary.outcome){
      ret$plug.in.CI[,1] <- pmax(ret$plug.in.CI[,1], 0)
      ret$plug.in.CI[,2] <- pmin(ret$plug.in.CI[,2], 1)
    }
    ret$plug.in.CI <- data.frame(ci.ll=ret$plug.in.CI[,1], ci.ul=ret$plug.in.CI[,2])
  }
  
  if("binary.plug.in" %in% tolower(scale.type)) {
    ret$binary.plug.in.kappa <- kappas$binary.plug.in.kappa
    binary.plug.in.scale <-  (4 * abs(derivs) * kappas$binary.plug.in.kappa)^{1/3}
    ret$binary.plug.in.CI <- cbind(ret$theta.hat - binary.plug.in.scale * quantile / n^{1/3}, ret$theta.hat + binary.plug.in.scale * quantile / n^{1/3})
    if(binary.outcome){
      ret$binary.plug.in.CI[,1] <- pmax(ret$binary.plug.in.CI[,1], 0)
      ret$binary.plug.in.CI[,2] <- pmin(ret$binary.plug.in.CI[,2], 1)
    }
    ret$binary.plug.in.CI <- data.frame(ci.ll=ret$binary.plug.in.CI[,1], ci.ul=ret$binary.plug.in.CI[,2])
  }
  if("dr" %in% tolower(scale.type)) {
    cat("\nEstimating doubly robust confidence interval...\n")
    ret$dr.kappa <- kappas$dr.kappa
    dr.scale <-  (4 * abs(derivs) * kappas$dr.kappa)^{1/3}
    ret$dr.CI <- cbind(ret$theta.hat - dr.scale * quantile / n^{1/3}, ret$theta.hat + dr.scale * quantile / n^{1/3})
    if(binary.outcome){
      ret$dr.CI[,1] <- pmax(ret$dr.CI[,1], 0)
      ret$dr.CI[,2] <- pmin(ret$dr.CI[,2], 1)
    }
    ret$dr.CI <- data.frame(ci.ll=ret$dr.CI[,1], ci.ul=ret$dr.CI[,2])
  }
  if(verbose) cat("\n")
  cat("\nInference done!\n")
  return(ret)
}

estimate.deriv <- function(u.knots, slope.knots, n) {
  pts <- u.knots[c(-1, -length(u.knots))] # BUG fixed
  # pts <- u.knots
  
  vals <- sapply(pts, function(u0) slope.knots[max(which(u.knots < u0))]) # BUG fixed
  # vals <- sapply(pts, function(u0) ifelse(length(which(u.knots < u0)) == 0, 0, slope.knots[max(which(u.knots < u0))]) )
  
  
  bw <- try(KernSmooth::dpill(pts, vals), silent=TRUE)
  while(class(bw) == "try-error" | is.na(bw)) {
    pts <- seq(0,1, length.out=length(pts) + 1)
    vals <- sapply(pts[-c(1,length(pts))], function(u0) slope.knots[max(which(u.knots < u0))])
    vals <- c(vals[1], vals, vals[length(vals)])
    bw <- try(KernSmooth::dpill(pts, vals), silent=TRUE)
  }
  fit <- KernSmooth::locpoly(pts, vals, drv=1, bandwidth=bw)
  return(list(u=fit$x, psi.prime.hat=fit$y))
  
  
}

estimate.kappa <- function(Y, X, W, mu.hats, g.hats, mu, g, mu.means, fit, scale.type=c('plug.in', 'DR'), x0.vals, sigma.sq=NULL, kern = function(x) .75 * (1 - x^2) * (abs(x) <= 1), verbose=FALSE,  binary.outcome=FALSE) {
  n <- length(Y)
  ret <- NULL
  x.ecdf <- ecdf(X)
  U <- x.ecdf(X)
  u.vals <- unique(sort(U))
  
  x0w.df <- expand.grid.df(data.frame(X=x0.vals), data.frame(W))
  x0w.df$U <- x.ecdf(x0w.df$X)
  
  W.4pred <- as.data.frame(x0w.df[,2:(ncol(W) + 1)])
  colnames(W.4pred) <- colnames(W)
  
  ### fix the case when crossfit TRUE
  if(is.list(g) & is.list(mu)){
    x0w.df$muhat <- rowMeans(sapply(1:length(mu), function(u){c(mu[[u]](x0w.df$X, W.4pred))}))
    
    x0w.df$ghat <- rowMeans(sapply(1:length(g), function(u){c(g[[u]](x0w.df$X, W.4pred))}))
  }else{
    x0w.df$muhat <- mu(x0w.df$X, W.4pred)
    x0w.df$ghat <- g(x0w.df$X, W.4pred)
  }
  if("plug.in" %in% tolower(scale.type)) {
    if(is.null(sigma.sq)) {
      warning("Scale type 'plug.in' specified but sigma.sq not supplied; no plug-in interval computed.")
    } else {
      if(verbose) cat("scale type 'plug.in'... ")
      x0w.df$sigma.sq.hat <- sigma.sq(x0w.df$X,  x0w.df[,2:(ncol(W) + 1)])
      kappa.hat <- ddply(x0w.df, .(X), function(df) data.frame(kappa.hat = mean(df$sigma.sq.hat / ghat)))
      kappa.hat <- sapply(x0.vals, function(x0) kappa.hat$kappa.hat[kappa.hat$X == x0])
      ret$plug.in.kappa <- kappa.hat
    }
    # ret$plug.in.kappa <- kappa.hat # BUG fixed
  }
  if("binary.plug.in" %in% tolower(scale.type) & binary.outcome) {
    if(length(setdiff(unique(as.numeric(Y)), 0:1)) != 0) {
      warning("Interval type 'binary.plug.in' specified but Y not detected to be binary; no binary plug-in interval computed.")
    } else {
      if(verbose) cat("scale type 'binary.plug.in'... ")
      kappa.hat <- ddply(x0w.df, .(X), function(df) data.frame(kappa.hat = mean(df$muhat *(1 - df$muhat) / df$ghat)))
      kappa.hat <- sapply(x0.vals, function(x0) kappa.hat$kappa.hat[kappa.hat$X == x0])
      ret$binary.plug.in.kappa <- kappa.hat
    }
  }
  if("dr" %in% tolower(scale.type)) {
    if(verbose) cat("\nscale type is 'dr'... \n")
    # fold <- sample(rep(1:10, length.out=n), n, replace=FALSE)
    hs <- exp(seq(-log(2*n), log(5), length.out=50))
    theta.hat.Uvals <- sapply(U, function(u0) fit$slope.knots[max(which(fit$u.knots < u0))])
    func.vals <- c(((Y - mu.hats) / g.hats +mu.means-theta.hat.Uvals )^2)
    loss.h <- sapply(1:length(hs), function(j) {
      h <- hs[j]
      if(j %% 5 == 0 & verbose) cat(2*j, "% ")
      all.kerns <- kern(outer(U, u.vals, "-")/h) * func.vals / h
      kappa.n.h <- colMeans(all.kerns)
      gamma.n.h <- mean(kappa.n.h)
      mean((func.vals - gamma.n.h)^2)
    })
    
    h.hat <- hs[which.min(loss.h)]
    
    kappa.hat <- sapply(x0.vals, function(x0) {
      mean(kern((U - x.ecdf(x0)) / h.hat) * func.vals) / h.hat
    })
    
    ret$dr.kappa <- kappa.hat
    
  }
  return(ret)
}

chernoff_quantiles <- structure(list(p = c(0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56,
                                           0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67,
                                           0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78,
                                           0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89,
                                           0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.975, 0.98, 0.99,
                                           0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999,
                                           0.9999),
                                     Finv = c(0, 0.013187, 0.026383, 0.039595, 0.05283, 0.066096,
                                              0.079402, 0.092757, 0.106168, 0.119645, 0.133196, 0.146831, 0.16056,
                                              0.174393, 0.188342, 0.202418, 0.216633, 0.230999, 0.24553, 0.260242,
                                              0.275151, 0.290274, 0.305629, 0.321238, 0.337123, 0.353308, 0.369821,
                                              0.386694, 0.403959, 0.421656, 0.439828, 0.458525, 0.477804, 0.497731,
                                              0.518383, 0.539855, 0.562252, 0.585706, 0.610378, 0.636468, 0.664235,
                                              0.694004, 0.726216, 0.761477, 0.800658, 0.845081, 0.896904, 0.960057,
                                              0.99818, 1.04303, 1.17153, 1.18981, 1.2099, 1.23224, 1.2575,
                                              1.28666, 1.32137, 1.36464, 1.42302, 1.51666, 1.78406)), .Names = c("p",
                                                                                                                 "Finv"), row.names = c(NA, -61L), class = "data.frame")