dr.ctseff <- function (y, a, x, bw.seq, a.vals, mu, g, se=FALSE, limited.mem=FALSE, verbose=TRUE) {
  require(KernSmooth)
  kern <- function(t) {
    dnorm(t)
  }
  
  x <- as.matrix(x)
  ord <- order(a)
  y <- y[ord]
  
  a <- a[ord]
  n <- length(y)
  if(limited.mem) {
    x <- data.frame(x)
    if(verbose) cat("\nComputing conditional means... \n")
    muhat.obs <- mu(a, x)
    mhat.obs <- sapply(1:n, function(i) mean(mu(rep(a[i], n), x)))
  } else {
    x.new <- data.frame(x[rep(1:n, n), ])
    x <- data.frame(x)
    # colnames(x) <- colnames(x.new)
    colnames(x.new) <- colnames(x)
    a.new <- rep(a, each = n)
    if(verbose) cat("\nComputing conditional means... \n")
    muhat.mat <- matrix(mu(a.new, x.new), byrow = TRUE, nrow = n, ncol = n)
    muhat.obs <- diag(muhat.mat)
    mhat.obs <- rowMeans(muhat.mat)
    mhat.mat <- matrix(mhat.obs, byrow = TRUE, ncol = n, nrow = n)
  }
  if(verbose) cat("\nComputing propensities... \n")
  ghat.obs <- g(a, x)
  pseudo.out <- (y - muhat.obs)/(ghat.obs) + mhat.obs
  hatvals <- function(bw) {
    w.avals <- NULL
    for (a.val in a) {
      a.std <- (a - a.val)/bw
      kern.std <- kern(a.std)/bw
      w.avals <- c(w.avals, mean(a.std^2 * kern.std) *
                     (kern(0)/bw)/(mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2))
    }
    return(w.avals/n)
  }
  cts.eff.fn <- function(out, bw) {
    approx(locpoly(a, out, bandwidth = bw), xout = a)$y
  }
  risk.fn <- function(h) {
    hats <- hatvals(h)
    mean(((pseudo.out - cts.eff.fn(pseudo.out, bw = h))/(1 - hats))^2, na.rm=TRUE)
  }
  if(verbose) cat("\nSelecting bandwidth... \n")
  risk.est <- sapply(bw.seq, risk.fn)
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  if(verbose) cat("\nEstimating dose-response curve... \n")
  est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
  res <- data.frame(a.vals, est)
  if(se) {
    se <- NULL
    if(verbose) {
      cat("\nCalculating standard errors... \n")
      pts <- (1:10) * 10
      pcts <- 100 * (1:length(a.vals)) / length(a.vals)
      a.val.pts <- sapply(pts, function(pt) a.vals[min(which(pcts >= pt))])
    }
    for (a.val in a.vals) {
      a.std <- (a - a.val)/h.opt
      kern.std <- kern(a.std)/h.opt
      beta <- coef(lm(pseudo.out ~ a.std, weights = kern.std))
      Dh <- matrix(c(mean(kern.std), mean(kern.std * a.std),
                     mean(kern.std * a.std), mean(kern.std * a.std^2)), nrow = 2)
      if(!limited.mem) kern.mat <- matrix(kern((a - a.val)/h.opt)/h.opt, byrow = FALSE, ncol = n, nrow = n)
      g2 <- (a - a.val)/h.opt
      if(limited.mem) {
        int1 <- int2 <- rep(NA, n)
        for(i in 1:n) {
          x.new <- as.data.frame(matrix(as.numeric(x[i,]), nrow=n, ncol=ncol(x), byrow=TRUE))
          names(x.new) <- names(x)
          muhat.i <- mu(a, x.new)
          int1[i] <- mean(kern.std * (muhat.i - mhat.obs))
          int2[i] <- mean(g2 * kern.std * (muhat.i - mhat.obs))
        }
      } else {
        int1 <- colMeans(kern.std * (muhat.mat - mhat.obs))
        int2 <- colMeans(g2 * kern.std * (muhat.mat - mhat.obs))
      }
      inf.fn <- t(solve(Dh) %*% rbind(kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int1,
                                      a.std * kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int2))
      sigma <- cov(inf.fn)
      se <- c(se, sqrt(sigma[1, 1]))
      if(verbose) {
        if(a.val %in% a.val.pts) cat(pts[which(a.val.pts == a.val)], "%... ")
      }
    }
    ci.ll <- est - 1.96 * se/sqrt(n)
    ci.ul <- est + 1.96 * se/sqrt(n)
    res$se <- se
    res$ci.ll <- ci.ll
    res$ci.ul <- ci.ul
  }
  return(list(res = res, bw.risk = bw.risk, h.opt=h.opt))
}