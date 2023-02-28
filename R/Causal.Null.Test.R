causalNullTest <- function(Y, A, W, p=2, control = list()) {

  if(!is.data.frame(W)) W <- as.data.frame(W)

  call <- match.call(expand.dots = TRUE)
  control <- do.call("causalNullTest.control", control)

  .check.input(Y=Y, A=A, W=W, p=p, control=control)

  library(mvtnorm)
  n <- length(Y)
  a.vals <- sort(unique(A))

  if(control$cross.fit & is.null(control$folds)) {
    control$folds <- sample(rep(1:control$V, length.out = n), replace=FALSE)
  }

  if(is.null(control$mu.hat)) {
    library(SuperLearner)
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating outcome regressions...")
      control$mu.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("fold", v, "...")
        if(length(setdiff(Y, c(0,1))) == 0) {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(A, W)[control$folds != v,], SL.library = control$mu.SL.library, family = 'binomial', method = 'method.NNloglik')
        } else {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(A, W)[control$folds != v,], SL.library = control$mu.SL.library, family = 'gaussian', method = 'method.NNLS')
        }
        function(a, w) c(predict(fit, newdata = cbind(A = a, w),onlySL = TRUE)$pred)
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating outcome regression...")
      if(length(setdiff(Y, c(0,1))) == 0) {
        mu.fit <- SuperLearner(Y = Y, X = cbind(A, W), SL.library = control$mu.SL.library, family = 'binomial', method = 'method.NNloglik')
      } else {
        mu.fit <- SuperLearner(Y = Y, X = cbind(A, W), SL.library = control$mu.SL.library, family = 'gaussian', method = 'method.NNLS')
      }
      control$mu.hat <- function(a, w) c(predict(mu.fit, newdata = cbind(A = a, w),onlySL = TRUE)$pred)
      if(control$verbose) cat("\n")
    }
  }

  if(is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating propensities..")
      control$g.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("fold", v)
        fit <- cmdSuperLearner(A = A[control$folds != v], W = W[control$folds != v,,drop=FALSE], newA = A[control$folds == v], newW = W[control$folds == v,,drop=FALSE], control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose, saveFitLibrary = FALSE))
        c(fit$SL.densities)
        #function(a, w) c(predict.cmdSuperLearner(fit, newA = a, newW = w))
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating propensity...")
      g.fit <- cmdSuperLearner(A = A, W = W, control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose, saveFitLibrary = FALSE))
      control$g.hat <- g.fit$SL.densities
      rm(g.fit)
      #control$g.hat <- function(a, w) c(predict.cmdSuperLearner(g.fit, newA = a, newW = w))

      if(control$verbose) cat("\n")
    }
  }

  if(control$verbose) cat("Computing Omega...")
  if(!control$cross.fit) {
    ord <- order(A)
    A <- A[ord]
    Y <- Y[ord]
    W <- W[ord,,drop=FALSE]
    if(inherits(control$g.hat, "function")) {
      g.hats <- control$g.hat(A, W)
      control$g.hat <- NULL
    }
    else g.hats <- control$g.hat
    if(any(g.hats < control$g.trunc)) {
      warning("Truncating g.hats below. Possible positivity issues.")
      g.hats[g.hats < control$g.trunc] <- control$g.trunc
    }
    a.ecdf <- ecdf(A)
    a.weights <- sapply(a.vals, function(a0) mean(A == a0))
    A.a.val <- sapply(A, function(a0) which(a.vals == a0))
    u.vals <- a.ecdf(a.vals)
    mu.hats.a.vals <- sapply(a.vals, function(a0) control$mu.hat(a0, W)) #rows index W, columns index a.vals
    control$mu.hat <- NULL
    mu.hats <- mu.hats.a.vals[,A.a.val]
    theta.a.vals <- colMeans(mu.hats.a.vals)
    theta.A <- theta.a.vals[A.a.val]
    mu.hats.data <- diag(mu.hats)
    partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / n
    gamma.hat <- mean(mu.hats)
    Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A <= a0) * theta.A)) - gamma.hat * u.vals

    IF.vals <- sapply(a.vals, function(a0) {
      if(any(A <= a0))  mumean.vals <- partial.mu.means[,max(which(A <= a0))]
      else mumean.vals <- 0
      (as.numeric(A <= a0) - a.ecdf(a0)) * ((Y - mu.hats.data) / g.hats + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,n] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
    })

    Omega.hat <- colMeans(IF.vals) + Omega.a.vals

    if(control$verbose) cat("\nComputing covariance...\n")

    Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
      mean(IF.vals[,s] * IF.vals[,t])
    }))
  }
  else {

    fold.Omega.hats <- matrix(NA, nrow = control$V, ncol = length(a.vals))
    IF.vals <- vector(length=control$V, mode='list')
    for(j in 1:control$V) {
      if(control$verbose) cat("fold", j, "...")
      Nv <- sum(control$folds == j)
      A.test <- A[control$folds == j]
      Y.test <- Y[control$folds == j]
      W.test <- W[control$folds == j,, drop=FALSE]
      ord <- order(A.test)
      A.test <- A.test[ord]
      Y.test <- Y.test[ord]
      W.test <- W.test[ord,, drop=FALSE]
      if(inherits(control$g.hat[[j]], "function")) {
        g.hats.test <- control$g.hat[[j]](a = A.test, w = W.test)
      }
      else g.hats.test <- control$g.hat[[j]]

      if(any(g.hats.test < control$g.trunc)) {
        warning("Truncating g.hats below. Possible positivity issues.")
        g.hats.test[g.hats.test < control$g.trunc] <- control$g.trunc
      }
      a.ecdf <- ecdf(A.test)
      a.weights <- sapply(a.vals, function(a0) mean(A.test == a0))
      A.a.val <- sapply(A.test, function(a0) which(a.vals == a0))
      u.vals <- a.ecdf(a.vals)
      mu.hats.a.vals <- sapply(a.vals, function(a0) control$mu.hat[[j]](a=a0, w=W.test)) #rows index W, columns index a.vals
      mu.hats <- mu.hats.a.vals[,A.a.val]
      theta.a.vals <- colMeans(mu.hats.a.vals)
      theta.A <- theta.a.vals[A.a.val]
      mu.hats.data <- diag(mu.hats)
      partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / Nv
      gamma.hat <- mean(mu.hats)
      Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A.test <= a0) * theta.A)) - gamma.hat * u.vals

      IF.vals[[j]] <- sapply(a.vals, function(a0) {
        if(any(A.test <= a0)) mumean.vals <- partial.mu.means[,max(which(A.test <= a0))]
        else mumean.vals <- 0
        (as.numeric(A.test <= a0) - a.ecdf(a0)) * ((Y.test - mu.hats.data) / g.hats.test + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,ncol(partial.mu.means)] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
      })

      fold.Omega.hats[j,] <- colMeans(IF.vals[[j]]) + Omega.a.vals
    }

    Omega.hat <- colMeans(fold.Omega.hats)
    if(control$verbose) cat("\nComputing covariance...\n")
    Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
      mean(unlist(lapply(IF.vals, function(IF) mean(IF[,s] * IF[,t]))))
    }))
  }

  if(control$verbose) cat("Simulating paths...\n")

  paths <- rmvnorm(control$n.sim, sigma=Sigma.hat)

  if(control$verbose) cat("Computing statistics...\n")

  a.weights <- sapply(a.vals, function(a) mean(A == a))
  ret <- t(sapply(p, function(pp) {
    stat <- ifelse(pp < Inf, (sum(abs(Omega.hat )^pp * a.weights))^{1/pp}, max(abs(Omega.hat)))

    if(pp < Inf) {
      stats <- (apply(abs(paths)^pp, 1, function(row) sum(row * a.weights)))^{1/pp}
    } else {
      stats <- apply(abs(paths), 1, max)
    }

    p.val <- mean(stats / sqrt(n) > stat)

    q <- quantile(stats, (1 - (1-control$conf.level) / 2))
    ci.ll <- max(stat - q / sqrt(n), 0)
    ci.ul <- stat + q / sqrt(n)

    res <- c(stat, p.val, ci.ll, ci.ul)

    res

  }))
  ret.df <- data.frame(p = p, obs.stat = ret[,1], p.val = ret[,2], ci.ll = ret[,3], ci.ul = ret[,4])
  ret.list <- list(test = ret.df)
  if(control$return.Omega) {
    ret.list <- data.frame(c(ret.list, list(Omega.hat = data.frame(a=a.vals, Omega.hat), IF.vals = IF.vals, paths = paths)))
  }
  if(control$save.nuis.fits) {
    ret.list <- data.frame(c(ret.list, mu.hat = control$mu.hat, g.hat = control$g.hat))
    if(control$cross.fit) ret.list <- c(ret.list, folds = control$folds)
  }

  return(ret.list)
}

causalNullTest.control <- function(mu.SL.library = NULL,
                                   g.SL.library = NULL,
                                   g.n.bins = 2:(length(unique(A))/50),
                                   cross.fit = TRUE,
                                   V = 10,
                                   folds = NULL,
                                   save.nuis.fits = FALSE,
                                   mu.hat = NULL,
                                   g.hat = NULL,
                                   g.trunc = .001,
                                   n.sim = 1e4,
                                   return.Omega = FALSE,
                                   conf.level = .95,
                                   verbose = FALSE) {
  list(mu.SL.library = mu.SL.library, g.SL.library = g.SL.library, g.n.bins = g.n.bins, cross.fit = cross.fit, V = V,
       folds = folds, mu.hat = mu.hat,  g.hat = g.hat, g.trunc = g.trunc, n.sim = n.sim, return.Omega = return.Omega,
       save.nuis.fits = save.nuis.fits, conf.level = conf.level, verbose = verbose)
}

.check.input <- function(Y, A, W, p, control) {
  if(length(Y) != length(A)) stop("Y and A must have the same length")
  if(!is.null(control$mu.hat)) {
    if(control$cross.fit) {
      if(is.null(control$folds)) {
        stop("mu.hat provided and cross.fit=TRUE, but folds not provided.")
      }
      if(length(control$mu.hat) < length(unique(control$folds))) {
        stop("mu.hat provided and cross.fit=TRUE, but mu.hats is not a list with the same length as the number of folds.")
      }
    } else {
      if(!is.null(control$folds)) {
        stop("mu.hat provided and cross.fit=FALSE, but folds were provided.")
      }
    }
  }
  if(!is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(is.null(control$folds)) {
        stop("g.hat provided and cross.fit=TRUE, but folds not provided.")
      }
      if(length(control$g.hat) < length(unique(control$folds))) {
        stop("g.hat provided and cross.fit=TRUE, but g.hats is not a list with the same length as the number of folds.")
      }
    } else {
      if(!is.null(control$folds)) {
        stop("g.hat provided and cross.fit=FALSE, but folds were provided.")
      }
    }
  }
  if(is.null(control$mu.SL.library)) {
    if(is.null(control$mu.hats)) {
      stop("mu.hats must be provided if mu.SL.library is not specified.")
    }
  } else {
    if(is.null(control$mu.hats)) {
      if(control$cross.fit & is.null(control$folds) & is.null(control$V)) {
        stop("cross.fit = TRUE, but number of folds not specified.")
      }
    }
  }
  if(is.null(control$g.SL.library)) {
    if(is.null(control$g.hats)) {
      stop("g.hats must be provided if g.SL.library is not specified.")
    }
  } else {
    if(is.null(control$g.hats)) {
      if(control$cross.fit & is.null(control$folds) & is.null(control$V)) {
        stop("cross.fit = TRUE, but number of folds not specified.")
      }
    }
  }
  if(any(is.na(Y) | is.na(A) | is.na(W))) {
    stop("Missing outcome, treatment, or confounders detected; missing data not allowed.")
  }
}
