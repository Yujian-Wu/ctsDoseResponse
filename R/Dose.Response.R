causalDoseResponse.control <- function(A, method = "loclin",
                                       var.method = "IF",
                                       mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                       g.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                       g.n.bins = 2:(length(unique(A))/50),
                                       cross.fit = TRUE,
                                       V = 10,
                                       folds = NULL,
                                       save.nuis.fits = FALSE,
                                       mu.hat = NULL,
                                       g.hat = NULL,
                                       # g.hat.fun=NULL,
                                       conf.level = 0.95,
                                       verbose = TRUE,
                                       limited.mem = length(A) > 1e4) {
  # limited.mem = length(A) > 1e1) {
  list(method = method, var.method = var.method, mu.SL.library = mu.SL.library, g.SL.library = g.SL.library, g.n.bins = g.n.bins, cross.fit = cross.fit, V = V, folds = folds, save.nuis.fits = save.nuis.fits, mu.hat = mu.hat, g.hat = g.hat, conf.level = conf.level, verbose = verbose, limited.mem = limited.mem)
}

cmdSuperLearner.control <- function (n.bins = 2:floor(length(unique(A))/50),
                                     SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                     saveFitLibrary = TRUE, verbose = FALSE) {
  list(n.bins = n.bins, SL.library = SL.library, saveFitLibrary = saveFitLibrary, verbose = verbose)
}

cmdSuperLearner.CV.control <- function (V = 10L, shuffle = TRUE, validRows = NULL) {
  V <- as.integer(V)
  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("validRows must be a list of length V containing the row numbers for the corresponding validation set")
    }
    if (!identical(V, length(validRows))) {
      stop("V and length(validRows) must be identical")
    }
  }
  list(V = V, shuffle = shuffle, validRows = validRows)
}

cmdSuperLearner.onebin <- function(A, W, newA=A, newW=W, b, SL.library, verbose, validRows, saveFitLibrary) {
  n.folds <- length(validRows)
  a.ecdf <- ecdf(A)
  U <- a.ecdf(A)
  n <- nrow(W)
  m <- nrow(newW)
  W <- as.data.frame(W)
  newW <- as.data.frame(newW)
  U <- as.numeric(U)

  tab <- table(U)
  un.U <- as.numeric(names(tab))
  un.U.frac <- as.numeric(tab) / length(U)
  if(b <= 1) stop("Number of bins must be > 1")
  if(length(un.U) < b) stop("Number of bins must not be larger than number of unique values of U.")
  if(length(un.U) == b) {
    mass.pts <- un.U
    bins <- data.frame(bin = 1:b, lower = un.U, upper = un.U, bin.length = 0, mass.pt = TRUE)
  }
  if(length(un.U) > b) {
    if(any(un.U.frac >= 1/b)) {
      mass.pts <- un.U[un.U.frac >= 1/b]
      n.mass.pts <- length(mass.pts)
      mass.pt.lowers <- round(sapply(mass.pts, function(x) max(c(U[x - U > 1/(10*n)], 0))), 7)
      mass.intervals <- lapply(1:n.mass.pts, function(j) {
        interval(mass.pt.lowers[j], mass.pts[j], bounds="(]")
      })
      cont.intervals <- data.frame(lower=c(0,mass.pts), upper=c(mass.pt.lowers, 1))
      cont.intervals$length <- cont.intervals$upper - cont.intervals$lower
      cont.intervals <- subset(cont.intervals, length > 0)
    }
    else {
      mass.pts <- NULL
      n.mass.pts <- 0
      mass.intervals <- NULL
      cont.intervals <- data.frame(lower=0, upper=1, length=1)
    }

    n.cont.bins <- b - n.mass.pts
    if(n.cont.bins > 0) {
      delta <- sum(cont.intervals$length) / n.cont.bins
      delta <- round(delta, digits=ceiling(log10(n)) + 2)
      cont.bin.endpts <- matrix(NA, nrow=n.cont.bins, ncol=2)

      for(j in 1:n.cont.bins) {
        if(j == 1) start <- cont.intervals$lower[1]
        else start <- end
        start.interval <- max(which(cont.intervals$lower <= start + 1e-6 & start <= cont.intervals$upper + 1e-6))
        if(start == cont.intervals$upper[start.interval]) {
          start.interval <- start.interval + 1
          start <- cont.intervals$lower[start.interval]
        }
        end <- start + delta
        end.interval <- start.interval
        if(!(all.equal(end, cont.intervals$upper[end.interval]) == TRUE) && end > cont.intervals$upper[end.interval]) {
          length.used <- cont.intervals$upper[end.interval] - start
          length.left <- delta - length.used
          end.interval <- end.interval + 1
          end <- cont.intervals$lower[end.interval] + length.left
        }
        while(!(all.equal(end, cont.intervals$upper[end.interval]) == TRUE) && end > cont.intervals$upper[end.interval]) {
          length.used <- length.used + cont.intervals$upper[end.interval] - cont.intervals$lower[end.interval]
          length.left <- delta - length.used
          end.interval <- end.interval + 1
          end <- cont.intervals$lower[end.interval] + length.left
        }
        end <- round(end, digits=ceiling(log10(n)) + 3)
        if(j == n.cont.bins) end <- cont.intervals$upper[nrow(cont.intervals)]
        cont.bin.endpts[j,] <- c(start, end)
      }

      cont.intervals <- lapply(1:n.cont.bins, function(j) {
        if(j == 1) int <- interval(cont.bin.endpts[j, 1], cont.bin.endpts[j, 2], bounds="[]")
        else int <- interval(cont.bin.endpts[j, 1], cont.bin.endpts[j, 2], bounds="(]")
        if(n.mass.pts > 0) {
          for(k in 1:n.mass.pts) {
            int <- interval_complement(mass.intervals[[k]], int)
          }
        }
        return(int)
      })
    } else {
      cont.intervals <- list()
    }

    bins <- c(mass.intervals, cont.intervals)
  }

  bin.sizes <- unlist(lapply(bins, interval_measure))

  disc.U <- .find.bin(U, bins)

  U.new <- a.ecdf(newA)
  disc.U.new <- .find.bin(U.new, bins)

  bin.fracs <- sapply(1:b, function(j) mean(disc.U == j))

  bin.fits <- NULL
  bin.probs <- matrix(NA, nrow=m, ncol=b)
  for(bin in 1:b) {
    if(verbose) cat("bin", bin, "...")
    capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, newX=newW, family='binomial', SL.library = SL.library, method='method.NNloglik', control = list(saveFitLibrary=saveFitLibrary), cvControl = list(V=n.folds, validRows=validRows)), silent=TRUE))
    if(class(bin.fit) == "try-error") {
      capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, newX=newW, family='binomial', SL.library = SL.library, method='method.NNLS',control = list(saveFitLibrary=saveFitLibrary), cvControl = list(V=n.folds, validRows=validRows)), silent=TRUE))
    }
    if(class(bin.fit) == "try-error") {
      capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, newX=newW, family='binomial', SL.library = SL.library, method='method.NNLS2', control = list(saveFitLibrary=saveFitLibrary), cvControl = list(V=n.folds, validRows=validRows)), silent=TRUE))
    }
    if(class(bin.fit) != "try-error") {
      bin.fits[[paste0("bin", bin, ".SL")]] <- bin.fit
      bin.probs[,bin] <- bin.fit$SL.predict
    } else {
      bin.mean <- mean(as.numeric(disc.U==bin))
      if(class(SL.library) == "character") n.algs <- length(SL.library)
      else n.algs <- sum(unlist(lapply(SL.library, function(sl) length(sl) - 1)))
      bin.fits[[paste0("bin", bin, ".SL")]] <- list(Z = matrix(bin.mean, nrow=n,ncol=n.algs), library.predict = matrix(bin.mean, nrow=m,ncol=n.algs))
      bin.probs[,bin] <- bin.mean
    }
  }



  #SL.bin.probs <- .make.doubly.stochastic(bin.probs, row.sums = rep(1, m), col.sums = bin.fracs * m)
  SL.bin.probs <- bin.probs / rowSums(bin.probs)

  SL.densities <- SL.bin.probs[cbind(1:m, disc.U.new)] / bin.sizes[disc.U.new]

  n.alg <- ncol(bin.fits[["bin1.SL"]]$Z)
  cv.library.densities <-  matrix(NA, nrow=n, ncol=n.alg)
  library.densities <- matrix(NA, nrow=m, ncol=n.alg)
  for (j in 1:n.alg) {
    cv.bin.probs <- matrix(NA, nrow = n, ncol = b)
    library.bin.probs <- matrix(NA, nrow = m, ncol = b)
    for (bin in 1:b) {
      cv.bin.probs[, bin] <- bin.fits[[paste0("bin", bin, ".SL")]]$Z[,j]
      library.bin.probs[, bin] <- bin.fits[[paste0("bin", bin, ".SL")]]$library.predict[,j]
    }
    if(any(is.na(cv.bin.probs)) | any(colSums(cv.bin.probs) == 0) | any(rowSums(cv.bin.probs) == 0)) {
      cv.library.densities[,j] <- rep(NA, n)
    } else {
      #cv.bin.probs <- .make.doubly.stochastic(cv.bin.probs, row.sums = rep(1, n), col.sums = bin.fracs * n)
      cv.bin.probs <- cv.bin.probs / rowSums(cv.bin.probs)
      cv.library.densities[,j] <- cv.bin.probs[cbind(1:n, disc.U)] / bin.sizes[disc.U]
    }

    if(any(is.na(library.bin.probs)) | any(colSums(library.bin.probs) == 0) | any(rowSums(library.bin.probs) == 0)) {
      library.densities[,j] <- rep(NA, m)
    } else {
      #library.bin.probs <- .make.doubly.stochastic(library.bin.probs, row.sums = rep(1, n), col.sums = bin.fracs * n)
      library.bin.probs <- library.bin.probs / rowSums(library.bin.probs)
      library.densities[,j] <- library.bin.probs[cbind(1:m, disc.U.new)] / bin.sizes[disc.U.new]
    }
  }

  alg.names <- paste0(bin.fits[["bin1.SL"]]$libraryNames, "_", b, "bins")

  ret <- list(bins = bins,  a.ecdf = a.ecdf, SL.bin.probs = SL.bin.probs, SL.densities = SL.densities, cv.library.densities = cv.library.densities, library.densities = library.densities, alg.names = alg.names)
  if(saveFitLibrary) ret$bin.fits <- bin.fits
  return(ret)

}

causalDoseResponse <- function(Y, A, W, control = list()) {

  if(!is.data.frame(W)) W <- data.frame(W)

  control$A <- A
  control <- do.call("causalDoseResponse.control", control)
  n <- length(Y)

  if(control$cross.fit & is.null(control$folds)) {
    control$folds <- sample(rep(1:control$V, length.out = n), replace=FALSE)
  }

  if(is.null(control$mu.hat)) {

    Fn <- ecdf(A)
    U <- Fn(A)
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating outcome regressions...")
      control$mu.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("\nfold", v, "...")
        if(length(setdiff(Y, c(0,1))) == 0) {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(U, W)[control$folds != v,], SL.library = control$mu.SL.library, family = 'binomial', method = 'method.NNloglik')
        } else {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(U, W)[control$folds != v,], SL.library = control$mu.SL.library, family = 'gaussian', method = 'method.NNLS')
        }
        function(a, w) c(predict(fit, newdata = cbind(U = Fn(a), w),onlySL = TRUE)$pred)
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating outcome regression...\n")
      if(length(setdiff(Y, c(0,1))) == 0) {
        mu.fit <- SuperLearner(Y = Y, X = cbind(U, W), SL.library = control$mu.SL.library, family = 'binomial', method = 'method.NNloglik')
      } else {
        mu.fit <- SuperLearner(Y = Y, X = cbind(U, W), SL.library = control$mu.SL.library, family = 'gaussian', method = 'method.NNLS')
      }
      control$mu.hat <- function(a, w) c(predict(mu.fit, newdata = cbind(U = Fn(a), w),onlySL = TRUE)$pred)
      if(control$verbose) cat("\n")
    }
  }

  if(is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(control$verbose) cat("\nEstimating propensities..\n")
      control$g.hat.fun <- vector(mode = "list", length = control$V)
      control$g.hat <- rep(NA, n)
      for(v in 1:control$V){
        if(control$verbose) cat("\nfold", v)
        g.fit <- cmdSuperLearner(A = A[control$folds != v], W = W[control$folds != v,,drop=FALSE], newA = A[control$folds == v], newW = W[control$folds == v,,drop=FALSE], control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose))
        # c(fit$SL.densities)
        control$g.hat[control$folds == v] <- g.fit$SL.densities
        control$g.hat.fun[[v]] <- function(a, w) c(predict.cmdSuperLearner(g.fit, newA = a, newW = w))
      }
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("\nEstimating propensity...\n")
      g.fit <- cmdSuperLearner(A = A, W = W, control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose))
      control$g.hat <- g.fit$SL.densities
      control$g.hat.fun <- function(a, w) c(predict.cmdSuperLearner(g.fit, newA = a, newW = w))
      # rm(g.fit)

      if(control$verbose) cat("\n")
    }
  }

  if(control$verbose) cat("Computing conditional means... \n")

  if(control$limited.mem) {
    if(control$cross.fit) {
      muhat.obs <- mhat.obs <- rep(NA, n)
      for(v in unique(control$folds)) {
        muhat.obs[control$folds == v] <- control$mu.hat[[v]](A[control$folds == v], W[control$folds == v,,drop=FALSE])
        mhat.obs[which(control$folds == v)] <- sapply(which(control$folds == v), function(i) {
          mean(control$mu.hat[[v]](rep(A[i], sum(control$folds == v)), W[control$folds == v,,drop=FALSE]))
        })
      }
    } else {
      muhat.obs <- control$mu.hat(A, W)
      mhat.obs <- sapply(1:n, function(i) mean(control$mu.hat(rep(A[i], n), W)))
    }
    if(!control$save.nuis.fits) control$mu.hat <- NULL
  } else {
    if(control$cross.fit) {
      muhat.obs <- mhat.obs <- rep(NA, n)
      for(v in unique(control$folds)) {
        nv <- sum(control$folds == v)
        test.set <- control$folds == v
        W.new <- data.frame(W[test.set,,drop=FALSE][rep(1:nv, nv), , drop=FALSE])
        A.new <- rep(A[test.set], each = nv)
        muhat.mat <- matrix(control$mu.hat[[v]](A.new, W.new), byrow = TRUE, nrow = nv, ncol = nv)
        muhat.obs[test.set] <- diag(muhat.mat)
        mhat.obs[test.set] <- rowMeans(muhat.mat)
        # mhat.mat <- matrix(mhat.obs, byrow = TRUE, ncol = n, nrow = n)
        # muhat.obs[control$folds == v] <- control$mu.hat[[v]](A[control$folds == v], W[control$folds == v,,drop=FALSE])
        # mhat.obs[which(control$folds == v)] <- sapply(which(control$folds == v), function(i) {
        #   mean(control$mu.hat[[v]](rep(A[i], ), W[control$folds == v,,drop=FALSE]))
        # })
      }
    } else {
      W.new <- data.frame(W[rep(1:n, n), , drop=FALSE])
      A.new <- rep(A, each = n)
      # muhat.mat <- matrix(mu(a.new, x.new), byrow = TRUE, nrow = n, ncol = n) # BUG
      muhat.mat <- matrix(control$mu.hat(A.new, W.new), byrow = TRUE, nrow = n, ncol = n)
      muhat.obs <- diag(muhat.mat)
      mhat.obs <- rowMeans(muhat.mat)
      # mhat.mat <- matrix(mhat.obs, byrow = TRUE, ncol = n, nrow = n) # BUG
    }
  }

  control$mu.predicted <- muhat.obs
  control$m.means.predicted <- mhat.obs

  if('loclin' %in% control$method) {
    ghat.obs <- rep(NA, n)
    if(control$cross.fit) {
      for(v in 1:control$V) {
        ghat.obs[control$folds == v] <- control$g.hat[[v]]
      }
    } else {
      ghat.obs <- control$g.hat
    }
  }
  return(control)
}

cmdSuperLearner <- function(A, W, newA = A, newW = W, control = list(), cvControl = list()) {
  # num.libraries <- length(control$SL.library)
  n <- nrow(W)

  call <- match.call(expand.dots = TRUE)
  control <- do.call("cmdSuperLearner.control", control)
  cvControl <- do.call("cmdSuperLearner.CV.control", cvControl)

  validRows <- cmdCVFolds(n = n, cvControl = cvControl)

  fits <- NULL
  for(b in control$n.bins) {
    if(control$verbose) cat("\nEstimating models with", b, "bins... \n")
    fits[[paste0('dens.fit.', b, 'bins')]] <- cmdSuperLearner.onebin(A, W, newA=newA, newW = newW, b=b, SL.library = control$SL.library, verbose = control$verbose, validRows = validRows, saveFitLibrary = control$saveFitLibrary)
  }

  algs.per.bin <- ncol(fits[[1]]$cv.library.densities)
  n.algs <- length(control$n.bins) * algs.per.bin
  cv.library.densities <- matrix(NA, nrow=n, ncol=n.algs)
  library.densities <- matrix(NA, nrow=length(newA), ncol=n.algs)
  library.names <- NULL
  start.col <- 1
  for(b in control$n.bins) {
    end.col <- start.col + algs.per.bin - 1
    cv.library.densities[,start.col:end.col] <- fits[[paste0('dens.fit.', b, 'bins')]]$cv.library.densities
    library.densities[,start.col:end.col] <- fits[[paste0('dens.fit.', b, 'bins')]]$library.densities
    library.names <- c(library.names, fits[[paste0('dens.fit.', b, 'bins')]]$alg.names)
    start.col <- end.col + 1
  }

  if(control$verbose) cat("\nOptimizing model weights...\n")

  # Remove algs with errors in cv predictions
  errors.in.library <- apply(cv.library.densities, 2, function(col) any(is.na(col)))
  if(any(errors.in.library)) warning(paste0("Errors in the following candidate algorithms: ", library.names[which(errors.in.library)]))
  n.include <- sum(!errors.in.library)

  # Do SL log-likelihood optimization
  cv_risk <- function(beta) -mean(log(cv.library.densities[,!errors.in.library] %*% beta))
  capture.output(solnp_solution <- solnp(rep(1/n.include, n.include), cv_risk, eqfun=sum, eqB=1, ineqfun=function(beta) beta, ineqLB=rep(0,n.include), ineqUB=rep(1, n.include)))
  coef <- rep(0, n.algs)
  coef[!errors.in.library] <- solnp_solution$pars
  if(control$verbose) {
    num.libraries <- length(library.names)
    if (num.libraries >= 5){
      cat("Top 5 learners by weight: \n")
      for(j in 1:5) {
        cat(library.names[order(coef, decreasing = TRUE)[j]], " (weight ", sort(coef, decreasing = TRUE)[j], ")\n", sep='')
      }
    }else{
      cat(sprintf("Top %d learners by weight: \n", num.libraries))
      for(j in 1:num.libraries) {
        cat(library.names[order(coef, decreasing = TRUE)[j]], " (weight ", sort(coef, decreasing = TRUE)[j], ")\n", sep='')
      }
    }
  }
  SL.density <- c(library.densities[,!errors.in.library,drop=FALSE] %*% solnp_solution$pars)

  return(list(fits = fits, cv.library.densities = cv.library.densities, library.densities = library.densities, SL.densities = SL.density, coef = coef, library.names = library.names, a.ecdf = ecdf(A), control=control, cvControl = cvControl))
}

predict.cmdSuperLearner <- function(fit, newA, newW, threshold = .001) {
  newW <- as.data.frame(newW)
  new.U <- fit$a.ecdf(newA)
  trunc.coef <- fit$coef
  trunc.coef[trunc.coef < threshold] <- 0
  trunc.coef <- trunc.coef / sum(trunc.coef)
  nonzero <- which(trunc.coef > 0)
  lib.name.splits <- strsplit(fit$library.names, "_")
  lib.name.nbins <- unlist(lapply(lib.name.splits, function(l) as.numeric(strsplit(l[3], "bins")[[1]])))
  lib.name.alg <- unlist(lapply(lib.name.splits, function(l) paste0(l[1:2], collapse="_")))
  bins.to.fit <- unique(lib.name.nbins[nonzero])
  pred.densities <- matrix(NA, nrow=length(newA), ncol=length(fit$library.names))
  for(bin in bins.to.fit) {
    ind <- which(fit$control$n.bins == bin)
    if(length(ind) > 0){  ### Bug fixed
      new.bins <- .find.bin(new.U, bins = fit$fits[[ind]]$bins)
      bin.sizes <- unlist(lapply(fit$fits[[ind]]$bins, interval_measure))
      pred.probs <- matrix(NA, nrow = length(new.U), ncol = length(unique(lib.name.alg)))
      for(k in 1:length(fit$fits[[ind]]$bin.fits)) {
        if(any(new.bins == k)) {
          pred.probs[new.bins == k,] <- predict.SuperLearner(fit$fits[[ind]]$bin.fits[[k]], newdata = newW[new.bins == k,, drop=FALSE], onlySL = TRUE)$library.predict
        }
      }
      pred.probs <- pred.probs / rowSums(pred.probs)
      pred.densities[,which(lib.name.nbins == bin)] <- pred.probs / bin.sizes[new.bins]
    }
  }
  c(pred.densities[,nonzero,drop=FALSE] %*% trunc.coef[nonzero])
}

cmdCVFolds <- function (n, cvControl) {
  if (!is.null(cvControl$validRows)) return(cvControl$validRows)
  stratifyCV <- cvControl$stratifyCV
  shuffle <- cvControl$shuffle
  V <- cvControl$V
  if (shuffle) {
    validRows <- split(sample(1:n), rep(1:V, length = n))
  }
  else {
    validRows <- split(1:n, rep(1:V, length = n))
  }

  return(validRows)
}

.find.bin <- function(x, bins) {
  mat <- t(sapply(x-1e-10, function(x0) {
    unlist(lapply(bins, function(bin) {
      interval_contains_element(bin, x0)
    }))
  }))
  if(any(rowSums(mat) > 1)) stop("Overlapping bins")
  if(any(rowSums(mat) == 0)) stop("Element outside all bins")

  apply(mat, 1, function(row) which(row))
}
