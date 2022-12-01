# Debiased kernel regression
# Originally implemented as "nprobust" by Calonico et al, 
# modified by Kenta Takatsu

kern = function(u, kernel="epa"){
  if (kernel=="epa") w <- 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w <- 0.5*(abs(u)<=1)
  if (kernel=="tri") w <- (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w <- dnorm(u)
  return(w)
}

qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}

my.lprobust <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)
  
  Estimate <- matrix(NA, neval, 4)
  colnames(Estimate) <- c("eval", "theta.hat", "mu.hat", "b.hat")
  c.2 <- integrate(function(t){t^2*kern(t, kernel.type)}, -Inf, Inf)$value
  
  for (i in 1:neval){
    # prevents noninvertible D matrix 
    bw.min   <- sort(abs(x-eval.pt[i]))[10] 
    h0     <- max(h, bw.min)
    b0     <- max(b, bw.min)
    
    k.h   <- kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- kern((x-eval.pt[i])/b0, kernel.type)/b0
    
    ind.h <- k.h>0;  ind.b <- k.b>0               
    ind   <- ind.b
    if (h>b) ind <- ind.h  # choose index for wider bandwidth 
    eN  <- sum(ind); eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind] 
    
    W.b <- matrix(NA, eN, 4) # (1, u, u^2, u^3)
    for (j in 1:4)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2) # (1, u) 
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    
    invD.b  <- qrXXinv((sqrt(K.b)*W.b)) 
    invD.h  <- qrXXinv((sqrt(K.h)*W.h)) 
    
    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y
    beta.bias <- (h/b)^2*invD.b%*%crossprod(W.b*K.b, eY)*c.2
    tau <- beta.ll[1,1]
    tau.bias <- beta.bias[3,1]
    
    Estimate[i,] <- c(eval.pt[i], tau-tau.bias, tau, tau.bias) 
  }
  Estimate
}

my.hatmatrix <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x); n <- length(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)
  
  c.2 <- integrate(function(t){t^2*kern(t, kernel.type)}, -Inf, Inf)$value
  e1 <- matrix(c(1, 0), ncol=2)
  e3 <- matrix(c(1, 0, 0, 0), ncol=4)
  w.h.zero   <- kern(0, kernel.type)/h  
  w.b.zero   <- kern(0, kernel.type)/b
  
  hat.mat <- rep(0, neval)
  for (i in 1:neval){
    # prevents noninvertible D matrix 
    bw.min   <- sort(abs(x-eval.pt[i]))[10] 
    h0     <- max(h, bw.min)
    b0     <- max(b, bw.min)
    
    k.h   <- kern((x-eval.pt[i])/h0, kernel.type)/h0 
    k.b   <- kern((x-eval.pt[i])/b0, kernel.type)/b0 
    ind.h <- k.h>0;  ind.b <- k.b>0               
    N.h   <- sum(ind.h);  N.b <- sum(ind.b)
    
    ind   <- ind.b
    if (h>b) ind <- ind.h  
    eY  <- y[ind]; eX  <- x[ind] 
    K.h <- k.h[ind]; K.b <- k.b[ind] 
    W.b <- matrix(NA,sum(ind),4) # (1, u, u^2, u^3) 
    for (j in 1:4)  W.b[,j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h  <- matrix(NA,sum(ind),2) # (1, u) 
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    
    invD.b  <- qrXXinv((sqrt(K.b)*W.b)) # (R.q^T W.b R.q)^-1
    invD.h  <- qrXXinv((sqrt(K.h)*W.h)) # (R.p^T W.h R.p)^-1
    hat.mat[i] <- (invD.h%*%t(e1*w.h.zero))[1,] - 
      ((h/b)^2*c.2*invD.b%*%t(e3*w.b.zero))[3,]
  }
  approx(eval.pt, hat.mat, xout = x)$y
}

robust.loocv <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x)
  x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x); n <- length(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  hat.val <- my.hatmatrix(x, y, h, b, eval.pt=eval.pt, kernel.type=kernel.type)
  est <- my.lprobust(x, y, h, b, eval.pt=eval.pt, kernel.type=kernel.type)
  est.fn <- approx(eval.pt, est[,"theta.hat"], xout = x)$y
  mean(((y - est.fn) / (1 - hat.val))^2, na.rm=TRUE)
}

my.locpoly <- function(x, y, h, eval=NULL, kernel.type="epa", degree=1){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval)){
    x.max <- max(x); x.min <- min(x)
    eval <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval)
  Estimate <- matrix(NA, neval, 2)
  colnames(Estimate) <- c("eval", "theta.hat")
  
  for (i in 1:neval){
    # prevents noninvertible D matrix 
    bw.min   <- sort(abs(x-eval.pt[i]))[10] 
    h0     <- max(h, bw.min)
    
    k.h   <- kern((x-eval[i])/h0, kernel.type)/h0
    ind <- k.h>0; eY  <- y[ind]; eX  <- x[ind] 
    K.h <- k.h[ind] 
    W.h <- matrix(NA,sum(ind), degree+1) 
    for (j in 1:degree+1)  W.h[,j] <- ((eX-eval[i])/h0)^(j-1)
    
    invD.h  <- qrXXinv((sqrt(K.h)*W.h)) 
    beta.ll  <- invD.h%*%crossprod(K.h*W.h, eY) 
    Estimate[i,] <- c(eval[i], beta.ll[1,1]) 
  }
  Estimate
}

my.loclinear <- function(x, y, h, eval=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval)){
    x.max <- max(x); x.min <- min(x)
    eval <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval)
  Estimate <- matrix(NA, neval, 2)
  colnames(Estimate) <- c("eval","theta.hat")
  
  # for each evaluation point
  for (i in 1:neval){
    # prevents noninvertible D matrix 
    bw.min   <- sort(abs(x-eval[i]))[10] 
    h0     <- max(h, bw.min)
    w.h   <- kern((x-eval[i])/h0, kernel=kernel.type)/h0 
    
    ind <- w.h>0  # choose index for wider bandwith (more index). same if rho = 1
    eY  <- y[ind]; eX  <- x[ind] 
    W.h <- w.h[ind] # kernelized X with h
    
    R.p <- matrix(NA, sum(ind), 2) # (1, u) for p = 1
    for (j in 1:2)  R.p[,j] <- ((eX-eval[i])/h0)^(j-1)
    
    invG.p  <- qrXXinv((sqrt(W.h)*R.p)) # (R.p^T W.h R.p)^-1
    beta.p  <- invG.p%*%crossprod(R.p*W.h, eY) #R.p^T W.h Y
    tau <- beta.p[1,1]
    Estimate[i,] <- c(eval[i], tau) 
  }
  Estimate
}