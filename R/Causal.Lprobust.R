# Debiased kernel regression with confidence intervals.
# Based on Calonico et al "on the effect of bias estimation on coverage 
# accuracy in nonparametric inference".
# Originally implemented as "nprobust" by Calonico et al, 
# modified by Kenta Takatsu

# a collection of kernel functions
W.fun = function(u, kernel){
  if (kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w = 0.5*(abs(u)<=1)
  if (kernel=="tri") w = (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w = dnorm(u)
  return(w)
}

# compute (X'X)^-1 
qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}

# compute the estimate of conditional residual 
lprobust.res = function(eX, eY, nnmatch, edups, edupsid) {
  n = length(eY)
  res = matrix(NA,n,1)
  for (pos in 1:n) {
    rpos = edups[pos] - edupsid[pos] 
    lpos = edupsid[pos] - 1 
    while (lpos+rpos < min(c(nnmatch,n-1))) {
      if (pos-lpos-1 <= 0) rpos = rpos + edups[pos+rpos+1]
      else if (pos+rpos+1>n) lpos = lpos + edups[pos-lpos-1]
      else if ((eX[pos]-eX[pos-lpos-1]) > (eX[pos+rpos+1]-eX[pos])) rpos = rpos + edups[pos+rpos+1]
      else if ((eX[pos]-eX[pos-lpos-1]) < (eX[pos+rpos+1]-eX[pos])) lpos = lpos + edups[pos-lpos-1]
      else {
        rpos = rpos + edups[pos+rpos+1]
        lpos = lpos + edups[pos-lpos-1]
      }
    }
    ind.J = (pos-lpos):min(c(n,(pos+rpos)))
    y.J   = sum(eY[ind.J])-eY[pos] # sum of neighbors 
    Ji = length(ind.J)-1 # number of neighbors
    res[pos,1] = sqrt(Ji/(Ji+1))*(eY[pos] - y.J/Ji) # residual
  }
  return(res)
}

# simplified version of nprobust library
# returns explicit bias estimation 
my.lprobust <- function(x, y, h, b, eval=NULL, kernel.type="epa"){
  # sort input
  ind <- order(x)
  x <- x[ind]; y <- y[ind]
  if (is.null(eval)){
    x.max <- max(x); x.min <- min(x); n <- length(x)
    eval <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval)
  
  deriv <- 0; p <- 1; q <- 3 # estimate first moment (p) with bias correction (q)
  if (length(h) != neval){
    h <- rep(h,neval); b <- rep(b,neval); rho <- h/b
  }
  N <- length(x); dups <- dupsid <- 0;
  for (j in 1:N) {
    dups[j]=sum(x==x[j]) # counting duplicates
  }
  j=1
  while (j<=N) {
    dupsid[j:(j+dups[j]-1)] <- 1:dups[j]
    j <- j+dups[j]
  }
  Estimate=matrix(NA,neval,8)
  colnames(Estimate)=c("eval", "h", "b", "N", "mu.hat", "b.hat", "se.us", "se.rb")
  
  # for each evaluation point
  for (i in 1:neval){
    bw.min   <- sort(abs(x-eval[i]))[31] # this prevents non-invertible X'T X
    h[i]     <- max(h[i], bw.min)
    b[i]     <- max(b[i], bw.min)
    
    w.h   <- W.fun((x-eval[i])/h[i], kernel.type)/h[i] # apply kernel.type with bandwith h
    w.b   <- W.fun((x-eval[i])/b[i], kernel.type)/b[i] # apply kernel.type with bandwith b
    ind.h <- w.h>0;  ind.b <- w.b>0               # select non zero entry
    N.h   <- sum(ind.h);  N.b <- sum(ind.b)
    
    ind   <- ind.b
    if (h[i]>b[i]) ind <- ind.h  # choose index for wider bandwith (more index). same if rho = 1
    eN  <- sum(ind) # number of all points considered
    eY  <- y[ind] # original Y
    eX  <- x[ind] # original X
    W.h <- w.h[ind] # kernelized X with h
    W.b <- w.b[ind] # kernelized X with b
    
    edups   <- dups[ind]; edupsid <- dupsid[ind]
    u   <- (eX-eval[i])/h[i]
    R.q <- matrix(NA,eN,(q+1)) # (1, u, u^2) for p = 1
    for (j in 1:(q+1))  R.q[,j] <- (eX-eval[i])^(j-1)
    R.p <- R.q[,1:(p+1)] # (1, u) for p = 1
    
    invG.q  <- qrXXinv((sqrt(W.b)*R.q)) # (R.q^T W.b R.q)^-1
    invG.p  <- qrXXinv((sqrt(W.h)*R.p)) # (R.p^T W.h R.p)^-1
    e.p1    <- matrix(0,(q+1),1); e.p1[p+2]=1  # (0,0,1)
    e.v     <- matrix(0,(p+1),1); e.v[deriv+1]=1 # (0,1)
    L <- crossprod(R.p*W.h,u^(p+1)) # Lambda.p
    
    Q.q.bias <- t(h[i]^(p+1)*(L%*%t(e.p1))%*%t(t(invG.q%*%t(R.q))*W.b)) # h^2  
    beta.p  <- invG.p%*%crossprod(R.p*W.h, eY) #R.p^T W.h Y
    beta.bias <- invG.p%*%crossprod(Q.q.bias, eY) 
    tau <- beta.p[(deriv+1),1]
    tau.bias <- beta.bias[(deriv+1),1]
    
    # confidence interval 
    nnmatch <- 3
    # nonparametric estimate of conditional residual via KNN
    res.h <- lprobust.res(eX, eY, nnmatch, edups, edupsid); res.b <- res.h
    V.Y.cl <- invG.p%*%crossprod(c(res.h)*as.matrix(R.p*W.h))%*%invG.p
    V.Y.bc <- invG.p%*%crossprod(c(res.h)*as.matrix(R.p*W.h-Q.q.bias))%*%invG.p
    se.cl  <- sqrt(factorial(deriv)^2*V.Y.cl[deriv+1,deriv+1])
    se.rb  <- sqrt(factorial(deriv)^2*V.Y.bc[deriv+1,deriv+1])
    
    # evaluation points, bandwidth h, bandwidth b, data instances used after kernel, 
    # estimated mean, estimated bias, original CI, bias-corrected CI
    Estimate[i,] <- c(eval[i], h[i], b[i], eN, tau, tau.bias, se.cl, se.rb) 
  }
  Estimate
}

my.locpoly <- function(x, y, h, eval=NULL, kernel.type="epa"){
  # sort input
  ind <- order(x); x <- x[ind]; y <- y[ind]
  # by default use 30 equidistant evaluation points 
  if (is.null(eval)){
    x.max <- max(x); x.min <- min(x); n <- length(x)
    eval <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval)
  
  if (length(h) != neval){
    h <- rep(h,neval);
  }
  Estimate=matrix(NA, neval, 3)
  colnames(Estimate)=c("eval", "h", "mu.hat")
  
  # for each evaluation point
  for (i in 1:neval){
    bw.min <- sort(abs(x-eval[i]))[31] # this prevents non-invertible X'T X
    h[i] <- max(h[i], bw.min)

    w.h   <- W.fun((x-eval[i])/h[i], kernel.type)/h[i] # apply kernel.type with bandwith h

    ind <- w.h>0  # choose index for wider bandwith (more index). same if rho = 1
    eY  <- y[ind] # original Y
    eX  <- x[ind] # original X
    W.h <- w.h[ind] # kernelized X with h

    R.p <- matrix(NA,sum(ind),2) # (1, u) for p = 1
    for (j in 1:2)  R.p[,j] <- (eX-eval[i])^(j-1)
    
    invG.p  <- qrXXinv((sqrt(W.h)*R.p)) # (R.p^T W.h R.p)^-1
    beta.p  <- invG.p%*%crossprod(R.p*W.h, eY) #R.p^T W.h Y
    tau <- beta.p[1,1]
   
    Estimate[i,] <- c(eval[i], h[i], tau) 
  }
  Estimate
}

my.loclinear <- function(x, y, h, eval=NULL, kernel.type="epa"){
  # sort input
  ind <- order(x); x <- x[ind]; y <- y[ind]
  # by default use 30 equidistant evaluation points 
  if (is.null(eval)){
    x.max <- max(x); x.min <- min(x); n <- length(x)
    eval <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval)
  
  if (length(h) != neval){
    h <- rep(h,neval);
  }
  Estimate=matrix(NA, neval, 3)
  colnames(Estimate)=c("eval", "h", "mu.hat")
  
  # for each evaluation point
  for (i in 1:neval){
    bw.min <- sort(abs(x-eval[i]))[31] # this prevents non-invertible X'T X
    h[i] <- max(h[i], bw.min)
    
    w.h   <- W.fun((x-eval[i])/h[i], kernel.type)/h[i] # apply kernel.type with bandwith h
    
    ind <- w.h>0  # choose index for wider bandwith (more index). same if rho = 1
    eY  <- y[ind] # original Y
    eX  <- x[ind] # original X
    W.h <- w.h[ind] # kernelized X with h
    
    R.p <- matrix(NA,sum(ind),2) # (1, u) for p = 1
    for (j in 1:2)  R.p[,j] <- (eX-eval[i])^(j-1)
    
    invG.p  <- qrXXinv((sqrt(W.h)*R.p)) # (R.p^T W.h R.p)^-1
    beta.p  <- invG.p%*%crossprod(R.p*W.h, eY) #R.p^T W.h Y
    tau <- beta.p[1,1]
    
    Estimate[i,] <- c(eval[i], h[i], tau) 
  }
  Estimate
}
