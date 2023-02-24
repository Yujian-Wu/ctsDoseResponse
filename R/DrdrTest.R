drdrtest.base <- function(y,a,pi,varpi,mu,ma,arange,h=NULL,b=1000,dist='TwoPoint',a.grid.size=401){
  pseudo.out <- (y-mu)/(pi/varpi)+ma
  n <- length(y)
  if(is.null(h)){
    h <- 0.9*n^(-1/5)*min(stats::sd(a),stats::IQR(a)/1.34)
  }
  test.out <- wildboot(a,pseudo.out,h,b,arange,dist,a.grid.size)
  return(list(p.value = mean(test.out$T.boot>=test.out$T.obs),
              test.stat = test.out$T.obs,
              Bootstrap.samples = test.out$T.boot,
              loc.fit = test.out$loc.fit,
              bandwidth=h))
}

drdrtest <- function(y,a,l,arange,pifunc,mufunc,h=NULL,b=1000,dist='TwoPoint',pi.low=0.01,a.grid.size=401){
  approx.fn <- function(x,y,z){
    stats::predict(stats::smooth.spline(x,y),x=z)$y
  }
  
  
  n <- length(y)
  la <- data.frame(l,a=a)
  a.vals <-seq(arange[1],arange[2],length.out=a.grid.size)
  
  rep.l <- data.frame(l[rep(1:n,length(a.vals)),],a=rep(a.vals,each=n))
  colnames(rep.l)[1:length(ncol(l))] <- colnames(l)
  la.new <- rbind(la, rep.l)
  
  l.new <- la.new[,-ncol(la.new)]
  a.new <- la.new[,ncol(la.new)]
  
  l.new <- as.data.frame(l.new)
  colnames(l.new) <- colnames(l)
  
  pihat.vals <- pifunc(a.new, l.new)
  pihat <- pmax(pihat.vals[1:n],pi.low)
  pihat.mat <- matrix(pihat.vals[-(1:n)],nrow = n , ncol= length(a.vals))
  varpihat <- pmax(approx.fn(a.vals,apply(pihat.mat,2,mean),a),pi.low)
  
  muhat.vals <- mufunc(a.new,l.new)
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)],nrow=n,ncol=length(a.vals))
  mhat <- approx.fn(a.vals,apply(muhat.mat,2,mean),a)
  
  return(drdrtest.base(y,a,pihat,varpihat,muhat,mhat,arange,h,b,dist,a.grid.size))
}

rrademachar <- function(n){
  ## intput:
  ## n: number of observations
  ##
  ## output:
  ## a vector containing the simulated observations        
  x <- c(1+sqrt(5),1-sqrt(5))/2
  p <- (sqrt(5)-1)/(2*sqrt(5))
  return(sample(x,size=n,replace = TRUE,prob=c(p,1-p)))
}

rtwopoint <- function(n){
  ## intput:
  ## n: number of observations
  ##
  ## output:
  ## a vector containing the simulated observations        
  x <- c(1+sqrt(5),1-sqrt(5))/2
  p <- (sqrt(5)-1)/(2*sqrt(5))
  return(sample(x,size=n,replace = TRUE,prob=c(p,1-p)))
}

wildboot <- function(x,y,h,b,xrange,dist='TwoPoint',x.grid.size){
  n <- length(y)
  yhat <- mean(y)
  res <- y-yhat
  if(dist=='TwoPoint'){
    estar <- matrix(rtwopoint(length(y)*b),nrow=n)*res
  }else if(dist=='Rademachar'){
    estar <- matrix(rrademachar(length(y)*b),nrow=n)*res
  }else{
    stop("dist should be either 'TwoPoint' or 'Rademachar'")
  }
  yboot <- yhat+estar
  tboot.list <- rep(0,length.out=b)
  for(i in 1:b){
    loc.wild <- KernSmooth::locpoly(x,yboot[,i],drv=0,bandwidth = h, range.x=xrange,gridsize = x.grid.size)
    linear.est <- mean(yboot[,i])
    Tboot <- sum((loc.wild$y-linear.est)^2)*(loc.wild$x[2]-loc.wild$x[1])*n*sqrt(h)
    tboot.list[i] <- Tboot
  }
  loc.obs <- KernSmooth::locpoly(x,y,drv=0,bandwidth =h, range.x=xrange,gridsize = x.grid.size)
  Tobs <- sum((loc.obs$y-yhat)^2)*(loc.obs$x[2]-loc.obs$x[1])*n*sqrt(h)
  
  return(list(T.obs = Tobs, T.boot = tboot.list,loc.fit = loc.obs))
}