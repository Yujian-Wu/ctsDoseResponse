# Bias-corrected inference of continuous treatment effect
# The implementation is based "npcausal" and "nprobust"
# modified by Kenta Takatsu
# @param y outcome of interest.
# @param a continuous treatment.
# @param x covariate matrix.
# @param bw.seq sequence of bandwidth values.
# @param eval.pts sequence of evaluation points.
# @param mu funciton for estimated outcome regession Pn(y | a, x)
# @param g funciton for estimated conditional density Pn(a | x)
# @param tau the ratio parameter (h/b) between bandwidths 
#   for mean (h) and estimated second derivative (b). default=1
debiased.ctseff <- function (y, a, x, bw.seq, eval.pts, mu, g, 
                             tau=1, 
                             kernel.type="epa", 
                             verbose=FALSE) {
  # compute point estimates and point-wise CIs of ATE 
  
  x.names <- colnames(x)
  
  # process data
  n <- length(a); ord <- order(a)
  Y <- y[ord]; X <- x[ord,]; A <- a[ord] # sort by the order of A
  a.max <- max(A); a.min <- min(A); n <- length(A)
  
  X <- data.frame(X); colnames(X) <- x.names
  X.new <- data.frame(X[rep(1:n, n), ])  # repeat a matrix over rows 
  colnames(X.new) <- x.names
  A.new <- rep(A, each = n) # repeat 
  
  # compute pseudo-outcome adjusted for confounding
  if(verbose) cat("\nComputing pseudo outcome\n")
  muhat.mat <- matrix(mu(A.new, X.new), byrow = TRUE, nrow = n, ncol = n)
  muhat.obs <- diag(muhat.mat)    
  mhat.obs <- rowMeans(muhat.mat) 
  ghat.obs <- g(A, X)
  # cat("\nMu is", muhat.obs)
  # cat("\nM is", mhat.obs)
  # cat("\nG is", ghat.obs)
  pseudo.out <- (Y - muhat.obs)/(ghat.obs) + mhat.obs

  if(verbose) cat("\nComputing bandwidth\n")
  # Data-driven bandwidth selection (section 3.5 of Kennedy)
  kern <- function(t) {
    if(kernel.type=="gau"){return (dnorm(t))}
    if(kernel.type=="epa"){return (0.75*(1-t^2)*(abs(t)<=1))}
    if(kernel.type=="uni"){return (0.5*(abs(t)<=1))}
    if(kernel.type=="tri"){return ((1-abs(t))*(abs(t)<=1))}
  }
  if(verbose) cat("\nRun crossvalidation for bandwidth selection\n")
  
  # NOTE: There are many methods for data-driven bandwidth selection
  # Line 56 is based on Calonico et al 
  # Line 57 - 80 is based on leave-one-out CV
  h <- lpbwselect(pseudo.out,A, eval=eval.pts, bwselect="imse-dpi")$bws[,2]
  h <- rep(h, length(eval.pts)); b <- h/tau
  est <- my.lprobust(A, pseudo.out, h, b, eval = eval.pts, 
                     kernel.type=kernel.type) 
  
  if(verbose) cat("\nEstimate pointwise confidence interval with EIF\n")
  # estimate variance with an influence function-based method
  # See influence_functions.R 
  inf.fns <- mapply(function(a0, h.val, b.val){
    compute.rinfl.func(pseudo.out, A, a0, h.val, b.val, kern, 
                       muhat.mat, mhat.obs)
  }, eval.pts, h, b)
  se.robust <- apply(inf.fns, 2, sd)/sqrt(n)
  
  # build output
  # debiased point estimate is mu-b
  # 95% CI based on an IF is mu-b +/- 1.96*se.infl.robust
  # se.rb is a robust se based on Calonico et al
  out <- data.frame(x=eval.pts, # evalation points on X
                    mu=est[,"mu.hat"],  # the estimate of mean 
                    b=est[,"b.hat"],    # the estimate of bias
                    se.infl.robust=se.robust, # se based on an influence function
                    se.rb= est[,"se.rb"]) # se based on fixed-n 
  return(out)
}

