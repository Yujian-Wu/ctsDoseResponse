# compute influence function for bias-corrected estimator with fixed h and b
# res: residuals y - E[y|x]
# A: observed exposure
# a: the exposure for which we want to evalute an IF
# h: bandwidth for local linear
# b: bandwidth for local polynomial 
# kern: kernel function
compute.rinfl.func <- function(Y, A, a, h, b, tau, kern, muhat.mat, mhat.obs){
  n <- length(A)
  bw.min <- sort(abs(A-a))[21] # this prevents non-invertible X'X (is there a better way?)
  h <- max(h, bw.min); b <- max(b, bw.min)
  a.std.h <- (A - a)/h; kern.std.h <- kern(a.std.h)/h
  a.std.b <- (A - a)/b; kern.std.b <- kern(a.std.b)/b
  
  c0.h <- mean(kern.std.h)
  c1.h <- mean(kern.std.h * a.std.h)
  c2.h <- mean(kern.std.h * a.std.h^2)
  
  c0.b <- mean(kern.std.b)
  c1.b <- mean(kern.std.b * a.std.b)
  c2.b <- mean(kern.std.b * a.std.b^2)
  c3.b <- mean(kern.std.b * a.std.b^3)
  c4.b <- mean(kern.std.b * a.std.b^4)
  c5.b <- mean(kern.std.b * a.std.b^5)
  c6.b <- mean(kern.std.b * a.std.b^6)
  
  Dh <- matrix(c(c0.h, c1.h, 
                 c1.h, c2.h), nrow = 2)
  Db <- matrix(c(c0.b, c1.b, c2.b, c3.b,
                 c1.b, c2.b, c3.b, c4.b,
                 c2.b, c3.b, c4.b, c5.b,
                 c3.b, c4.b, c5.b, c6.b), nrow = 4)
  
  g2.h <- (A - a)/h
  g2.b <- (A - a)/b
  g3.b <- ((A - a)/b)^2
  g4.b <- ((A - a)/b)^3
  
  int1.h <- colMeans(kern.std.h * (muhat.mat - mhat.obs))
  int2.h <- colMeans(g2.h * kern.std.h * (muhat.mat - mhat.obs))
  
  int1.b <- colMeans(kern.std.b * (muhat.mat - mhat.obs))
  int2.b <- colMeans(g2.b * kern.std.b * (muhat.mat - mhat.obs))
  int3.b <- colMeans(g3.b * kern.std.b * (muhat.mat - mhat.obs))
  int4.b <- colMeans(g4.b * kern.std.b * (muhat.mat - mhat.obs))
  
  gamma.h <- coef(lm(Y ~ a.std.h, weights = kern.std.h))
  model.b <- lm(Y ~ poly(a.std.b, 3), weights = kern.std.b)
  res.h <- Y - (gamma.h[1] + gamma.h[2]* a.std.h)
  res.b <- Y - predict(model.b)
  
  inf.fn <- t(solve(Dh) %*% rbind(res.h * kern.std.h + int1.h,
                                  g2.h * res.h * kern.std.h + int2.h))
  inf.fn.robust <- t(solve(Db) %*% rbind(res.b * kern.std.b+ int1.b,
                                         g2.b * res.b * kern.std.b + int2.b,
                                         g3.b * res.b * kern.std.b + int3.b,
                                         g4.b * res.b * kern.std.b + int4.b))
  # c2 is the second-moment of Kernel
  # you can empirically estimate this by 
  # c2.mat <- crossprod(cbind(rep(1, n), a.std)*kern.std, (a.std)^2)
  # e3 <- matrix(0, 4,1); e3[3] <- 1
  # c2.hat <- (solve(Dh)%*% c2.mat %*%t(e3) /n)[1,3]
  # I notice this is more stable for empirical results.
  
  c2 <- integrate(function(u){u^2*kern(u)}, -Inf,Inf)$value
  return (inf.fn[,1]-tau^2*c2*inf.fn.robust[,3])
}
