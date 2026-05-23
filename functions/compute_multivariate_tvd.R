
compute_multivariate_tvd <- function(p = 10) {
  
  B <- 5*10^5

  I =  diag(1, p)
  mu = rep(0, p)
  
  # Sample from a proposal distribution.
  # Here we use g = N(0, I_p) as the proposal.
  x <- mvtnorm::rmvnorm(B, mean = mu, sigma = I)
  
  g_vals <- mvtnorm::dmvnorm(x, mean = mu, sigma = I)
    
  f_vals = NonNorMvtDist::dmvunif(x)
  
  # Since x ~ g, we estimate:
  # 0.5 * int |f-g| dx
  # = 0.5 * E_g[ |f(X)-g(X)| / g(X) ]
  Y <- 0.5*abs(f_vals - g_vals)/g_vals
  tv_mc = mean(Y)
  
  tv_mc
  
}

