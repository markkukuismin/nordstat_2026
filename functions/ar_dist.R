ar_dist = function(x, h = "silverman"){
  
  f0 = mvtnorm::dmvnorm(x)
  
  n = nrow(x)
  
  fhat = Rfast2::kernel(x, h = h)
  
  rhos = f0/fhat
  
  rhos[is.na(rhos)] = 0
  
  rhos[rhos > 1] = 1
  
  rho = mean(rhos)
  
  return(list(rho = rho, rhos = rhos))
  
}