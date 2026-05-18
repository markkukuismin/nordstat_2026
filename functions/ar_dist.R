ar_dist = function(x, h = "silverman", sim_dist = FALSE, M = 10^3){
  
  f0 = mvtnorm::dmvnorm(x)
  
  n = nrow(x)
  
  fhat = Rfast2::kernel(x, h = h)
  
  rhos = f0/fhat
  
  rhos[is.na(rhos)] = 0
  
  rhox = rhos
  
  rhos[rhos > 1] = 1
  
  rho = mean(rhos)
  
  rho_sim = NULL
  
  if(sim_dist){
    
    rho_sim = rep(0, M)
    
    for(i in 1:M){
      
      u = runif(n)
      
      rho_sim[i] = mean(rhox > u)
      
    }
    
  }
  
  return(list(rho = rho, rhos = rhos, rho_sim = rho_sim))
  
}