

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

set.seed(1234)

ddist = "Multivariate logistic" # Multivariate logistic, uniform, normal, t

p = 500
sample_sizes = c(10, 20, 40, 60, 100)

ar_stat = rep(0, length(sample_sizes))

rhos = vector(mode = "list", length = length(sample_sizes))

j = 1

df = 5

I = diag(1, p)
mu = rep(0, p)

for(n in sample_sizes){
  
  if(ddist == "Multivariate logistic"){
    
    x = NonNorMvtDist::rmvlogis(n,
                                parm1 = rep(0, p),
                                parm2 = rep(1, p))
    
  }
  
  if(ddist == "Multivariate normal"){
    
    x = MASS::mvrnorm(n, mu = mu, Sigma = I)
    
  }
  
  if(ddist == "Multivariate t"){
    
    x = mvtnorm::rmvt(n = n, delta = mu, sigma = I, df = df)
    
  }
  
  if(ddist == "Multivariate uniform"){
    x = matrix(runif(n*p, min = -4, max = 4), n, p)
  }
  
  rho_res = ar_dist(x)
  
  ar_stat[j] = rho_res$rho
  
  rhos[[j]] = rho_res$rhos
  
  cat("\r", j)
  
  j = j + 1
  
}

U = L = c()

for(j in 1:length(sample_sizes)){
  
  UL = poibin::qpoibin(c(0.005, 0.995), pp = rhos[[j]])/sample_sizes[j]
  
  U = c(U, UL[1])
  L = c(L, UL[2])
  
  
}

library(ggplot2)

df = data.frame(rho = ar_stat,
                U = U,
                L = L,
                x = 1:length(sample_sizes))

pl = ggplot(df, aes(x = rho, y = x)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmax = U, xmin = L), width = 0.2, linewidth = 1) +
  xlim(c(0, 1)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 1.5) +
  ylab("Sample size") +
  scale_y_continuous(labels=as.character(sample_sizes)) +
  ggtitle(paste0(ddist, ", p = ", p))

pl

ggsave(paste0("simulations/figures/p", p, "_", ddist, ".png"))
