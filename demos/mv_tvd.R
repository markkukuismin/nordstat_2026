
library(ggplot2)
library(mvtnorm)

source("functions/ar_dist.R")

set.seed(1)

B <- 10^6

p = 50

du = 5

ddist = "Multivariate logistic"

I =  diag(1, p)
mu = rep(0, p)

# Sample from a proposal distribution.
# Here we use g = N(0, I_p) as the proposal.
x <- rmvnorm(B, mean = mu, sigma = I)

g_vals <- dmvnorm(x, mean = mu, sigma = I)
if(ddist == "t") f_vals <- dmvt(x, delta = mu, sigma = I, df = du, log = FALSE)
if(ddist == "logis") f_vals = NonNorMvtDist::dmvlogis(x,parm1 = mu,parm2 = rep(1, p))


if(ddist == "norm"){
  
  f_vals = dmvnorm(x, mean = mu, sigma = I)
  
}

# Since x ~ g, we estimate:
# 0.5 * int |f-g| dx
# = 0.5 * E_g[ |f(X)-g(X)| / g(X) ]
Y <- 0.5*abs(f_vals - g_vals)/g_vals
tv_mc = mean(Y)

1 - tv_mc
mc_se <- sd(Y)/sqrt(B)

tv_mc
mc_se

##

sample_sizes = c(20, 40, 60, 80, 100)

ar_rhos = vector(mode = "list", length = length(sample_sizes))

j = 1

ar_stat = c()

for(n in sample_sizes){
  
  if(ddist == "Multivariate logistic"){
    
    x = NonNorMvtDist::rmvlogis(n,
                                parm1 = mu,
                                parm2 = rep(1, p))
    
  }
  
  if(ddist == "Multivariate normal"){
    
    x = MASS::mvrnorm(n, mu = mu, Sigma = I)
    
  }
  
  if(ddist == "Multivariate t"){
    
    x = mvtnorm::rmvt(n = n, delta = mu, sigma = I, df = du)
    
  }
  
  rho_res = ar_dist(x)
  
  ##
  
  ar_stat[j] = rho_res$rho
  
  ar_rhos[[j]] = rho_res$rhos
  
  j = j + 1
    
}

LU = matrix(0, length(sample_sizes), 2)

for(i in 1:length(sample_sizes)){
  
  LU[i, ] = poibin::qpoibin(c(0.025, 0.975), pp = ar_rhos[[i]])/sample_sizes[i]
  
}


df = data.frame(rho = ar_stat,
                L = LU[, 1],
                U = LU[, 2],
                n = sample_sizes)

##

M = 300

j = 1

for(k in 1:(M - 1)){
    
    for(n in sample_sizes){
      
      if(ddist == "Multivariate logistic"){
        
        x = NonNorMvtDist::rmvlogis(n,
                                    parm1 = mu,
                                    parm2 = rep(1, p))
        
      }
      
      if(ddist == "Multivariate normal"){
        
        x = MASS::mvrnorm(n, mu = mu, Sigma = I)
        
      }
      
      if(ddist == "Multivariate t"){
        
        x = mvtnorm::rmvt(n = n, delta = mu, sigma = I, df = du)
        
      }
      
      rho_res = ar_dist(x)
      
      ##
      
      ar_stat[j] = rho_res$rho
      
      ar_rhos[[j]] = rho_res$rhos
      
      j = j + 1
      
    }
  
  j = 1
  
  LU = matrix(0, length(sample_sizes), 2)
  
  for(i in 1:length(sample_sizes)){
    
    LU[i, ] = poibin::qpoibin(c(0.025, 0.975), pp = ar_rhos[[i]])/sample_sizes[i]
    
  }
  
  df2 = data.frame(rho = ar_stat,
                  L = LU[, 1],
                  U = LU[, 2],
                  n = sample_sizes)
  
  df = df + df2
  
}

df = df/M

pl = ggplot(df, aes(x = n, y = rho)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = L, ymax = U), width = 2, linewidth = 1) +
  ylim(c(0, 1)) +
  geom_hline(yintercept = tv_mc + c(-1, 1)*1.96*mc_se, linetype = 2, linewidth = 1.5) +
  geom_hline(yintercept = tv_mc, linewidth = 1.5) +
  xlab("Sample size") +
  ggtitle(paste0(ddist, ", p = ", p))

pl
