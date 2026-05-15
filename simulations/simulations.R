# Test that do samples come from multivariate
# normal dist.

library(energy)
library(mvnormalTest)
library(mvtnorm)
library(mnormt)
library(Rlab)
library(gTests)
library(expm)
library(MASS)
#library(normwhn.test)
library(MVN)

source("functions/chen_xia_2021/funcs.r")
source("functions/chen_xia_2021/funcs2.r")

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

set.seed(1)

n = 50

# Power when alternative is multivariate t, logis and normal 

alt_dist = "mvlogis" # mvt, norm, mvlogis

p = 40 # 40, 50, 70

df = 5

# Method of Chen & Xia, 2021

BB = 50  #repeat times for Step 2 in Algo 1
L = 1  #estimate p-values L times

M = 1000

ar_stat = energy_stat = Tn_stat = rep(0, M)

NEW_pvalue = matrix(0, M, 2)

I = diag(1, p)
mu = rep(0, p)

rhos = rhos_null = matrix(0, n, M)

for(j in 1:M){
  
  if(alt_dist == "mvlogis"){
    
    x = NonNorMvtDist::rmvlogis(n,
                                parm1 = rep(0, p),
                                parm2 = rep(1, p))
    
  }
  
  if(alt_dist == "norm"){
    
    x = MASS::mvrnorm(n, mu = mu, Sigma = I)
    
  }
  
  if(alt_dist == "mvt"){
    
    x = mvtnorm::rmvt(n = n, delta = mu, sigma = I, df = df)
    
  }
  
  rho_res = ar_dist(x)
  
  ar_stat[j] = rho_res$rho
  
  rhos[, j] = rho_res$rhos
  
  energy_stat[j] = energy::mvnorm.test(x,
                                       R = 2)$statistic
  
  if(p < n) Tn_stat[j] = as.numeric(mvnormalTest::mvnTest(x, B = 10)$mv.test["p-value"])
  
  temp = getp(x, L=L, BB=BB)
  
  NEW_pvalue[j,] = temp$Yp
  
  cat("\r", j)
  
}

# Under null

ar_stat_null = Tn_stat_null = 
  energy_stat_null = rep(0, M)

for(j in 1:M){
  
  x = MASS::mvrnorm(n, 
                    mu = mu,
                    Sigma = I)
  
  energy_stat_null[j] = energy::mvnorm.test(x,
                                            R = 2)$statistic
  
  rho_res = ar_dist(x)
  
  ar_stat_null[j] = rho_res$rho
  
  rhos_null[, j] = rho_res$rhos
  
  #if(p < n) Tn_stat_null[j] = as.numeric(mvnormalTest::mvnTest(x, B = 2)$mv.test["Tn"])
  
  cat("\r", j)
  
}

qer = quantile(energy_stat_null, 0.95, na.rm = TRUE)
beta_energy = mean(energy_stat > qer, na.rm = TRUE)

#qTn = quantile(Tn_stat_null, 0.05)
beta_Tn = mean(Tn_stat < 0.05)

qar = quantile(ar_stat_null, 0.05)
beta_ar = mean(ar_stat < qar)

beta_NEW = mean(NEW_pvalue[, 1] < 0.05)

b = c(round(beta_ar, 2),
      round(beta_NEW, 2),
      round(beta_energy, 2),
      round(beta_Tn, 2))

names(b) = c("AR", "NEW", "E", "Tn")
b
paste0(b, collapse = " & ")

hist(ar_stat_null, 
     probability = T,
     xlim = c(min(c(ar_stat, ar_stat_null)),
              max(c(ar_stat, ar_stat_null))))

hist(ar_stat, 
     probability = T,
     col = "lightblue",
     add = T)

legend("topleft",
       col = c("gray", "lightblue"),
       pch = 15,
       cex = 1.5,
       legend = c("Null", "Alt"))

D1 = data.frame(method = names(b),
                pwr = b,
                dist = alt_dist,
                n = n,
                p = p)

D2 = data.frame(rhos = c(rhos),
                rhos_null = c(rhos_null),
                n = n,
                p = p,
                sim = rep(1:M, each = n))

fpath1 = "simulations/demo_res.txt"
fpath2 = "simulations/demo_res_extra.txt"

write.table(D1, 
            file = fpath1,
            append = TRUE,
            row.names = FALSE)

write.table(D2, 
            file = fpath2,
            append = TRUE,
            row.names = FALSE)

x = matrix(rnorm(1000*6), 1000, 6)
system.time(d <- ar_dist(x, thumb = "silverman"))

x = matrix(runif(1000*6), 1000, 6)
system.time(d <- ar_dist(x, thumb = "silverman"))
