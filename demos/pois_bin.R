set.seed(1)

source("functions/compute_multivariate_tvd.R")

p = 5

Sigma = diag(1, p)
mu = rep(0, p)

N = c(100, 120, 150)

M = 10^5

i = 1

Q2 = data.frame()

par(mfrow = c(1, 3), 
    cex.axis=1.5,
    mar = c(5.1, 5.1, 4.1, 2.1))

for(n in N){
  
  x = NonNorMvtDist::rmvunif(n = n,
                             dim = p)
  
  fhat = Rfast2::kernel(x)
  
  rhox = mvtnorm::dmvnorm(x)/fhat
  
  rhox[is.na(rhox)] = 0
  
  rhox = ifelse(rhox > 1, 1, rhox)
  
  rho = mean(rhox)
  
  U = matrix(runif(n*M), n, M)
  
  f = function(x) rhox > x
  
  I = apply(U, 2, f)
  
  d = colMeans(I)
  
  a = table(d)/M
  
  a = as.numeric(names(a))
  
  plot(table(d)/M, 
       type = "h",
       xlab = expression(rho),
       #ylab = "Probability",
       ylab = " ",
       main = paste0("n = ", n),
       lwd = 4,
       cex.lab = 2,
       col = "darkgray",
       xaxt = "n",
       xlim = c(0, 0.15))
  
  axis(1,
       at = a,
       labels = as.character(round(a, 2)))
  
  xx = 0:n
  
  probpb = poibin::dpoibin(xx, pp = rhox)
  
  lines(xx/n + 0.002, 
        probpb, 
        type = "h", 
        lwd = 4, 
        col = rgb(red=0, green=0, blue=1, alpha=0.5))
  
  abline(v = rho, lty = 2, lwd = 4, col = "red")
  
  if(i == 1){
    tvd_u = compute_multivariate_tvd(p=p) 
  }
  
  abline(v = 1 - tvd_u, lwd = 3)
  
  i = i + 1
  
}

tvd_u
1 - tvd_u
