
library(ggplot2)
library(tidyverse)

source("functions/read.libsvm.R")
source("functions/ar_dist.R")

set.seed(1)

data_raw = read.libsvm("data/heart.txt",
                       dimensionality = 13)

# target: Presence/absence of heart disease
# age:	Age (years)
# sex:	Sex (1 = male, 0 = female)
# cp:	Chest pain type
# trestbps:	Resting blood pressure (mm Hg)
# chol:	Serum cholesterol (mg/dL)
# fbs:	Fasting blood sugar > 120 mg/dL (1 = true, 0 = false)
# restecg:	Resting electrocardiographic results
# thalach:	Maximum heart rate achieved
# exang:	Exercise-induced angina (1 = yes, 0 = no)
# oldpeak:	ST depression induced by exercise relative to rest
# slope:	Slope of the peak exercise ST segment
# ca:	Number of major vessels colored by fluoroscopy
# thal:	Thalassemia status

y = data_raw[, 1]
X = data_raw[, -1]
colnames(X) <- c(
  "age",
  "sex",
  "cp",
  "trestbps",
  "chol",
  "fbs",
  "restecg",
  "thalach",
  "exang",
  "oldpeak",
  "slope",
  "ca",
  "thal"
)
Y = scale(X)

d = ar_dist(x = Y, sim_dist = TRUE)

p = table(d$rho_sim)

dd = as.numeric(names(p))
dd = round(dd, 3)

df = data.frame(Prob = as.vector(p)/sum(as.vector(p)),
                rho = as.factor(dd))

ggplot(df) + 
  geom_bar(aes(x = rho, y = Prob), stat="identity") + 
  geom_vline(
    xintercept = 1, 
    linewidth = 1,
    linetype = 2) +
  theme(
    axis.text.x = element_text(angle = 45)
  )

d$rho
energy::mvnorm.test(Y, R = 1000)

##

Y = X[, c(1, 4, 5, 8, 10)]
Y = scale(Y)

d = ar_dist(x = Y, sim_dist = TRUE)

p = table(d$rho_sim)

dd = as.numeric(names(p))
idx = which.min(abs(dd - d$rho)) + 1
dd = round(dd, 3)

df = data.frame(Prob = as.vector(p)/sum(as.vector(p)),
                rho = as.factor(dd))

ggplot(df) + 
  geom_bar(aes(x = rho, y = Prob), stat="identity") + 
  geom_vline(
    xintercept = idx, 
    linewidth = 1,
    linetype = 2,
    color = "red") +
  theme(
    axis.text.x = element_text(angle = 45)
  )

d$rho
energy::mvnorm.test(Y, R = 1000)

plot(as.data.frame(Y))

##

Y = X[, c(1, 4, 5, 8)]
Y = scale(Y)

d = ar_dist(x = Y, sim_dist = TRUE)

p = table(d$rho_sim)

dd = as.numeric(names(p))
idx = which.min(abs(dd - d$rho)) + 1
dd = round(dd, 3)

df = data.frame(Prob = as.vector(p)/sum(as.vector(p)),
                rho = as.factor(dd))

ggplot(df) + 
  geom_bar(aes(x = rho, y = Prob), stat="identity") + 
  geom_vline(
    xintercept = idx, 
    linewidth = 1,
    linetype = 2) +
  theme(
    axis.text.x = element_text(angle = 45)
  )

d$rho
energy::mvnorm.test(Y, R = 1000)

##

hist(Y[, 1])
Y = huge::huge.npn(Y)
hist(Y[, "age"])

d = ar_dist(x = Y, sim_dist = TRUE)

p = table(d$rho_sim)

dd = as.numeric(names(p))
idx = which.min(abs(dd - d$rho)) + 1
dd = round(dd, 3)

df = data.frame(Prob = as.vector(p)/sum(as.vector(p)),
                rho = as.factor(dd))

ggplot(df) + 
  geom_bar(aes(x = rho, y = Prob), stat="identity") + 
  geom_vline(
    xintercept = idx, 
    linewidth = 1,
    linetype = 2,
    color = "red") +
  theme(
    axis.text.x = element_text(angle = 45)
  )

d$rho  
energy::mvnorm.test(Y, R = 1000)

##

library(huge)

res = huge(x = Y, 
           nlambda = 100,
           method = "glasso")

res_sel = huge.select(res, 
                      criterion = "ric")

R = -cov2cor(res_sel$opt.icov)
colnames(R) = rownames(R) = colnames(Y)

round(R, 3)

A = res_sel$opt.icov
A[A != 0] = 1
diag(A) = 0
colnames(A) = colnames(Y)
rownames(A) = colnames(Y)

G = igraph::graph_from_adjacency_matrix(A,
                                        mode = "undirected")

plot(G)

###

M = 10^3

rho_null = rep(0, M)

n = nrow(Y)

mu = rep(0, ncol(Y))

for(i in 1:M){
  
  XX = mvtnorm::rmvnorm(n = n, mean = mu)
  
  rho_null[i] = ar_dist(x = XX)$rho
  
}

(sum(rho_null <= d$rho) + 1)/(M + 1)

hist(rho_null, probability = TRUE)
abline(v = d$rho, lty = 2, lwd = 2)

###

f = function(x) ar_dist(x)$rho

p = ncol(Y)

bootobj <- boot::boot(Y, 
                      statistic = f, 
                      R = M, 
                      sim = "parametric", 
                      ran.gen = function(x, y){
                        return(matrix(rnorm(n*p), 
                                      nrow = n, ncol = p))
                      })

mean(bootobj$t < bootobj$t0)

hist(bootobj$t, add = TRUE, probability = TRUE, col = "red")

data_long = as.data.frame(Y) %>% 
  pivot_longer(colnames(Y)) %>% 
  as.data.frame()

ggp <- ggplot(data_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density))) + 
  geom_density(col = "#1b98e0", linewidth = 2) + 
  facet_wrap(~ name, scales = "free")
ggp

sc = ifelse(y == 1, "red", "blue")

plot(as.data.frame(Y),
     col = sc,
     pch = 16)

#trestbps: Patient's level of blood pressure at resting
#          mode in mm/HG

# Many repeated values -> Evidence against H0

###

Y = X[, c("age", "chol", "thalach")]
Y = huge::huge.npn(Y)

d = ar_dist(Y)

d$rho

bootobj <- boot::boot(Y, 
                      statistic = f, 
                      R = M, 
                      sim = "parametric", 
                      ran.gen = function(x, y){
                        return(matrix(rnorm(n*p), 
                                      nrow = n, ncol = p))
                      })

mean(bootobj$t < bootobj$t0)
(sum(rho_null < d$rho) + 1)/(M + 1)
energy::mvnorm.test(Y, R = M)

plot(as.data.frame(Y), 
     col = sc,
     pch = 16)

## Examine only individuals that have a heart disease

idx = which(y == 1)

Y = Y[idx, ]

d = ar_dist(Y)

d$rho

bootobj <- boot::boot(Y, 
                      statistic = f, 
                      R = M, 
                      sim = "parametric", 
                      ran.gen = function(x, y){
                        return(matrix(rnorm(n*p), 
                                      nrow = n, ncol = p))
                      })

mean(bootobj$t < bootobj$t0)
(sum(rho_null < d$rho) + 1)/(M + 1)
energy::mvnorm.test(Y, R = M)

plot(as.data.frame(Y),
     pch = 16)

## Examine only healty individuals

Y = X[, c("age", "chol", "thalach")]
Y = huge::huge.npn(Y)
idx = which(y == -1)

Y = Y[idx, ]

d = ar_dist(Y)

d$rho

bootobj <- boot::boot(Y, 
                      statistic = f, 
                      R = M, 
                      sim = "parametric", 
                      ran.gen = function(x, y){
                        return(matrix(rnorm(n*p), 
                                      nrow = n, ncol = p))
                      })

mean(bootobj$t < bootobj$t0)
(sum(rho_null < d$rho) + 1)/(M + 1)
energy::mvnorm.test(Y, R = M)

plot(as.data.frame(Y),
     pch = 16)
