# Quinn-Martin Replication
library(dplyr)
library(rstan)


setwd('~/Documentos/R/Quinn-Martin-Replication/')



data <- read.csv('bigdata.txt', sep='\t', stringsAsFactors=FALSE)


colnames(data) <- as.character(1:29)




yy <- matrix(NA, ncol=3, nrow=(3449*29))
colnames(yy) <- c("k","n","y")
control <- 1
for (i in 1:3449) {
  for (j in 1:29) {
    yy[control,1] <- i

        yy[control,2] <- j

        yy[control,3] <- data[i,j]
    
        control <- control + 1
    
  }
}


yy <- data.frame(yy)
yy <- filter(yy, y!=9)
yy <- na.omit(yy)

L <- length(unique(yy$n))
K <- length(unique(yy$k))

data = list(N=nrow(yy), L=L, K=K, jj=yy$n, kk=yy$k, y=yy$y)

stanstr <-
  '
data {
int<lower=0> N;
int<lower=0> L; // legisladores
int<lower=0> K; // rollcalls
int<lower=0> jj[N];
int<lower=0> kk[N];
int<lower=0, upper=1> y[N];
}
parameters {
real beta0[K];
real beta1[K];
vector[L] x;
}

model {
for (i in 1:N){
y[i] ~ bernoulli(Phi_approx(x[jj[i]]*beta1[kk[i]] - beta0[kk[i]]));
}

beta1 ~ normal(0, 1);
beta0 ~ normal(0, 1);
x ~ normal(0,1);  

}
'
fit <- stan(model_code = stanstr, data=data, iter=12000, warmup=2000, thin=10, chains=3)

/* for (i in 1:40){ x[i] ~ normal(0,1); }
for (i in 42:65){ x[i] ~ normal(0,1); }
for (i in 67:100){ x[i] ~ normal(0,1); }
x[41] ~ normal(-1, .0001);
x[66] ~ normal(1, .0001);
*/ 
  
