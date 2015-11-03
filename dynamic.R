####### REPLICATION ####
##### Martin, Andrew D., Quinn, Kevin M. 
##### "Dynamic Ideal Point Estiamtion via Markov Chain Monte Carlo for the US Supreme Court, 1953-1999"
#####
#####
#####

#### LOADING PACKAGES AND SETTING FOLDER
library(dplyr)
library(rstan)
setwd('~/Documents/R/Quinn-Martin-Replication/')




##################
# PREPARING DATA #
#################

df_votes <- data.frame(matrix(NA, nrow=3450, ncol=29)) # This dataframe will contain all the data
df_votes$T <- NA # Define a column T for term


# It's necessary to load 'adj.txt' data file to know how many justices  per term
justices_per_term <- read.csv('adj.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE) 
colnames(justices_per_term) <- 1:29
row.names(justices_per_term) <- 1:47



list_terms_files <- list.files('./vote/', pattern='txt')
initial_df_line <- 1
setwd('./vote/')
# THere is an error on 1956.txt file.
# Adj matrix refers to 10 justices, but the file contains 11 columns.
# We're gonna correct it excluding the last column of the file.
# Verify with authors


for (i in 1:length(list_terms_files)) {
    file <- list_terms_files[i]
    cases_term <- read.csv(file, header=FALSE, sep='\t',stringsAsFactors=FALSE) 
    current_justices <- which(justices_per_term[i,]!=-9)
    final_df_line <- initial_df_line + nrow(cases_term) - 1
    df_votes[initial_df_line:final_df_line,current_justices] <- cases_term
    df_votes$T[initial_df_line:final_df_line] <- i
    initial_df_line <- final_df_line + 1
}


# It's not easy to get along with missing data on Stan as it is on Jags so we must transform data
yy <- matrix(NA, ncol=4, nrow=(3450*29))
colnames(yy) <- c("cas","jus","y", "ter")
control <- 1
for (i in 1:3450) {
   for (j in 1:29) {
    yy[control,1] <- i
    yy[control,2] <- j
    yy[control,3] <- df_votes[i,j]
    yy[control,4] <- df_votes[i, 30]
    control <- control + 1
  }
}
yy <- data.frame(yy)
yy <- filter(yy, y!=9)
yy <- na.omit(yy)


# Creating index for parameters
xx <- data.frame(unique(cbind(yy$jus, yy$ter)))
colnames(xx) <- c("jus","ter")
xx <- arrange(xx, jus,ter)
par_function <- Vectorize(function(jus, ter) {
  return(rownames(xx[(xx[,1]==jus & xx[,2]==ter),]))
})
yy <- mutate(yy, par=par_function(jus, ter))
# yy$par is a index for each parameter x such that x_it, with i relating to justice and t to term



# define a matrix zz, with initial term and final term
zz <- matrix(NA, nrow=length(unique(yy$jus)), ncol=2)
colnames(zz) <- c("jus", "interval")
for ( j in 1:length(unique(yy$jus))) {
  k <- unique(yy$jus)[j]
  zz[j,1] <- k
  temp <- filter(yy, jus==k) %>%
    arrange(ter)
  zz[j,2] <- tail(temp,1)$ter - head(temp,1)$ter
  }
zz <- as.data.frame(zz)
vec <- 0
for (i in 1:(2*nrow(zz))) { # Verdadeira gambiarra
  print(i)
  if ( ( i %% 2 )==1 ) {
    vec <- c(vec, (tail(vec,1)+1))
  }
  else {
    vec <- c(vec, (tail(vec,1) + zz$interval[i/2]))    
  }
}
vec <- vec[2:length(vec)]
vec <- as.data.frame(t(matrix(vec, nrow=2)))
colnames(vec) <- c("in_ter", "fi_ter")


n_cas <- length(unique(yy$cas))
n_jus <- length(unique(yy$jus))
n_par <- length(unique(yy$par)) # There is a parameter for term and justice



yy$par <- as.numeric(yy$par)


###################################
### Defining Stan data and model ##
###################################
data = list(N=nrow(yy), ncas=n_cas, njus=n_jus, npar=n_par, 
            cas=yy$cas, par=yy$par, in_ter=vec$fi_ter, fi_ter= vec$fi_ter, y=yy$y )

stanstr <-
'
data {
int<lower=0> N;
int<lower=0> ncas; // number of cases
int<lower=0> njus; // number of justices 
int<lower=0> npar; // number of parameters
int<lower=0> cas[N];
int<lower=0> par[N];
int<lower=0> in_ter[njus];
int<lower=0> fi_ter[njus];
int<lower=0, upper=1> y[N];
}
parameters {
real beta0[ncas];
real beta1[ncas];
vector[npar] x;
real sigma[njus];
}

model {

  for (i in 1:N){
  y[i] ~ bernoulli(Phi_approx(x[par[i]]*beta1[cas[i]] - beta0[cas[i]]));
  }
  
  
  for (i in 1:ncas) { 
  beta1 ~ normal(0, 1);
  }
  for (i in 1:ncas) {
  beta0 ~ normal(0, 1);
  }
  for (j in 1:njus) {
    if ( j == 3 ) {
      x[in_ter[j]] ~ normal(-3,.1); // Justice Douglas, term=0, set prior to normal(-3,.01)
    }
    else if ( j == 1 ) {
      x[in_ter[j]] ~ normal(1,.1); // Justice Harlam, term=0, set prior to normal(1,.01)
    }
    else if ( j == 5 ) {
      x[in_ter[j]] ~ normal(-2,.1); // Justice Marshall, term=0, set prior to normal(-2,.01)
    }
    else if ( j == 6 ) {
      x[in_ter[j]] ~ normal(-2,.1); // Justice Brennan, term=0, set prior to normal(-2,.01)
    }
    else if ( j == 10 ) {
      x[in_ter[j]] ~ normal(1,.1); // Justice Frankfurter, term=0, set prior to normal(1,.01)
    }
    else if ( j == 14 ) {
      x[in_ter[j]] ~ normal(-1,.1); // Justice Fortas, term=0, set prior to normal(-1,.01)
    }
    else if ( j == 21 ) {
      x[in_ter[j]] ~ normal(2,.1); // Justice Rehnquist, term=0, set prior to normal(2,.01)
    }
    else if ( j == 24 ) {
      x[in_ter[j]] ~ normal(2.5,.1); // Justice Scalie, term=0, set prior to normal(2.5,.01)
    }
    else if ( j == 27 ) {
      x[in_ter[j]] ~ normal(2.5,.1); // Justice Thomas, term=0, set prior to normal(2.5,.01)
    }

    else { 
      x[in_ter[j]] ~ normal(0,1);
    }
    if ( in_ter[j]!=fi_ter[j]) { 
      for (k in (in_ter[j]+1):fi_ter[j]) {
        x[k] ~ normal(x[k-1],.01);
      }
    }
  }
  
}
'
fit <- stan(model_code = stanstr, data=data, iter=1000, warmup=20, thin=10, chains=1)


la <- extract(fit, permuted = TRUE)

x <- la$x
x <- apply(x, 2, mean)


results <- matrix(NA, ncol=length(unique(yy$ter)), nrow=n_jus)

for (i in 1:29) {
  for (j in 1:47) {
    pare <- yy$par[yy$jus==i & yy$ter==j][1]
    if ( length(pare)!=0) 
      results[i,j] <- x[pare]
  }
}


harlam <- results[1,]
plot(index=1:47, harlam, type="l")

rehnquist <- results[21,]
plot(index=1:47, rehnquist, type="l")

harlam <- results[24,]
plot(index=1:47, harlam, type="l")

harlam <- results[1,]
plot(index=1:47, harlam, type="l")

harlam <- results[1,]
plot(index=1:47, harlam, type="l")
