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


T <- length(unique(yy$ter)) # number of terms
J <- length(unique(yy$jus)) # number of justices
K <- length(unique(yy$cas)) # number of cases


theta0 <- rep(0, J) # mean of priors for justices in term 0
theta0[1] <- 1 # Justice Harlam, term=0, set prior to normal(1,.01)
theta0[3] <- -3 # Justice Douglas, term=0, set prior to normal(-3,.01)
theta0[5] <- -2 # Justice Marshall, term=0, set prior to normal(-2,.01)
theta0[6] <- -2 # Justice Brennan, term=0, set prior to normal(-2,.01)
theta0[10] <- 1 # Justice Frankfurter, term=0, set prior to normal(1,.01)
theta0[14] <- -1 # Justice Fortas, term=0, set prior to normal(-1,.01)
theta0[21] <- 2 # Justice Rehnquist, term=0, set prior to normal(2,.01)
theta0[24] <- 2.5 # Justice Scalia, term=0, set prior to normal(2.5,.01)
theta0[27] <- 2.5 # Justice Thomas, term=0, set prior to normal(2.5,.01)

delta0 <- rep(1, J) # mean of priors for justices in term 0
delta0[c(1,3,5,6,10,14,21,24,27)] <- .1 # Setting variances to .1 for 
# justices Harlam, Douglas, Marshall, Brennan, Frankfurter, Fortas, Rehnquist, Scalia, Thomas

delta <- matrix(.1, ncol=J, nrow=(T-1)) # Variance=.1 for all justices and 1:T terms
delta[,3] <- .001 # Variance=.001 - Justice Douglas

b0 <- c(0,0)
B0 <- matrix(c(1,0,0,1), nrow=2,ncol=2)


###################################
### Defining Stan data and model ##
###################################
data = list(T=T,J=J,K=K,
            N=nrow(yy), tt=yy$ter, jj=yy$jus, kk=yy$cas, z=yy$y,
            theta0=theta0, delta0=delta0, delta=delta, b0=b0, B0=B0)
  
stanstr <-
  '
  /**
   *  Martin and Quinn (2002) Dynamic Ideal Point Estimation via
   *  Markov Chain Monte Carlo for the U.S. Supreme Court, 1953--1999.
  *
  *  Model section of paper is silent on fixed (?) priors for params.
  */
  data {
  int<lower=0> T;                 // number of terms
  int<lower=0> J;                 // number of justices
  int<lower=0> K;                 // total number of cases
  
  int<lower=0> N;                 // num votes
  int<lower=1> tt[N];             // term for vote n
  int<lower=1> jj[N];             // justice for vote n
  int<lower=1> kk[N];             // case for vote n
  int<lower=0, upper=1> z[N];     // vote (yea = 1, nay = 0)
  
  vector[J] theta0;               // prior mean (M&Q: mu0)
  vector<lower=0>[J] delta0;      // initial t sd (M&Q: sqrt of C[0,])
  vector<lower=0>[J] delta[T-1];  // "evolution sd"  (M&Q: sqrt of Delta[t,j])
  vector[2] b0;                   // prior mean for pos and discrim
  cov_matrix[2] B0;               // prior covar for pos and discrim
  }
  parameters {
  vector[2] alpha_beta[K];        // position and discrim of case k
  vector[J] theta[T];             // position in term t of justice j
  }
  model {
  for (n in 1:N)
  z[n] ~ bernoulli_logit(alpha_beta[kk[n], 1]
  + alpha_beta[kk[n], 2] * theta[tt[n], jj[n]]);
  
  alpha_beta ~ multi_normal(b0, B0);
  
  theta[1] ~ normal(theta0, delta0);
  for (t in 2:T)
  theta[t] ~ normal(theta[t - 1], delta[t - 1]);
  }
  '
  fit <- stan(model_code = stanstr, data=data, iter=1000, warmup=20, thin=10, chains=1)
