/*  Variable naming:
 N          = total number of observations (length of data)
 S          = number of sample ids
 T          = max timepoint (number of timepoint ids)
 X          = number of covariates
 s          = sample id for each obs
 t          = timepoint id for each obs
 event      = integer indicating if there was an event at time t for sample s
 covars     = matrix of real-valued covariates at time t for sample n [N, X]
*/
// Jacqueline Buros Novik <jackinovik@gmail.com>

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> X;
  int<lower=1, upper=N> s[N];
  int<lower=1, upper=T> t[N];
  //int<lower=0, upper=1> obs[N];
  int<lower=0, upper=1> event[N];
  matrix[N,X] covars;
}
transformed data {
  real r;
  real c;
  r <- 0.1;
  c <- 0.01;
}
parameters {
  vector<lower=0>[T] baseline; // unstructured baseline hazard for each timepoint t
  //real sample_frailty[S]; 
  vector[X] beta; // beta for each covariate
}
transformed parameters {
  vector<lower=0>[N] hazard;
  
  for (n in 1:N) {
    hazard[n] <- exp(covars[n,]*beta)*baseline[t[n]];
  }
}
model {
  // baseline ~ gamma(r * c, c);
  baseline ~ normal(0, 1);
  //sample_frailty ~ normal(0, 1);
  beta ~ normal(0, 1);
  event ~ poisson(hazard);
}
