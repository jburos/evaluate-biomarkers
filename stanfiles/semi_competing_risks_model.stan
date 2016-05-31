
/*  notes about input data &/or the model

 -- dimensions are capitalized
 N              = total number of observations (length of data)
 S              = number of sample ids
 T              = max timepoint (number of timepoint ids)
 X              = number of covariates

 -- data elements are lower-case 
 s_id       = sample id for each obs
 t_id       = timepoint id for each obs
 t_dur      = duration of time for this timepoint
 t_time     = ending timepoint (time since t0) for this timepoint
 ev1        = integer (1:yes, 0:no) indicating if non-terminating event occurred at t_time
 ev2        = integer (1:yes, 0:no) indicating if terminating event occurred at t_time
 post_ev1   = integer (1:yes, 0:no) indicating if ev1 had happened previously (for h3)
 x          = matrix of numeric covariates, available to any sub-model
 
 -- indicating which covariates to use for which submodel
 x1         = integers (1:yes, 0:no) mapping covariates to submodel 1
 x2         = integers (1:yes, 0:no) mapping covariates to submodel 2
 x3         = integers (1:yes, 0:no) mapping covariates to submodel 3

*/
// Jacqueline Buros Novik <jackinovik@gmail.com>

data {
  // dimensions
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> X;
  
  // observed data
  int<lower=1, upper=S> s_id[N];
  int<lower=1, upper=T> t_id[N];
  real<lower=0> t_dur[N];
  real<lower=0> t_time[N];
  int<lower=0, upper=1> ev1[N];
  int<lower=0, upper=1> ev2[N];
  int<lower=0, upper=1> post_ev1[N];
  matrix[N,X] x;
  
  // covariate mapping to submodels
  int<lower=0, upper=1> x1[X];
  int<lower=0, upper=1> x2[X];
  int<lower=0, upper=1> x3[X];
}
transformed data {
  real t_id_dur[T];
  real t_id_end[T];
  real hprior_frailty_mean_loc;
  real hprior_frailty_sd_loc;
  real r;
  real c;
  
  c <- 0.001; // param for baseline hazard
  r <- 0.1;  // param for baseline hazard 
  
  hprior_frailty_mean_loc <- 0;
  hprior_frailty_sd_loc <- 1;
  
  for (n in 1:N) {
    t_id_dur[t_id[n]] <- t_dur[n]; // update with input t_dur (assume these are consistent for each t_id!)
    t_id_end[t_id[n]] <- t_time[n];
  }
}
parameters {
  real<lower=0> frailty_mean;
  real<lower=0> frailty_sd;
  real<lower=0> subject_frailty[S];
  vector[X] beta_m1;
  vector[X] beta_m2;
  vector[X] beta_m3;
  real<lower=0> h0_m1[T];
  real<lower=0> h0_m2[T];
  real<lower=0> h0_m3[T];
}
transformed parameters {
  vector[N] linpred1;
  vector[N] linpred2;
  vector[N] linpred3;
  vector[X] b1;
  vector[X] b2;
  vector[X] b3;
  vector[N] h1;
  vector[N] h2;
  vector[N] h3;
  
  // limit each beta to appropriate submodel
  for (k in 1:X) {
    b1[k] <- beta_m1[k] * x1[k];
    b2[k] <- beta_m2[k] * x2[k];
    b3[k] <- beta_m3[k] * x3[k];
  }

  for (n in 1:N) {
    linpred1[n] <- exp(x[n,] * b1);
    linpred2[n] <- exp(x[n,] * b2);
    linpred3[n] <- exp(x[n,] * b3);
    h1[n] <- subject_frailty[s_id[n]] * h0_m1[t_id[n]] * linpred1[n];
    h2[n] <- subject_frailty[s_id[n]] * h0_m2[t_id[n]] * linpred2[n];
    h3[n] <- subject_frailty[s_id[n]] * h0_m3[t_id[n]] * linpred3[n];
  }
}
model {
  // priors baseline hazard functions 
  for (t in 1:T) {
    h0_m1[t] ~ gamma(r * t_id_dur[t] * c, c);
    h0_m2[t] ~ gamma(r * t_id_dur[t] * c, c);
    h0_m3[t] ~ gamma(r * t_id_dur[t] * c, c);
  }
  
  // clinical event submodel
  beta_m1 ~ normal(0, 1);
  beta_m2 ~ normal(0, 1);
  beta_m3 ~ normal(0, 1);
  
  // subject-level frailty terms
  frailty_mean ~ normal(hprior_frailty_mean_loc, 1);
  frailty_sd ~ normal(hprior_frailty_sd_loc, 1);
  subject_frailty ~ normal(frailty_mean, frailty_sd);
  
  for (n in 1:N) {
    ev1[n] ~ poisson(h1[n]);
    if (post_ev1[n] == 0) {
      ev2[n] ~ poisson(h2[n]);
    } else {
      ev2[n] ~ poisson(h3[n]);
    }
  }
}
