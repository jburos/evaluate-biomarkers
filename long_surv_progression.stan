

/*  notes about input data 

submodels : 
   1: risk for progression
   2: risk for mortality | no progression
   3: risk for mortality | pr_progression

Variable naming:
 N            = total number of observations (length of data)
 S            = number of sample ids
 T            = max timepoint (number of timepoint ids)
 X            = number of covariates
 s            = sample id for each obs
 t            = timepoint id for each obs
 (obs)        = tumor size observed at time t for sample s
 event        = integer indicating if there was a terminating event at time t for sample s
 progression  = integer indicating if there was progression at time t for sample s
 pr_progress  = number of prior progressions by time for sample s (used to infer post-progression risk)
 covars       = matrix of real-valued covariates at time t for sample n [N, X]
                for now, all covars are used in each submodel

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
  int<lower=0, upper=1> progression[N]; 
  int<lower=0, upper=T> pr_progress[N];
  matrix[N,X] covars;
}
/* 
transformed data {
  real r;
  real c;
  r <- 0.1;
  c <- 0.01;
}
*/
parameters {
  real risk1;  // multiplier determining overall rate of progression
  real risk2;  // multiplier determining overall rate of mortality | no progression
  real risk3;  // multiplier determining overall rate mortality | progression 
  vector<lower=0>[T] base1; // baseline hazard for each timepoint t, for each submodel 
  vector<lower=0>[T] base2; // 
  vector<lower=0>[T] base3; // 
  vector[X] beta1;               // betas for each of X covariates, for each of 3 submodels
  vector[X] beta2;               // 
  vector[X] beta3;               // 
  vector[X] beta_shared;         // overall means for each beta; common between the submodels
  vector<lower=0>[S] frailty;             // subject-level frailty term
}
model {
  risk1 ~ cauchy(0, 2);
  risk2 ~ cauchy(0, 2);
  risk3 ~ cauchy(0, 3);
    
  // priors on baseline hazard function(s)
  base1 ~ normal(risk1, 1);
  base2 ~ normal(risk2, 1);
  base3 ~ normal(risk3, 1);
  
  // priors on betas (shared component)
  beta_shared ~ normal(0, 1);
  // per submodel (deviates from beta_shared)
  for (x in 1:X) {
    beta1[x] ~ normal(beta_shared[x], 1);
    beta2[x] ~ normal(beta_shared[x], 1);
    beta3[x] ~ normal(beta_shared[x], 1);
  }
  for (n in 1:N) {
    progression[n] ~ poisson(exp(covars[n,]*beta1)*base1[t[n]]*frailty[s[n]]);
    if (pr_progress[n] == 0) {
      event[n] ~ poisson(exp(covars[n,]*beta2)*base2[t[n]]*frailty[s[n]]);
    }
    if (pr_progress[n] > 0) {
      event[n] ~ poisson(exp(covars[n,]*beta3)*base3[t[n]]*frailty[s[n]]);
    }
  }
}
