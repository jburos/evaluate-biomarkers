
/*  notes about input data &/or the model

 -- measurement/growth submodel
 N_obs        = int; number of measurement occasions (length of data)
 obs_t        = int; timepoint (id) for each measurement occasion
 obs_s        = int; sample id for each measurement occasion
 obs_size     = real; tumor measurement at time t (use 0 when no measurement)

 -- event submodel
 N              = total number of observations (length of data)
 S              = number of sample ids
 T              = max timepoint (number of timepoint ids)
 X              = number of covariates
 s              = sample id for each obs
 t              = timepoint id for each obs
 event          = integer indicating if there was an event at time t for sample s
 event_covars   = matrix of real-valued covariates at time t for sample n [N, X]

*/
// Jacqueline Buros Novik <jackinovik@gmail.com>

functions { 
  /** 
   * System State is represented by the current volume (TV) of the tumor 
   * 
   * The system has parameters 
   * 
   *   theta = (growth, max_size) 
   * 
   * where 
   * 
   *   growth:   inherent growth rate 
   *   max_size: max tolerable size for the tumor (constant)
   * 
   * The time derivative is 
   * 
   *   d.TV / d.t  =  TV*growth*(1-TV/max_size)
   * 
   * @param t time at which derivative is evaluated. 
   * @param TV tumor volume at time of evaluation.
   * @param theta parameters for system. 
   * @param x_r real constants for system (stan defaults; empty). 
   * @param x_i integer constants for system (stan defaults; empty). 
   */ 
  real[] tumor_growth(real t, real[] TV, real[] theta, 
                           real[] x_r, int[] x_i) { 
    real growth_rate;
    real max_size;

    real dTV_dt[1]; 
 
    growth_rate <- theta[1]; 
    max_size <- theta[2]; 

    dTV_dt[1] <- (TV[1]*growth_rate*(1-TV[1]/max_size)); 
    
    return dTV_dt; 
  } 

}
data {
  // observed tumor size submodel
  int<lower=0> N_obs;
  int<lower=1> obs_s[N_obs];
  int<lower=0> obs_t[N_obs]; // timepoint ids
  real<lower=0> obs_size[N_obs];

  // clinical event submodel
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> X;
  int<lower=1, upper=S> s[N];
  int<lower=1, upper=T> t[N];
  int<lower=0, upper=1> event[N];
  matrix[N,X] event_covars;
}
transformed data { 
  real x_r[0];                 // no real data for ODE system 
  int x_i[0];                  // no integer data for ODE system 
  real seq_t[T];
  int t0;
  t0 <- 0;
  for (time in 1:T)
    seq_t[time] <- time;
}
parameters {
  real<lower=0> max_size[S];
  real<lower=0> growth_rate[S];         // rate of growth
  real<lower=0> meas_error;             // global measurement error of diameters
  real<lower=0> init_vol[S];            // estimated initial volume
  real beta_tumor_vol;
  real intercept;                       // overall hazard 
  vector[T] baseline;                   // baseline hazard for each timepoint t
  vector[X] beta;                       // beta for each covariate
}
transformed parameters {
  vector<lower=0>[N] hazard;
  matrix[S,T] tumor_vol;                // estimated tumor volumes for each sample / timepoint
  matrix[S,T] tumor_diam;               // estimated tumor diameters for each sample / timepoint
  real obs_tumor_diam[N_obs];
  
  // populate inferred vol & diam for each timepoint / sample combination
  for (sample in 1:S) {
    real sample_tumor_vol[T, 1];  // inferred tumor size
    real sample_theta[2];
    real sample_init_state[1];
    
    sample_theta[1] <- growth_rate[sample];
    sample_theta[2] <- max_size[sample];
    sample_init_state[1] <- init_vol[sample];
    
    // observed tumor volumes given parameters
    sample_tumor_vol <- integrate_ode(tumor_growth
                                    , sample_init_state
                                    , t0
                                    , seq_t
                                    , sample_theta
                                    , x_r
                                    , x_i
                                    );
    
    tumor_vol[sample,] <- to_row_vector(sample_tumor_vol[,1]);
    for (time in 1:T)
      tumor_diam[sample, time] <- 2*cbrt(tumor_vol[sample, time] * (0.75) / pi());
  }
  
  // hazard estimate for each event observation
  for (n in 1:N) {
    hazard[n] <- inv_logit(intercept + baseline[t[n]] + event_covars[n,]*beta + tumor_vol[s[n],t[n]]*beta_tumor_vol);
  }
  
  // diameter estimate for each sample & diameter measurement
  for (n in 1:N_obs) {
    obs_tumor_diam[n] <- tumor_diam[obs_s[n],obs_t[n]];
  }

}
model {
  // observed size submodel
  meas_error ~ cauchy(0, 1);
  init_vol ~ normal(0, 1);
  growth_rate ~ normal(0, 1);
  max_size ~ normal(0, 10);
  obs_size ~ normal(obs_tumor_diam, meas_error);
  
  // clinical event submodel
  beta_tumor_vol ~ normal(0, 1);
  intercept ~ normal(-13, 1);
  baseline ~ normal(0, 1);
  beta ~ normal(0, 1);
  event ~ bernoulli(hazard);
}
generated quantities {
  int y_hat[N];
  real loglik[N];
  
  for (n in 1:N) {
    if (!is_nan(hazard[n]) && hazard[n] <= pow(2,30)) {
      y_hat[n] <- bernoulli_rng(hazard[n]);
      loglik[n] <- bernoulli_log(event[n], hazard[n]);
    } else {
      y_hat[n] <- 0;
      loglik[n] <- not_a_number();
    }
  }
}
