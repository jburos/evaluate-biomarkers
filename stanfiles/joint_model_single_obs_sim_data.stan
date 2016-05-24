
/*  notes about input data &/or the model

 -- measurement/growth submodel
 N_obs        = int; number of measurement occasions (length of data)
 obs_t        = int; timepoint (id) for each measurement occasion
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
  real<lower=1> obs_s[N_obs];
  real<lower=0> obs_t[N_obs]; // timepoint ids
  real<lower=0> obs_size[N_obs];

  // clinical event submodel
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> X;
  int<lower=1, upper=N> s[N];
  int<lower=1, upper=T> t_id[N];
  real<lower=0> t[N];
  int<lower=0, upper=1> event[N];
  matrix[N,X] event_covars;
}
transformed data { 
  real x_r[0];                 // no real data for ODE system 
  int x_i[0];                  // no integer data for ODE system 
  int t0;
  real r;
  real c;
  
  t0 <- 0;
  r <- 0.1;     // priors for event submodel
  c <- 0.01;    // priors for event submodel
}
parameters {
  real<lower=0> max_size;
  real<lower=0> growth_rate;         // rate of growth
  real<lower=0> meas_error;             // global measurement error of diameters
  real<lower=0> init_vol;            // estimated initial volume
  real beta_tumor_vol;
  real intercept;
  vector[T] baseline; // unstructured baseline hazard for each timepoint t
  //real sample_frailty[S]; 
  vector[X] beta; // beta for each covariate
}
transformed parameters {
  real obs_tumor_vol[N_obs, 1];  // inferred tumor size
  real obs_tumor_diam[N_obs];
  vector<lower=0>[N] hazard;
  real hz_tumor_vol[N, 1];
  real obs_theta[2];
  real obs_init_state[1];
  real hz_theta[2];
  real hz_init_state[1];
  
  obs_theta[1] <- growth_rate;
  obs_theta[2] <- max_size;
  obs_init_state[1] <- init_vol;
  
  // observed tumor volumes given parameters
  obs_tumor_vol <- integrate_ode(tumor_growth
                                , obs_init_state
                                , t0
                                , obs_t
                                , obs_theta
                                , x_r
                                , x_i
                                );
  
  // lp for observed tumor diameters
  for (n in 1:N_obs) 
    obs_tumor_diam[n] <- 2*cbrt(obs_tumor_vol[n,1] * (0.75) / pi());
  
  // inferred tumor volumes at hazard times, given parameters
  hz_init_state[1] <- init_vol;
  hz_theta[1] <- growth_rate;
  hz_theta[2] <- max_size;
  hz_tumor_vol <- integrate_ode(tumor_growth
                                , hz_init_state
                                , t0
                                , t
                                , hz_theta
                                , x_r
                                , x_i
                                );
  for (n in 1:N) {
    hazard[n] <- inv_logit(intercept + baseline[t_id[n]] + event_covars[n,]*beta + hz_tumor_vol[n,1]*beta_tumor_vol);
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
