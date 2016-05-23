
/*  notes about input data &/or the model

Variable naming:
 N            = total number of observations (length of data)
 S            = number of sample ids
 T            = max timepoint (number of timepoint ids)
 X            = number of covariates
 s            = sample id for each obs
 t            = timepoint id for each obs
 failure      = integer indicating if there was a failure event at time t for sample s
 progression  = integer indicating if there was progression at time t for sample s
 risk_covars  = real matrix of risk factors for failure/progression
 
 -- tumor-measurement timepoints --
 N_obs        = int; number of measurement occasions (length of data)
 obs_s        = int; sample id for each measurement
 obs_t        = int; timepoint (id) for each measurement occasion
 obs_size     = real; tumor measurement at time t (use 0 when no measurement)

( size_covars  = real matrix of covars potentially influencing rate of growth )

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
   *   d.TV / d.t  =  (TV*growth*(max_size-TV)/max_size)
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

    dTV_dt[1] <- (TV[1]*growth_rate*(max_size-TV[1])/max_size); 
    
    return dTV_dt; 
  } 

}
data {
  
  // inputs to the event model
//  int<lower=1> N;
    int<lower=1> S;
//  int<lower=1> T;
//  int<lower=0> X;
//  int<lower=1, upper=N> s[N];
//  int<lower=1, upper=T> t[N];
//  int<lower=0, upper=1> event[N];
//  int<lower=0, upper=1> progression[N]; 
//  matrix[N,X] risk_covars;

  // inputs to the growth model
  int<lower=0> N_obs;
  int<lower=0> T_obs;
  real<lower=0> unique_T[T_obs];
  int<lower=0> obs_s[N_obs];
  int<lower=0> obs_t[N_obs];
  real<lower=0> obs_size[N_obs];
}
transformed data { 
  real x_r[0];                 // no real data for ODE system 
  int x_i[0];                  // no integer data for ODE system 
  int t0;
  t0 <- 1;
}
parameters {
  real<lower=0> max_size;            // max tolerable size for each tumor
  real<lower=0> growth_rate[S];      // each patient has avg rate of growth over time
  real<lower=0> meas_error;          // global measurement error of diameters
  real<lower=0> init_vol[S];
 // real<lower=0> risk[S];             // each patient has local risk of failure / progression
 // real<lower=0> baseline_surv[T];    // baseline hazard for failure model
 // real<lower=0> baseline_progression[T]; // baseline hazard for progression model
}
transformed parameters {
  real<lower=0> tumor_vol[N_obs];  // inferred tumor size
  real<lower=0> vol_change[N_obs];  // change tumor size
  //real<lower=0> tumor_diam[N_obs]; // inferred tumor diam
  
  for (n in 1:N_obs) {
    real vol_change;
    real curr_state[1];
    real theta[2];
    theta[1] <- growth_rate[obs_s[n]];
    theta[2] <- max_size;
    if (obs_t[n]==0) {
      curr_state[1] <- init_vol[obs_s[n]];
      vol_change[n] <- integrate_ode(tumor_growth  
                           , curr_state
                           , t0
                           , unique_T
                           , theta
                           , x_r
                           , x_i
                           )[1, obs_t[n]];
      tumor_vol[n] <- init_vol[obs_s[n]] + vol_change[n];
    }
    else { // assume data are sorted!
      curr_state[1] <- tumor_vol[n-1];
      vol_change[n] <- integrate_ode(tumor_growth  
                                 , curr_state
                                 , unique_T[obs_t[n-1]]
                                 , unique_T
                                 , theta
                                 , x_r
                                 , x_i
                                 )[1, obs_t[n]];
      tumor_vol[n] <- tumor_vol[n-1] + vol_change[n];
    }
    //tumor_diam[n] <- cbrt(tumor_vol[n] * (0.75) / pi());  // not sure if this transformation is allowed
  }
}
model {
  meas_error ~ normal(0, 1);
  obs_size ~ normal(tumor_vol, meas_error);
}
