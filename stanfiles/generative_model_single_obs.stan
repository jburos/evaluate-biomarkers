
/*  notes about input data &/or the model

 N_obs        = int; number of measurement occasions (length of data)
 obs_t        = int; timepoint (id) for each measurement occasion
 obs_size     = real; tumor measurement at time t (use 0 when no measurement)

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

    dTV_dt[1] <- (TV[1]*growth_rate*(max_size-TV[1])); 
    
    return dTV_dt; 
  } 

}
data {
  int<lower=0> N_obs;
  real<lower=0> obs_t[N_obs]; // timepoint ids
  real<lower=0> obs_size[N_obs];
}
transformed data { 
  real x_r[0];                 // no real data for ODE system 
  int x_i[0];                  // no integer data for ODE system 
  int t0;
  t0 <- 0;
}
parameters {
  real<lower=0> max_size_raw;        // max tolerable size for each tumor
  real<lower=0> growth_rate;         // rate of growth
  real<lower=0> meas_error;          // global measurement error of diameters
  real<lower=0> init_vol;            // estimated initial volume
}
transformed parameters {
  real tumor_vol[N_obs, 1];  // inferred tumor size
  real theta[2];
  real init_state[1];
  real max_size;
  max_size <- 4000 + max_size_raw;
  theta[1] <- growth_rate;
  theta[2] <- max_size;
  init_state[1] <- init_vol;
  
  
  tumor_vol <- integrate_ode(tumor_growth
                             , init_state
                             , t0
                             , obs_t
                             , theta
                             , x_r
                             , x_i
                             );
}
model {
  meas_error ~ cauchy(0, 1);
  init_vol ~ normal(0, 1);
  growth_rate ~ normal(0, 1);
  max_size_raw ~ normal(0, 10);
  for (n in 1:N_obs)
    obs_size[n] ~ normal(tumor_vol[n,1], meas_error);
}
