
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
  real[] tumor_growth(real t, real[] y, real[] theta, 
                           real[] x_r, int[] x_i) { 
    real growth_rate;
    real max_size;

    real dy_dt[1]; 
 
    growth_rate <- theta[1]; 
    max_size <- theta[2]; 

    dy_dt[1] <- y[1]*growth_rate*((max_size-y[1])/max_size); 
    
    return dy_dt; 
  } 

}
data {
  // N obs to simulate
  int<lower=0> N_obs;
  real<lower=0> obs_t[N_obs]; // timepoint ids
  real<lower=0> init_vol;
  real<lower=0> growth_rate;
  real<lower=0> max_size;
}
transformed data { 
  real x_r[0];                 // no real data for ODE system 
  int x_i[0];                  // no integer data for ODE system 
  real t0;
  t0 <- 0;
}
model {
  
}
generated quantities {
  real tumor_vol[N_obs, 1];  // inferred tumor size
  real theta[2];
  real init_state[1];
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

