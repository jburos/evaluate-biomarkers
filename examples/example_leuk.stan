## from : https://github.com/stan-dev/example-models/blob/master/bugs_examples/vol1/leuk/leuk.stan

data {
  int<lower=0> N;       // number of samples 
  int<lower=0> NT;      // number of unique failure timepoints
  int<lower=0> obs_t[N];   // observed timepoints for each sample i
  int<lower=0> t[NT + 1]; // distinct failure timepoints (plus last timepoint)
  int<lower=0> fail[N]; // whether sample i failed or not at obs_t
  real Z[N];     // covariate (I assume this is a treatment indicator)
}
transformed data {
  int Y[N, NT];  
  int dN[N, NT]; 
  real c;
  real r; 
  for(i in 1:N) {
    for(j in 1:NT) {
      Y[i, j] <- int_step(obs_t[i] - t[j] + .000000001); // was observation observed / at risk at timepoint j for subject i 
      dN[i, j] <- Y[i, j] * fail[i] * int_step(t[j + 1] - obs_t[i] - .000000001); // failure at time j for subject i
    }
  }
  c <- 0.001; // param for baseline hazard
  r <- 0.1;  // param for baseline hazard 
}
parameters {
  real beta;  // coefficients 
  real<lower=0> dL0[NT];  // baseline hazard for each timepoint
} 
model {
  beta ~ normal(0, 1000);
  for(j in 1:NT) {
    dL0[j] ~ gamma(r * (t[j + 1] - t[j]) * c, c); // prior : gamma(r * duration * c, c)
    for(i in 1:N) {
      if (Y[i, j] != 0)  
        increment_log_prob(poisson_log(dN[i, j], Y[i, j] * exp(beta * Z[i]) * dL0[j])); 
        // yhat : failure at time j for patient i
        // linpred: if-at-risk * exp(X*beta) * baseline_hazard 
    }     
  }
}
generated quantities {
  real S_placebo[NT];
  real S_treat[NT];

  for (j in 1:NT) {
    // Survivor function = exp(-Integral{l0(u)du})^exp(beta*z)
    real s;
    s <- 0;
    for (i in 1:j)
      s <- s + dL0[i];
    S_treat[j] <- pow(exp(-s), exp(beta * -0.5));
    S_placebo[j] <- pow(exp(-s), exp(beta * 0.5));      
  }
  
}
