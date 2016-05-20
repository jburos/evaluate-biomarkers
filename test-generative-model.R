rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstan)
library(survival)
library(arm)
library(lme4)
rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(),3))

## ------ prep & review simulated data -------

source('simulate-data.function.R')
source('prep-data.function.R')
source('make-data-plots.function.R')

## get simulated data 
data <- simulate_data(n = 100, max_size = 4000, max_t = 50, failure_threshold = 4, progression_threshold = 3)

## review graphics for simulated data
make_data_plots(data)

## prep data for analysis
res <- prep_data(data) 
survd <- res$per_patient ## summarized per patient; appropriate for typical survival analysis
adata <- res$per_observation ## denormalized; appropriate for longitudinal analysis
rm(res)

## ------ review survival model (os) using coxph -------

survfit <- coxph(
  formula = Surv(first_failure, failure_status) ~ rescaled_init_size
  , data = survd
  )
print(survfit)

## ------ review survival model (pfs) using coxph -------

survfit <- coxph(
  formula = Surv(first_failure_or_progression, failure_or_progression_status) ~ rescaled_init_size
  , data = survd
  )
print(survfit)

## ------ review growth model using lmer -------

library(lme4)
library(arm)
growthfit <- lmer(rescaled_patient_observed_size ~ t + rescaled_init_size + 
                    (1 + t | patid)
                  , data = adata
                  )
display(growthfit)        

## ------ first version of generative model using stan -------


sample_data <- adata %>% semi_join(adata %>% sample_n(1) %>% dplyr::select(patid), by = 'patid')
### testing model for change in tumor size over time 
## measured with error as "diameter" of the tumor
growthdata <- list(
  N_obs = nrow(sample_data)
  , obs_s = sample_data$patid
  , obs_t = sample_data$t
  , obs_size = sample_data$tumor_size
)
testfit <- stan('generative_model2.stan', data = growthdata, chains = 1, iter = 10)
