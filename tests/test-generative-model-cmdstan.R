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


## ------ test out ode using gen quant block ------ 

sample_data <- adata %>% semi_join(adata %>% sample_n(1) %>% dplyr::select(patid), by = 'patid')

sample_params <- list(
  N_obs = nrow(sample_data)
  , obs_t = sample_data$t
  , init_vol = unique(sample_data$init_size)
  , growth_rate = unique(sample_data$growth_rate)
  , max_size = 4000
)
testfit <- stan('generative_model_sim_data.stan', data = sample_params, chains = 1, iter = 100, algorithm = 'Fixed_param')
print(testfit)

## ------ first version of generative model using stan -------

### testing model for change in tumor size over time 
## measured with error as "diameter" of the tumor
growthdata <- list(
  N_obs = nrow(sample_data)
  , obs_t = sample_data$t
  , obs_size = sample_data$tumor_size
  , max_size = 4000
)
## call using cmdstan instead of rstan. more reliable for longer-running models
stan_home <- '/usr/local/Cellar/cmdstan/2.9.0'
modelpath <- file.path(getwd(),'generative_model_single_obs_more_params.stan')
modelpath <- gsub(modelpath, pattern = "(.*)\\.stan", replacement = '\\1')
modelname <- gsub(modelpath, pattern = ".*\\/([^\\/.]+)$", replacement = "\\1")
datafile <- file.path(getwd(),paste0(modelname,'.data.R'))
## write data to disk
with(growthdata, stan_rdump(names(growthdata), file = datafile))
## translate stan -> c++; compile to executable
system(paste0('(cd ',stan_home,' && make ',modelpath,')'))
## test on 10 iterations

system(paste0('./',modelname,' diagnose data file=',datafile))
system(paste0('./',modelname,' sample num_samples=5 num_warmup=5 random seed=12345 data file=',datafile))
#growthtest <- stan('generative_model_single_obs.stan', data = growthdata, chains = 1, iter = 10)
#growthfit <- stan('generative_model_single_obs.stan', data = growthdata, chains = 3, iter = 1000)

## call using Rstan 
growthfit_fixed <- stan(file = 'generative_model_single_obs_more_params.stan', chains = 1, iter = 10, seed = 12345, data = growthdata)

