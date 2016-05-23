rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstan)
library(purrr)
library(deSolve)
rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(),3))

## ------ prep simulated data -------

source('simulate-data.function.R')
source('prep-data.function.R')
source('make-data-plots.function.R')

## get simulated data 
## changes from previous versions:
##   1. remove all sources of noise
##   2. make hazard directly proportional to size of tumor
data <- simulate_data(n = 100
                      , max_size = 4000
                      , max_t = 50
                      , failure_threshold = 4
                      , progression_threshold = 3
                      , size_noise_fun = create_scalar(0)
                      , growth_rate_noise_fun = create_scalar(0)
                      , hazard_noise_fun = create_scalar(0)
                      , hazard_fun = function(row) {row$tumor_size} ## for now, hazard proportional to size
                      )

## review data for a few simulated points
plot_simulated_data(data, n = 6)

## prep data for analysis
res <- prep_data(data) 
survd <- res$per_patient ## summarized per patient; appropriate for typical survival analysis
adata <- res$per_observation ## denormalized; appropriate for longitudinal analysis
rm(res)

## ------ test growth model using lmer -------

library(lme4)
library(arm)
growthfit <- lmer(rescaled_patient_observed_size ~ t + rescaled_init_size + 
                    (1 + t || patid)
                  , data = adata
                  )
display(growthfit)        

## how well does this fit?

adata$pd <- predict(growthfit)
ggplot(adata %>% semi_join(adata %>% sample_n(1), by = 'patid')
       , aes(x = t, group = patid)
       ) +
  geom_line(aes(y = rescaled_patient_observed_size)) + 
  geom_point(aes(y = pd))

## need to transform data, since we are modeling a % growth rate
growthfit2 <- lmer(log1p(observed_size) ~ t + rescaled_init_size + 
                    (1 + t || patid)
                  , data = adata
)
display(growthfit2)        
adata$pd2 <- expm1(predict(growthfit2))
ggplot(adata %>% semi_join(adata %>% sample_n(1), by = 'patid')
       , aes(x = t, group = patid)
) +
  geom_line(aes(y = observed_size)) + 
  geom_point(aes(y = pd2))

## ------ test growth model using stan ------ 

## pick random patient
sample_data <- adata %>% semi_join(adata %>% sample_n(1) %>% dplyr::select(patid), by = 'patid')
plot_simulated_data(sample_data, n = NULL)

## what does simulated data look like according to these params?
## make sure data simulated according to R match those according to Stan
sample_params <- list(
  N_obs = nrow(sample_data)
  , obs_t = sample_data$t
  , init_vol = unique(sample_data$init_size)
  , growth_rate = unique(sample_data$growth_rate)
  , max_size = 4000
)
stangen <- stan('generative_model_sim_data.stan', data = sample_params, chains = 1, iter = 5, algorithm = 'Fixed_param')
print(stangen, pars = 'tumor_vol')
ppd_vol <- extract(stangen, 'tumor_vol')$tumor_vol
ppd_diam <- extract(stangen, 'tumor_diam')$tumor_diam
sample_data$vol_from_stan <- apply(ppd_vol, FUN = unique, MARGIN = 2)
sample_data$diam_from_stan <- apply(ppd_diam, FUN = unique, MARGIN = 2)
ggplot(sample_data, aes(x = t)) + 
  geom_line(aes(y = tumor_size, colour = 'simulated - R')) + 
  geom_line(aes(y = vol_from_stan, colour = 'simulated - stan'))

ggplot(sample_data, aes(x = t)) + 
  geom_line(aes(y = observed_size, colour = 'simulated - R')) + 
  geom_line(aes(y = diam_from_stan, colour = 'simulated - stan'))


## perfect agreement!
ggplot(sample_data, aes(x = tumor_size, y = vol_from_stan)) + 
  geom_point() + 
  ggtitle ('Comparing simulated data using Stan & R (volumes)')
  
ggplot(sample_data, aes(x = observed_size, y = diam_from_stan)) + 
  geom_point() + 
  ggtitle ('Comparing simulated data using Stan & R (diameters)')


## now, estimate model for these data, given parameters
standata <- list(
  N_obs = nrow(sample_data)
  , obs_t = sample_data$t
  , obs_size = sample_data$tumor_size
  , max_size = 4000
)

testfit <- stan('generative_model_single_obs_more_params.stan', data = standata, iter=10, chains = 1)

stanfit1 <- stan(fit = testfit, data = standata, iter = 1000, chains = 3)

print(stanfit1, pars = c('init_vol','growth_rate'))
print(sample_params$init_vol)
print(sample_params$growth_rate)

## can we estimate the 'max_size' from the data?
standata2 <- list(
  N_obs = nrow(sample_data)
  , obs_t = sample_data$t
  , obs_size = sample_data$tumor_size
)

testfit2 <- stan('generative_model_single_obs2.stan', data = standata2, iter=10, chains = 1)

stanfit2 <- stan(fit = testfit2, data = standata2, iter = 500, chains = 3)

print(stanfit2, pars = c('init_vol','growth_rate','max_size'))
print(sample_params$init_vol)
print(sample_params$growth_rate)

## how about the case where we measure diameter instead of volume?
standata3 <- list(
  N_obs = nrow(sample_data)
  , obs_t = sample_data$t
  , obs_size = sample_data$observed_size
)

test_init <- list()
for (chain in 1:1) {
  test_init[[chain]] <- list(
    meas_error = abs(rnorm(n = 1, mean = 0, sd = 0.1))
    , init_vol = abs(rnorm(n = 1, mean = 0, sd = 1))
    , growth_rate = abs(rbeta(n = 1, 10, 20))
    , max_size_raw = 0.1
  )
}
testfit3 <- stan('generative_model_single_obs_diam.stan', data = standata3, iter=10, chains = 1)

init3 <- list()
for (chain in 1:3) {
  init3[[chain]] <- list(
    meas_error = abs(rnorm(n = 1, mean = 0, sd = 0.1))
    , init_vol = abs(rnorm(n = 1, mean = 0, sd = 1))
    , growth_rate = abs(rbeta(n = 1, 10, 20))
    , max_size_raw = 0.1
  )
}
stanfit3 <- stan(fit = testfit3, data = standata3, iter = 1000, chains = 3)




