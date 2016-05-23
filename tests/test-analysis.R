rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(),3))

source('simulate-data.function.R')

## get simulated data 
data <- simulate_data(n = 100, max_size = 4000)
plot_simulated_data(data, n = 12)

## distribution of survival times
ggplot(data %>% 
         dplyr::filter(observed == 1) %>%
         group_by(patid) %>% 
         dplyr::summarise(surv_time  = max(t, na.rm = T)
                          , event = factor(last(failure_event, order_by = t), levels = c(0,1), labels = c('censor', 'failure'))
                          ) %>%
         ungroup()
       , aes(x = surv_time, group = event, fill = event, colour = event)) +
  geom_histogram(alpha = 0.2) +
  facet_wrap(~event)

## summarize inputs to survival analysis (using MLE)
survd <- data %>%
  dplyr::filter(observed == 1) %>%
  dplyr::group_by(patid) %>%
  dplyr::mutate(
    surv_time  = max(t, na.rm = T)
    , init_observed_size = first(observed_size, order_by = t)
    , last_observed_size = last(observed_size, order_by = t)
    , outcome = last(failure_event, order_by = t)
  ) %>%
  dplyr::filter(t == surv_time) %>%
  ungroup() %>%
  dplyr::mutate(avg_growth_rate = (last_observed_size - init_observed_size)/t
                , mean_growth_rate = mean(avg_growth_rate)
                , sd_growth_rate = sd(avg_growth_rate)
                , rescaled_growth_rate = (avg_growth_rate - mean_growth_rate)/sd_growth_rate
                , mean_init_size = mean(init_observed_size)
                , sd_init_size = sd(init_observed_size)
                , rescaled_init_size = (init_observed_size - mean_init_size)/sd_init_size
                ) %>%
  dplyr::select(patid, t, init_observed_size, last_observed_size, avg_growth_rate, outcome, starts_with('rescaled'))

survd %>%
  group_by(outcome) %>%
  tally()

## if only one obs in censor or failure group, revisit data 
## Stan code below will be problematic
if (nrow(survd %>% filter(outcome == 1))==0 || nrow(survd %>% filter(outcome == 0))==0) {
  stop('Regenerate simulated data - Stan code below will not work (yet) without some obs both censored & observed.')
}

## Fit standard cox-regression to these data
library(survival)
#aareg(formula = Surv(t, outcome) ~ init_observed_size, data = survd, nmin = 1)
coxph(formula = Surv(t, outcome) ~ rescaled_init_size, data = survd)

## imagine we had avg growth rate of these observations?
coxph(formula = Surv(t, outcome) ~ rescaled_init_size + rescaled_growth_rate, data = survd)

## test using simple survival model in stan 
## code from https://github.com/to-mi/stan-survival-shrinkage/blob/master/wei_hs.stan

## wei_gaus.stan
data_obs <- survd %>% filter(outcome == 1)
data_cens <- survd %>% filter(outcome == 0)
known_covars = c('rescaled_init_size')
biomarkers = c('rescaled_growth_rate')

M_bg = length(known_covars)
M_biom = length(biomarkers)

nchains <- 3
niter <- 1000
nwarmup <- niter/2

init1 <- list()
for (i in 1:nchains){
  init1[[i]] <- list(
    tau_s_bg_raw = 0.1*abs(rnorm(1)),
    tau_bg_raw = array(abs(rnorm(M_bg)), dim = M_bg),
    tau_s1_biom_raw = 0.1*abs(rnorm(1)),
    tau_s2_biom_raw = 0.1*abs(rnorm(1)),
    tau_biom_raw = array(abs(rnorm(M_biom)), dim = M_biom),
    tau1_biom_raw = array(abs(rnorm(M_biom)), dim = M_biom),
    tau2_biom_raw = array(abs(rnorm(M_biom)), dim = M_biom),
    alpha_raw = 0.01*rnorm(1),
    beta_bg_raw = array(rnorm(M_bg), dim = M_bg),
    beta_biom_raw = array(rnorm(M_biom), dim = M_biom),
    mu = rnorm(1)
  )
}

standata <- list(
  Nobs = nrow(data_obs)
  , Ncen = nrow(data_cens)
  , M_bg = M_bg
  , M_biom = M_biom
  , yobs = data_obs %>% dplyr::select(t) %>% unlist()
  , ycen = data_cens %>% dplyr::select(t) %>% unlist()
  , Xobs_bg = data_obs %>% dplyr::select(one_of(known_covars))
  , Xcen_bg = data_cens %>% dplyr::select(one_of(known_covars))
  , Xobs_biom = data_obs %>% dplyr::select(one_of(biomarkers))
  , Xcen_biom = data_cens %>% dplyr::select(one_of(biomarkers))
)
yobs = data_obs %>% dplyr::select(t) %>% unlist()
ycens = data_cens %>% dplyr::select(t) %>% unlist()

fit <- stan(file = "./stan-survival-shrinkage/wei_gau.stan"
            , data = standata
            , init = init1
            , chains = nchains, iter = niter
            , warmup = nwarmup
            , control = list(adapt_delta = 0.995)
            )



fit2 <- stan(file = "./stan-survival-shrinkage/wei_bg.stan"
            , data = standata
            , init = init1
            , chains = nchains, iter = niter
            , warmup = nwarmup
            , control = list(adapt_delta = 0.995)
)




standata_hs <- list(
  Nobs = nrow(data_obs)
  , Ncen = nrow(data_cens)
  , M_bg = M_bg
  , M_biom = M_biom
  , yobs = data_obs %>% dplyr::select(t) %>% unlist()
  , ycen = data_cens %>% dplyr::select(t) %>% unlist()
  , Xobs_bg = data_obs %>% dplyr::select(one_of(known_covars))
  , Xcen_bg = data_cens %>% dplyr::select(one_of(known_covars))
  , Xobs_biom = data_obs %>% dplyr::select(one_of(biomarkers))
  , Xcen_biom = data_cens %>% dplyr::select(one_of(biomarkers))
  , nu = 1
)
fit3 <- stan(file = "./stan-survival-shrinkage/wei_hs.stan"
             , data = standata_hs
             , init = init1
             , chains = nchains, iter = niter
             , warmup = nwarmup
             , control = list(adapt_delta = 0.995)
)

fit4 <- stan(file = "./stan-survival-shrinkage/wei_lap.stan"
             , data = standata
             , init = init1
             , chains = nchains, iter = niter
             , warmup = nwarmup
             , control = list(adapt_delta = 0.995)
)

