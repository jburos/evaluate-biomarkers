rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstan)
library(survival)
library(arm)
library(lme4)
rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(),3))

## ------ prep simulated data -------

source('simulate-data.function.R')

## get simulated data 
data <- simulate_data(n = 100, max_size = 4000, max_t = 50, failure_threshold = 4, progression_threshold = 3)

## look at distribution of "adverse" events in data
## >= 3 events count as "progression" events.
## >= 4 events count as "failure" (mortality)
(
  pl <- ggplot(
    data %>% 
      mutate(frequency = 1
             , event_type = factor(ifelse(failure == 1, 2
                                          , ifelse(progression == 1, 1
                                                   , 0))
                                   , levels = c(0, 1, 2)
                                   , labels = c('No event', 'Progression', 'Failure')
                                   )
             )
    , aes(y = frequency, x = events, fill = event_type)) + 
  stat_summary(geom = 'bar', fun.y = sum) +
  scale_x_continuous('adverse "events"', breaks = c(0,seq_len(max(data$events)))) +
  ggtitle("Event thresholds: \n Event type according to simulated severity")
)
if (!interactive()) {
  ggsave(pl, file = 'simulated_event_freq_showing_thresholds.png')
}

## example trajectories for a few products
(pl <- plot_simulated_data(data, n = 12) +
  ggtitle('Simulated data for 12 sample patients'))
if (!interactive()) {
  ggsave(pl, file = 'simulated_data_example_patients.png')
}

## distribution of survival times by failure status
(
  pl <- ggplot(data %>% 
                 distinct(patid) %>% 
                 dplyr::mutate(failure_status = factor(failure_status, levels = c(0,1), labels = c('censor', 'failure'))) %>%
                 ungroup()
               , aes(x = first_failure, group = failure_status, fill = failure_status, colour = failure_status)) +
  geom_histogram(alpha = 0.2, position = 'dodge')
)
if (!interactive()) {
  ggsave(pl, file = 'simulated_data_survival_time_distribution.png')
}

## prep data for survival analysis (typical approach)
survd <- data %>%
  dplyr::filter(observed == 1) %>%
  dplyr::group_by(patid) %>%
  dplyr::mutate(
    surv_time = max(t, na.rm = T)
    , init_observed_size = first(observed_size, order_by = t)
    , last_observed_size = last(observed_size, order_by = t)
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
  dplyr::select(patid
                , ends_with('status') ## patient-level survival/censor status
                , starts_with('first') ## patient-level time to survival/censor events
                , init_observed_size
                , last_observed_size
                , avg_growth_rate
                , starts_with('rescaled')
                , starts_with('sd')
                , starts_with('mean')
                )

## display event rate(s)
survd %>%
  mutate(total = n()) %>%
  group_by(failure_status, failure_or_progression_status, total) %>%
  summarise(n = n()
         , percent = unique(n() / total)
         ) %>%
  ungroup() %>%
  dplyr::select(-total) %>%
  mutate(percent = paste0(round(percent*100, 0),'%'))

## analysis data - survdata combined with original data
adata <- data %>%
  inner_join(survd 
             , by = "patid") %>%
  filter(observed == 1) %>%
  mutate(overall_mean_size = mean(observed_size, na.rm = T)
         , overall_sd_size = sd(observed_size, na.rm = T)
         ) %>%
  group_by(patid) %>%
  mutate(rescaled_tumor_size = (observed_size - overall_mean_size)/overall_sd_size
         , patient_mean_size = mean(observed_size)
         , patient_sd_size = sd(observed_size)
         , rescaled_patient_observed_size = (observed_size - patient_mean_size)/ patient_sd_size
  ) %>%
  ungroup() %>%
  dplyr::group_by(patid) %>%
  dplyr::arrange(t) %>%
  dplyr::mutate(cum_progression = cumsum(progression)
                , pr_progression = c(0, cum_progression[-n()]) ## take previous obs of cum_progression for current timepoint
  ) %>%
  ungroup()

## if only one obs in censor or failure group, revisit data 
## Stan code below will be problematic
if (nrow(survd %>% filter(outcome == 1))==0 || nrow(survd %>% filter(outcome == 0))==0) {
  stop('Regenerate simulated data - Stan code below will not work (yet) without some obs both censored & observed.')
}

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



standata <- list(
  N = nrow(newdata)
  , S = max(newdata$patid)
  , T = max(newdata$t)
  , X = 1
  , s = newdata$patid
  , t = newdata$t
  , event = newdata$failure_event
  , covars = newdata %>% dplyr::select(rescaled_init_size)
)

testfit <- stan('long_surv.stan', data = standata, chains = 1, iter = 10)

stanfit <- stan('long_surv.stan', data = standata, chains = 3, iter = 1000)

print(stanfit, 'beta')

## ------ semi-competing-risks survival model using stan -------

## look at joint model for disease progression & survival 
## called a semi-competing risks model since progression not possible after mortality
## using parameterization here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4744123/ 
## this allows for post-progression risk to be different from pre-progression risk

  
standata2 <- list(
  N = nrow(newdata2)
  , S = max(newdata2$patid)
  , T = max(newdata2$t)
  , X = 1
  , s = newdata2$patid
  , t = newdata2$t
  , event = newdata2$failure_event
  , covars = newdata2 %>% dplyr::select(rescaled_init_size)
  , progression = newdata2$progression
  , pr_progress = newdata2$pr_progression
)

testfit2 <- stan('long_surv_progression.stan'
                 , model_name = 'semi-competing risks model with 3 submodels & constrained coefs'
                 , data = standata2
                 , chains = 1
                 , iter = 10
                 )
stanfit2 <- stan('long_surv_progression.stan'
                 , data = standata2
                 , chains = 3
                 , iter = 1000
                 )

print(stanfit2, c('beta_shared'))

## variation allowing betas to vary across submodels
testfit2a <- stan('long_surv_progression_free_beta.stan'
                 , model_name = 'semi-competing risks model with 3 submodels & unconstrained coefs'
                 , data = standata2
                 , chains = 1
                 , iter = 10
)
stanfit2a <- stan('long_surv_progression_free_beta.stan'
                 , data = standata2
                 , chains = 3
                 , iter = 1000
)

print(stanfit2a, c('beta_shared', 'beta1', 'beta2', 'beta3'))

## variation of model using only 2 submodels.


standata3 <- list(
  N = nrow(newdata2)
  , S = max(newdata2$patid)
  , T = max(newdata2$t)
  , X = 1
  , s = newdata2$patid
  , t = newdata2$t
  , event = newdata2$failure_event
  , covars = newdata2 %>% dplyr::select(rescaled_init_size)
  , progression = newdata2$progression
  , pr_progress = newdata2$pr_progression
)

testfit3 <- stan('long_surv_progression2.stan'
                 , model_name = 'semi-competing risks model with 2 submodels'
                 , data = standata3
                 , chains = 1
                 , iter = 10
                 )
stanfit3 <- stan('long_surv_progression2.stan'
                 , data = standata3
                 , chains = 3
                 , iter = 1000
                 )

print(stanfit3, c('beta'))

## variation of model allowing betas to be different

standata4 <- list(
  N = nrow(newdata2)
  , S = max(newdata2$patid)
  , T = max(newdata2$t)
  , X = 1
  , s = newdata2$patid
  , t = newdata2$t
  , event = newdata2$failure_event
  , covars = newdata2 %>% dplyr::select(rescaled_init_size)
  , progression = newdata2$progression
  , pr_progress = newdata2$pr_progression
)

testfit4 <- stan('long_surv_progression2_free_beta.stan'
                 , model_name = 'semi-competing risks model with 2 submodels & unconstrained betas'
                 , data = standata4
                 , chains = 1
                 , iter = 10
                 )
stanfit4 <- stan('long_surv_progression2_free_beta.stan'
                 , data = standata4
                 , chains = 3
                 , iter = 1000
                 )

print(stanfit4, c('beta1','beta2'))

