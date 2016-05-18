rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(),3))

## ------ prep simulated data -------

source('simulate-data.function.R')

## get simulated data 
data <- simulate_data(n = 100, max_size = 4000, max_t = 50)
pl <- plot_simulated_data(data, n = 12, progress_threshold = 2)
ggsave(pl, file = 'simulated_data.png')

## distribution of survival times
pl <- ggplot(data %>% 
         dplyr::filter(observed == 1) %>%
         group_by(patid) %>% 
         dplyr::summarise(surv_time  = max(t, na.rm = T)
                          , event = factor(last(failure_event, order_by = t), levels = c(0,1), labels = c('censor', 'failure'))
                          ) %>%
         ungroup()
       , aes(x = surv_time, group = event, fill = event, colour = event)) +
  geom_histogram(alpha = 0.2) +
  facet_wrap(~event)
ggsave(pl, file = 'distribution of event times.png')

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
  dplyr::select(patid, t, init_observed_size, last_observed_size, avg_growth_rate, outcome, starts_with('rescaled'))  %>%
  dplyr::rename(last_t = t)

survd %>%
  group_by(outcome) %>%
  tally()

## create joint growth & survival model 
newdata <- data %>%
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
  ungroup()
  

## if only one obs in censor or failure group, revisit data 
## Stan code below will be problematic
if (nrow(survd %>% filter(outcome == 1))==0 || nrow(survd %>% filter(outcome == 0))==0) {
  stop('Regenerate simulated data - Stan code below will not work (yet) without some obs both censored & observed.')
}

## ------ try fitting using JM (skipped) -------

## Try to fit model using JM package, as described in vignette
## note: doesn't work. 
#
# See the following:
#  Error in optim(thetas, opt.survWB, gr.survWB, method = "BFGS", control = list(maxit = if (it <  : 
#   non-finite value supplied by optim
if (TRUE == FALSE) {
  library(JM)
  library(nlme)
  # linear mixed model fit without treatment effect
  fitLME.null <- lme(rescaled_patient_observed_size ~ t,
                     random = ~ 1 | patid, data = newdata %>% dplyr::filter(t < last_t))
  
  # cox model fit without treatment effect
  fitCOX.null <- coxph(Surv(last_t, outcome) ~ 1,
                       data = survd, x = TRUE)
  
  # joint model fit without treatment effect
  fitJOINT.null <- jointModel(fitLME.null, fitCOX.null,
                              timeVar = "t")
  
  # linear mixed model fit with treatment effect
  fitLME.alt <- lme(rescaled_patient_observed_size ~ t*rescaled_init_size,
                    random = ~ 1 | patid, data = newdata %>% dplyr::filter(t < last_t))
  # cox model fit with treatment effect
  fitCOX.alt <- coxph(Surv(last_t, outcome) ~ rescaled_init_size,
                      data = survd, x = TRUE)
  # joint model fit with treatment effect
  fitJOINT.alt <- jointModel(fitLME.alt, fitCOX.alt, timeVar = "t",
                             method = "weibull-PH-aGH")
  # likelihood ratio test for treatment effect
  anova(fitJOINT.null, fitJOINT.alt)
} 

## ------ test survival model using coxph -------

## for sanity, start with separate surv + longitudinal models.

library(survival)
#aareg(formula = Surv(t, outcome) ~ init_observed_size, data = survd, nmin = 1)
survfit <- coxph(formula = Surv(last_t, outcome) ~ rescaled_init_size, data = survd)

## ------ test growth model using lmer -------

library(lme4)
growthfit <- lmer(rescaled_patient_observed_size ~ t + rescaled_init_size + 
                    (1 + rescaled_init_size + t | patid)
                  , data = newdata
                  )
                  

## ------ test survival model (long version) using stan -------

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

newdata2 <- newdata %>%
  dplyr::group_by(patid) %>%
  dplyr::arrange(t) %>%
  dplyr::mutate(progression = ifelse(failure > 0, 1, 0)
                , cum_progression = cumsum(progression)
                , pr_progression = c(0, cum_progression[-n()]) ## take previous obs of cum_progression for current timepoint
                ) %>%
  ungroup()
  
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

testfit2 <- stan('long_surv_progression.stan', data = standata2, chains = 1, iter = 10)
stanfit2 <- stan('long_surv.stan', data = standata2, chains = 3, iter = 1000)

print(stanfit2, 'beta')

