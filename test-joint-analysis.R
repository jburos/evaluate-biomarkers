rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(),3))

## ------ prep simulated data -------

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

## TODO try SemiCompRisks (https://cran.r-project.org/web/packages/SemiCompRisks/SemiCompRisks.pdf)

## ------ test survival model using coxph -------

## for sanity, start with separate surv + longitudinal models.

library(survival)
#aareg(formula = Surv(t, outcome) ~ init_observed_size, data = survd, nmin = 1)
survfit <- coxph(formula = Surv(first_failure, failure_status) ~ rescaled_init_size, data = survd)
print(survfit)

## ------ test growth model using lmer -------

library(lme4)
library(arm)
growthfit <- lmer(rescaled_patient_observed_size ~ t + rescaled_init_size + 
                    (1 + t || patid)
                  , data = adata
                  )
display(growthfit)        

## ------ test survival model (long version) using stan -------

standata <- list(
  N = nrow(adata)
  , S = max(adata$patid)
  , T = max(adata$t)
  , X = 1
  , s = adata$patid
  , t = adata$t
  , event = adata$failure_status
  , covars = adata %>% dplyr::select(rescaled_init_size)
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
  N = nrow(adata)
  , S = max(adata$patid)
  , T = max(adata$t)
  , X = 1
  , s = adata$patid
  , t = adata$t
  , event = adata$failure_status
  , covars = adata %>% dplyr::select(rescaled_init_size)
  , progression = adata$progression
  , pr_progress = adata$pr_progression
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
  N = nrow(adata)
  , S = max(adata$patid)
  , T = max(adata$t)
  , X = 1
  , s = adata$patid
  , t = adata$t
  , event = adata$failure_status
  , covars = adata %>% dplyr::select(rescaled_init_size)
  , progression = adata$progression
  , pr_progress = adata$pr_progression
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
  N = nrow(adata)
  , S = max(adata$patid)
  , T = max(adata$t)
  , X = 1
  , s = adata$patid
  , t = adata$t
  , event = adata$failure_status
  , covars = adata %>% dplyr::select(rescaled_init_size)
  , progression = adata$progression
  , pr_progress = adata$pr_progression
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

