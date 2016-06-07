library(dplyr)
library(ggplot2)
library(survival)

## ------ prep & review simulated data -------

source('simulate-data.function.R')
source('prep-data.function.R')

simulate_survival_analysis <- function(n = 100, max_size = 4000, max_t = 50, failure_threshold = 4, progression_threshold = 3) {
  ## get simulated data 
  data <- simulate_data(n = n, max_size = max_size, max_t = max_t, failure_threshold = failure_threshold, progression_threshold = progression_threshold)

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
  ss <- summary(survfit)
  res <- list(method = 'coxph', outcome = 'failure', coef = 'rescaled_init_size', estimate = ss$conf.int['rescaled_init_size',], pval = ss$logtest['pvalue'], concordance = ss$concordance[1])
  rm(ss)

  ## ------ review survival model (pfs) using coxph -------

  survfit2 <- coxph(
    formula = Surv(first_failure_or_progression, failure_or_progression_status) ~ rescaled_init_size
    , data = survd
    )
 
  ss <- summary(survfit2)
  res2 <- list(method = 'coxph', outcome = 'failure_or_progression', coef = 'rescaled_init_size', estimate = ss$conf.int['rescaled_init_size',], pval = ss$logtest['pvalue'], concordance = ss$concordance[1])
  return(list(res, res2))
}

