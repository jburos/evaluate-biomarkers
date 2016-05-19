
#' Function to simulate data for simple survival model, as applied to onco-immunotherapy data
#' 
#' @description 
#' This posits a generic logistic growth model over time, where patient-level hazard is a function of tumor size. 
#' The function takes several functions as input, used to generate components of the simulation. The only fixed element is the 
#' growth over time, given two parameters : rate (varies over time), and most recent size of tumor. 
#' 
#' The rate of growth in tumor size (T) at any timepoint t, given current tumor size (base), max_size & current rate is:
#'  
#'  dTdt = (base*rate*(max_size-base)/max_size)
#' 
#' @param n sample size (number of units observed)
#' @param max_t max period of time per unit
#' @param max_size hypothesized max tolerable size of the tumor, at which growth is truncated.
#' @param plot (default TRUE) if true, plot the simulated data
#' 
#' @section parameters for tumor size over time
#' @param init_size_fun function yielding initial tumor size(s) per obs. Takes single param (n)
#' @param growth_rate_fun function yielding growth rate parameters per obs. Takes single param (n)
#' @param growth_rate_noise_fun function yielding variance in growth rate per obs*timepoint. Takes single param (n)
#' @param size_noise_fun function yielding measurement error in size of tumor, per obs*timepoint. Does not impact hazard, only observed values. Takes single param (n)
#' @param observed_size_fun function taking named vector of values for each obs*timepoint, returns observed tumor size.
#' 
#' @section parameters for estimated hazard
#' @param hazard_noise_fun function yielding hazard noise parameters. Takes single param (n)
#' @param hazard_coefs_fun function yielding nXc matrix of named coefs for each obs. Takes single param (n)
#' @param hazard_fun function yielding hazard estimate, given set of input params. Takes list of values 
#' 
#' @section parameters for censoring / behavior
#' @param censor_fun function yielding censor times. Takes single param (n)
#'
#' @import purrr
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom dplyr `%>%`
#'
#' @return data frame containing (long-version) of simulated data & parameters
#' 
#' @example 
#' simdt <- simulate_data()
#' 
simulate_data = function(
  n = 20
  , max_t = 50
  , max_size = 40000
  , prob_failure = 1/max_size
  , plot = TRUE
  , init_size_fun = create_rt(df = 5, ncp = 0, half = TRUE)
  , growth_rate_fun = create_rbeta(shape1 = 10, shape2 = 20)
  , growth_rate_noise_fun = create_rnorm(mean = 0, sd = 1)
  , size_noise_fun = create_rnorm(mean = 0, sd = 1)
  , observed_size_fun = function(row) {(vol2rad(row$tumor_size)*2 + row$size_noise)} ## diameter is observed, whereas volume determines hazard
  , hazard_noise_fun = create_rcauchy(location = 0, scale = 2, half = TRUE)
  , hazard_coefs_fun = function(n) { list(intercept = rnorm(n, mean = 0, sd = 1), beta_tumor_size = rnorm(n, mean = 3, sd = 1)) }
  , hazard_fun = function(row) { row$intercept + row$tumor_size*row$beta_tumor_size + row$hazard_noise }
  , censor_time_fun = create_scalar(value = max_t)   ## create_rt(df = 10, ncp = 20, half = TRUE)
  , failure_threshold = 3 ## >= this many events == FAILURE
  , progression_threshold = 2 ## >= this many events == PROGRESSION (disease progression)
  ) {

  ## simulate tumor growth over time at patient level
  simd <- 
    data.frame(patid = seq_len(n)) %>%
      mutate(growth_rate = growth_rate_fun(n = n)
             , init_size = init_size_fun(n = n)
             , censor_time = censor_time_fun(n = n)
             ) %>%
    bind_cols(hazard_coefs_fun(n = n))

  ## simulate data at patid*timepoint level
  simdt <- as.data.frame(expand.grid(patid = seq_len(n)
                                     , t = seq_len(max_t+1)-1
                                     )) %>%
    dplyr::mutate(hazard_noise = hazard_noise_fun(n = n())
                  , size_noise = size_noise_fun(n = n())
                  , growth_rate_noise = growth_rate_noise_fun(n = n())
                  ) %>%
    inner_join(simd, by = 'patid')

  simdt <- 
    simdt %>% 
    arrange(patid, t) %>%
    group_by(patid) %>%
    arrange(t) %>%
    dplyr::mutate(
      growth = ifelse(t == 0, init_size, ifelse(growth_rate + growth_rate_noise < 0, 0, growth_rate + growth_rate_noise))
      , tumor_size = purrr::accumulate(growth
                                , function(base, rate) {
                                    dTdt <- (base*rate*(max_size-base)/max_size)
                                    #print(paste0('base: ',round(base),'; rate: ',round(rate),'; dTdt'))
                                    base + dTdt
                                  }
                                )
      ) %>%
    ungroup() %>%
    dplyr::filter(t > 0) %>%
    mutate(observed = ifelse(t >= censor_time, 0, 1))
  
  simdt$observed_size <- 
    simdt %>%
    rowwise() %>%
    do(observed_size = observed_size_fun(.)) %>% 
    unlist()
  
  
  simdt$hazard <- 
    simdt %>%
    rowwise() %>%
    do(hazard = hazard_fun(.)) %>% 
    unlist()
  
  # simulate failure process
  simdt2 <- 
    simdt %>%
    rowwise() %>% 
    ## calc prob of failure (rowwise b/c each obs has different hazard value)
    dplyr::mutate(eff_hazard = round(ifelse(hazard >= 4000, 4000, ifelse(hazard < 0, 0, hazard)), digits = 0)
                  , events = rbinom(n = n(), size = eff_hazard, prob = prob_failure)
                  ) %>% 
    ungroup() %>% 
    group_by(patid) %>% 
    mutate(
      failure = ifelse(events >= failure_threshold, 1, 0)
      , failure_status = max(failure)
      , first_failure = min(ifelse(failure == 1, t, max_t + 1), na.rm = T) 
      , progression = ifelse(events >= progression_threshold, 1, 0)
      , progression_status = max(progression)
      , first_progression = min(ifelse(progression == 1, t, max_t + 1), na.rm = T)
      , failure_or_progression = ifelse(failure == 1, 1, ifelse(progression == 1, 1, 0))
      , failure_or_progression_status = max(failure_or_progression)
      , first_failure_or_failure = min(ifelse(failure_or_progression == 1, t, max_t + 1), na.rm = T)
    ) %>%
    ## marked post-failure events as unobserved
    dplyr::mutate(observed = ifelse(t > first_failure, 0, observed)) %>%
    ungroup()

  simdt <- simdt2
  rm(simdt2)

  simdt
}


#' Plot simulated data 
#' 
#' @param d simulated data, generated by simulate_data()
#' @param n (optional) number of simulated obs to plot. defaults to ALL
#' 
#' @returns ggplot grob
#' @import ggplot2
#' 
#' @example
#' simdt <- simulate_data()
#' plot_simualated_data(simdt)
#' 
plot_simulated_data <- function(d, n = NULL) {
  ## filter input data to restrict to N products
  ## sampled at random
  if (!is.null(n)) {
    d = d %>% 
      semi_join(d %>% 
                  dplyr::distinct(patid) %>%
                  dplyr::select(patid) %>%
                  dplyr::ungroup() %>%
                  dplyr::sample_n(n) 
                , by = 'patid'
                )
  }
  
  ## plot time-series data 
  pl1 <- 
    ggplot2::ggplot() + 
    ggplot2::geom_line(
      data = 
        d %>% 
        dplyr::filter(observed == 1) %>%
        tidyr::gather(var, value, hazard, tumor_size) %>% 
        dplyr::mutate(grt = paste('growth_rate:',round(growth_rate, digits = 1),sep=' ')
                      , init = paste('init_size:',round(init_size, digits = 1), sep=' ')
        )
      , mapping = ggplot2::aes(x = t, y = value, colour = var, group = var)) + 
    ggplot2::facet_wrap(~patid + grt + init, scale = 'free_y') +
    ggplot2::scale_colour_discrete("Measurement type")
    
  ## add censor times
  pl2 <- 
    pl1 +
    ggplot2::geom_vline(data = 
                 d %>% 
                 dplyr::filter(observed == 1) %>% 
                 dplyr::group_by(patid) %>% 
                 dplyr::filter(t == max(t)) %>% 
                 dplyr::ungroup() %>% 
                 dplyr::mutate(
                   failure_type = ifelse(failure == 1,'failure','censor')
                   , grt = paste('growth_rate: ',round(growth_rate, digits = 1),sep='')
                   , init = paste('init_size: ',round(init_size, digits = 1), sep='')
                 )
               , ggplot2::aes(xintercept = t, linetype = failure_type), colour = 'black') +
    ggplot2::scale_linetype_manual('Survival status', values = c('failure' = 'solid', 'censor' = 'dashed'))
  
  ## add "failure" (progression) events 
  pl <- 
    pl2 + 
    ggplot2::geom_point(
      data = 
        d %>%
        dplyr::filter(observed == 1 & progression == 1) %>%
        dplyr::mutate(
          failure_type = 'progression'
          , grt = paste('growth_rate: ',round(growth_rate, digits = 1),sep='')
          , init = paste('init_size: ',round(init_size, digits = 1), sep='')
        )
      , ggplot2::aes(x = t, y = hazard, shape = failure_type), colour = 'black', size = 1)
  pl
}

#' helper functional wrapping 'rep', intended for use with simulate_data
#'
#' @param val
#' 
#' @returns function taking parameter 'n' that repeats value n times
#' 
create_scalar <- function(value) {
  function(n) {
    rep(x = value, times = n)
  }
}

#' helper functional wrapping 'rnorm', inteded for use with simulate_data
#' 
#' @param mean
#' @param sd
#' 
#' @import purrr
#' 
#' @returns function taking parameter 'n' returning draws from normal distribution
#' 
create_rnorm <- function(mean, sd) {
  purrr::partial(rnorm, mean = mean, sd = sd)
} 

#' helper function to truncate a distribution. Intended for use with simulate_data
#' 
#' @param .val minimum value - draws are filtered to be greater than this value
#' @param .dist function yielding samples to be filtered
#' @param n number of obs - required parameter to .dist
#' @param ... params to .dist
#' 
#' @import purrr
#' 
#' @returns output from .dist (vector of length n), filtered so obs >= .val
#' 
left_truncate <- function(.dist, .val, n = n, ...) {
  .dist(n = n*10, ...) %>%
    purrr::keep(~ .x >= .val) %>%
    sample(., size = n, replace = TRUE)
}
  
#' helper functional for rt, intended for use with simulate_data
#' 
#' @param df degrees of freedom for t distribution
#' @param ncp (optional) param to rt. See rt for details
#' @param half (default FALSE) if true, truncates result to x > 0 
#' 
#' @import purrr
#' 
#' @returns function taking parameter 'n' returning draws from t distribution
create_rt <- function(df, half = FALSE, ncp = 0) {
  if (half == TRUE)
    purrr::partial(left_truncate, .dist = rt, .val = ncp, df = df, ncp = ncp)
  else {
    purrr::partial(rt, df = df, ncp = ncp)
  }
} 

#' helper functional for rcauchy, intended for use with simulate_data
#' 
#' @param location
#' @param scale
#' @param half (default FALSE) if TRUE, result is truncated at values > location
#' 
#' @import purrr
#' 
#' @returns function taking parameter 'n' returning draws from cauchy distribution
create_rcauchy <- function(location, scale, half = FALSE) {
  if (half == TRUE)
    purrr::partial(left_truncate, .dist = rcauchy, .val = location, location = location, scale = scale)
  else {
    purrr::partial(rcauchy, location = location, scale = scale)
  }
}

#' helper functional for rbeta, intended for use with simulate_data
#' 
#' @param shape1
#' @param shape2
#' @param half if TRUE, result is truncated at values > location
#' 
#' @import purrr
#' 
#' @returns function taking parameter 'n' returning draws from beta distribution
create_rbeta <- function(shape1, shape2) {
  purrr::partial(rbeta, shape1 = shape1, shape2 = shape2)
}

## v = (4/3) * pi * r^3
vol2rad <- function(volumes) {
  (volumes * (3/4) / pi)^(1/3)
}
