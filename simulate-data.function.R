
#' Function to simulate data for simple survival model (applied to onco-immunotherapy data)
#' 
#' 
#' @param n sample size (number of units observed)
#' @param max_t max period of time per unit
#' @param growth_rate_fun function yielding growth rate parameters. Takes single param (n)
#' @param hazard_noise_fun function yielding hazard noise parameters. Takes single param (n)
#' @param hazard_coefs_fun function yielding nXc matrix of named coefs for each obs. Takes single param (n)
#' @param hazard_fun function yielding hazard estimate, given set of input params. Takes list of values 
#' @param censor_fun function yielding censor times. Takes single param (n)
#' @param plot (default TRUE) if true, plot the simulated data
#'
#' @import purrr
#' @import dplyr
#' @import ggplot2
#'
#' @return data frame containing (long-version) of simulated data & parameters
#' 
simulate_data = function(
  n = 20
  , max_t = 50
  , max_size = 2000
  , init_size_fun = create_rt(df = 5, ncp = 0, half = TRUE)
  , growth_rate_fun = create_rbeta(10, 20)
  , growth_rate_noise_fun = create_rnorm(mean = 0, sd = 1)
  , size_noise_fun = create_rnorm(mean = 0, sd = 1)
  , observed_size_fun = function(row) {(row$tumor_size + row$size_noise)}
  , hazard_noise_fun = create_rcauchy(location = 0, scale = 5, half = TRUE)
  , hazard_coefs_fun = function(n) { list(intercept = rnorm(n, mean = 0, sd = 1), beta_tumor_size = rnorm(n, mean = 3, sd = 1)) }
  , hazard_fun = function(row) { row$intercept + row$tumor_size*row$beta_tumor_size + row$hazard_noise }
  , censor_time_fun = create_rt(df = 10, ncp = 10, half = TRUE)
  , plot = TRUE
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
      , tumor_size = accumulate(growth
                                , function(base, rate) {
                                    dTdt <- (base*rate*(max_size-base)/max_size)
                                    #print(paste0('base: ',round(base),'; rate: ',round(rate),'; dTdt'))
                                    base + dTdt
                                  }
                                )
      ) %>%
    ungroup() %>%
    dplyr::filter(t > 0)
  
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
  
  ## review simulated data 
  if (plot == TRUE) {
    ggplot(simdt %>% 
             gather(var, value, hazard, tumor_size) %>% 
             mutate(grt = paste('growth_rate:',round(growth_rate, digits = 1),sep='')
                    , init = paste('init_size:',round(init_size, digits = 1), sep='')
                    )
           , aes(x = t, y = value, colour = var, group = var)) + 
      geom_line() + 
      facet_wrap(~patid + grt + init, scale = 'free_y')
  }
  
  simdt
}

#' helper functional wrapping 'rep'
#'
#' @param val
#' 
#' @returns function taking parameter 'n' that repeats value n times
#' 
create_scalar <- function(value) {
  function(.n) {
    rep(x = value, times = .n)
  }
}

#' helper functional wrapping 'rnorm'
#' 
#' @param mean
#' @param sd
#' 
#' @returns function taking parameter 'n'
create_rnorm <- function(mean, sd) {
  purrr::partial(rnorm, mean = mean, sd = sd)
} 

#' helper function to truncate a dist 
#' 
#' @param .val value 
#' @param .dist
#' @param ... params to .dist
#' 
#' @returns 
left_truncate <- function(.dist, .val, n = n, ...) {
  .dist(n = n*10, ...) %>%
    purrr::keep(~ .x >= .val) %>%
    head(., n = n)
}
  
#' helper functional for rt
#' 
#' @param df degrees of freedom for t distribution
#' @param ncp (optional) param to rt. See rt for details
#' @param half (default FALSE) if true, truncates result to x > 0 
#' 
#' @returns function taking parameter 'n'
create_rt <- function(df, half = FALSE, ncp = NULL) {
  if (half == TRUE)
    purrr::partial(left_truncate, .dist = rt, .val = ncp, df = df, ncp = ncp)
  else {
    purrr::partial(rt, df = df, ncp = ncp)
  }
} 

#' helper functional for rcauchy
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

#' helper functional for rbeta
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