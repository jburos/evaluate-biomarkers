
library(dplyr)

prep_data <- function(data) {
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
    dplyr::mutate(total = n()) %>%
    dplyr::group_by(failure_status, failure_or_progression_status, total) %>%
    dplyr::summarise(n = n()
              , percent = unique(n() / total)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-total) %>%
    dplyr::mutate(percent = paste0(round(percent*100, 0),'%')) %>%
    print()
  
  ## analysis data - survdata combined with original data
  adata <- data %>%
    dplyr::inner_join(survd 
               , by = "patid") %>%
    dplyr::filter(observed == 1) %>%
    dplyr::mutate(overall_mean_size = mean(observed_size, na.rm = T)
           , overall_sd_size = sd(observed_size, na.rm = T)
    ) %>%
    dplyr::group_by(patid) %>%
    dplyr::mutate(rescaled_tumor_size = (observed_size - overall_mean_size)/overall_sd_size
           , patient_mean_size = mean(observed_size)
           , patient_sd_size = sd(observed_size)
           , rescaled_patient_observed_size = (observed_size - patient_mean_size)/ patient_sd_size
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(patid) %>%
    dplyr::arrange(t) %>%
    dplyr::mutate(cum_progression = cumsum(progression)
                  , pr_progression = c(0, cum_progression[-n()]) ## take previous obs of cum_progression for current timepoint
    ) %>%
    dplyr::ungroup()
  
  return(list(per_observation = adata, per_patient = survd))
}
