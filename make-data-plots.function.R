library(ggplot2)
library(dplyr)
source('plot-simulated-data.function.R')
make_data_plots <- function(data) {
  ## look at distribution of "adverse" events in data
  ## >= 3 events count as "progression" events.
  ## >= 4 events count as "failure" (mortality)
  print(
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
      scale_x_continuous('Simulated "events"', breaks = c(0,seq_len(max(data$events)))) +
      scale_fill_discrete('Clinical events') +
      ggtitle("Event thresholds: \n Event type according to simulated severity")
  )
  if (!interactive()) {
    ggsave(pl, file = 'simulated_event_freq_showing_thresholds.png')
  }
  
  ## example trajectories for a few products
  print(pl <- plot_simulated_data(data, n = 12) +
          ggtitle('Simulated data for 12 sample patients'))
  if (!interactive()) {
    ggsave(pl, file = 'simulated_data_example_patients.png')
  }
  
  ## distribution of survival times by failure status
  print(
    pl <- ggplot(data %>% 
                   distinct(patid) %>% 
                   dplyr::mutate(failure_status = factor(failure_status, levels = c(0,1), labels = c('censor', 'failure'))) %>%
                   ungroup()
                 , aes(x = first_failure, group = failure_status, fill = failure_status, colour = failure_status)) +
      geom_histogram(alpha = 0.2, position = 'dodge') + 
      ggtitle('Distribution of censor times for simulated data')
  )
  if (!interactive()) {
    ggsave(pl, file = 'simulated_data_survival_time_distribution.png')
  }
  
}