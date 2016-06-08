
library(dplyr)
library(tidyr)

#' Takes a data frame of survival data (one obs per patient), reformats for stan model
#' 
#' @param data input data frame, prepped according to prep_data()$per_patient
#' @param time_precision decimal points to use when rounding time (e.g. 0 = nearest integer)
#' 
#' @returns a data frame
format_data_semicomprisks <- function(data, time_precision = 2, event1 = 'progression', event2 = 'failure', allow_unequal_counts = FALSE) {
  
  subjd <- data %>%
    dplyr::mutate(subject_id = row_number())
  
  time_vars <- stringr::str_c('first', c(event1, event2), sep = '_')
  event_vars <- stringr::str_c(c(event1, event2), 'status', sep = '_')
  
  
  ## rearrange data so that progression & failure events are on separate rows
  d <- 
    subjd %>%
    rename_(.dots = time_vars %>% set_names(c('time1','time2'))) %>%
    rename_(.dots = event_vars %>% set_names(c('event1','event2'))) %>%
    tidyr::gather(var, val
                  , starts_with('time')
                  , starts_with('event')
                  ) %>%
    dplyr::mutate(event_type = gsub(var,pattern=".*(\\d)$",replacement="\\1")
                  , var = gsub(var,pattern="(.*)(\\d)$",replacement="\\1")
    ) %>%
    tidyr::spread(var, val) %>%
    dplyr::rename(outcome = event)
  
  ## round failure times to nearest X 
  d <- 
    d %>%
    dplyr::mutate(time = round(time, time_precision)) 
  
  ## identify unique failure timepoints
  d <- d %>%
    dplyr::mutate(timepoint_id = as.integer(factor(time, ordered = TRUE)))
  
  tps <- d %>%
    dplyr::distinct(timepoint_id) %>%
    dplyr::select(timepoint_id, time) %>%
    dplyr::arrange(timepoint_id) %>%
    dplyr::rename(time_to = time) %>%
    dplyr::mutate(time_from = dplyr::lag(time_to, n = 1, order_by = timepoint_id)
                  , duration = ifelse(timepoint_id == 1, 10^(-1*time_precision), time_to - time_from)
    )
  
  ## could be made *MUCH* more efficient -- don't need cum/functions
  longd <- expand.grid(list(subject_id = seq_len(nrow(data)), timepoint_id = seq_len(nrow(tps)))) %>%
    as.data.frame() %>%
    dplyr::inner_join(d %>% dplyr::rename(event_timepoint = timepoint_id), by = 'subject_id') %>%
    dplyr::mutate(progression = ifelse(event_type == 1 & outcome == 1 & timepoint_id == event_timepoint, 1, 0)
                  , failure = ifelse(event_type == 2 & outcome == 1 & timepoint_id == event_timepoint, 1, 0)
                  , censor = ifelse(event_type == 2 & outcome == 0 & timepoint_id == event_timepoint, 1, 0)
    ) %>%
    dplyr::group_by(subject_id, timepoint_id) %>%
    dplyr::mutate(progression = max(progression)
                  , failure = max(failure)
                  , censor = max(censor)
    ) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::arrange(timepoint_id) %>%
    dplyr::mutate(post_progression = cummax(progression) - progression
                  , post_failure = cummax(failure) - failure
                  , post_censor = cummax(censor) - censor
    ) %>%
    ungroup() %>%
    dplyr::filter(post_failure == 0 & post_censor == 0) %>%
    dplyr::select(-event_timepoint, -time, -outcome, -event_type) %>%
    unique()
  
  ## confirm no duplicates by subject/timepoint
  stopifnot(all(duplicated(cbind(longd$subject_id, longd$timepoint_id))==FALSE))
  
  ## confirm number of eventsof type1 the same as in original data 
  for (event in c(event1, event2)) {
    original_number_events <- data %>% dplyr::select(one_of(stringr::str_c(event,'status',sep='_'))) %>% unlist() %>% sum()
    revised_number_events <- longd %>% dplyr::select(one_of(event)) %>% unlist() %>% sum()
    if (original_number_events != revised_number_events) {
      if (allow_unequal_counts == FALSE) {
        print(stringr::str_c('Note: number of ',event,' events has changed from ',original_number_events,' to ',revised_number_events))
        print("Here is a sample of data for the records that have changed:")
        disjoint_event_records <- 
          longd %>% 
          dplyr::group_by(patid) %>% 
          dplyr::summarise_each(funs = funs(new_status = 'sum'), one_of(event)) %>% 
          ungroup() %>% 
          inner_join(data, by = 'patid') %>% 
          dplyr::filter_(stringr::str_c(event,' != ',event,'_status'))
        print(head(disjoint_event_records))
        stopifnot(original_number_events == revised_number_events)
      } else {
        warning(stringr::str_c('Number of ',event,' events has changed from ',original_number_events,' to ',revised_number_events))
      }
    }
  }
  
  longd %>% inner_join(tps, by = 'timepoint_id')
}
