if (!require('bigrquery'))
  install.packages("bigrquery")
library(bigrquery)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---- get data ----

project <- "pici-1286"
biospec_sql <- 'SELECT * FROM [isb-cgc:tcga_201510_alpha.Biospecimen_data] where Study = "BLCA" LIMIT 1000'
bio_d <- query_exec(biospec_sql, project = project)

clin_sql <- 'SELECT * FROM [isb-cgc:tcga_201510_alpha.Clinical_data] where Study = "BLCA" LIMIT 1000'
clin_d <- query_exec(clin_sql, project = project)

## ---- explore biospec data ----

## Several biospec samples per person
table(duplicated(bio_d$ParticipantBarcode))

## Seems split between Tumor & normal samples
table(bio_d$SampleType)

## among tumor samples, only 3 Participants have more than one data point
bio_d %>%
  dplyr::filter(SampleType == 'Primary solid Tumor') %>%
  dplyr::group_by(ParticipantBarcode) %>% 
  dplyr::mutate(duped = ifelse(n()>1,1,0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(ParticipantBarcode) %>%
  dplyr::group_by(duped) %>%
  tally()

## those that are duplicated appear to have been collected at the same time 
## where one is ffpe & the other is not
bio_d %>% 
  group_by(SampleType,ParticipantBarcode) %>% 
  dplyr::mutate(duped = ifelse(n()>1,1,0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(SampleType == 'Primary solid Tumor' & duped == 1) %>%
  dplyr::select(ParticipantBarcode, SampleType, is_ffpe, days_to_collection, days_to_sample_procurement) %>%
  dplyr::arrange(ParticipantBarcode, is_ffpe)

## Among the rest, all are ffpe samples
bio_d %>%
  dplyr::filter(SampleType == 'Primary solid Tumor') %>%
  dplyr::group_by(ParticipantBarcode) %>% 
  dplyr::mutate(duped = ifelse(n()>1,1,0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(duped == 0) %>%
  dplyr::select(is_ffpe) %>%
  unlist() %>%
  table()

tumor_samples <- 
  bio_d %>%
  dplyr::filter(SampleType == 'Primary solid Tumor' & is_ffpe == 'NO')

stopifnot(!duplicated(tumor_samples$patid))

## ---- explore summary stats for biospec data ----

## what is the difference between avg, max & min for each metric?
tumor_wide <- 
  tumor_samples %>%
  tidyr::gather('variable','value', starts_with('avg'), starts_with('max'), starts_with('min')) %>%
  dplyr::mutate(stat = gsub(variable, pattern='^(.{3})_.*', replacement = "\\1")
                , variable = gsub(variable, pattern='^(.{3})_(.*)', replacement = "\\2")
                ) %>%
  tidyr::spread(stat, value)

## what is the typical range of values, across many patients?
tumor_wide %>%
  dplyr::mutate(range = max - min) %>%
  dplyr::select(range) %>%
  summary()

## how many participants have range > 0?
tumor_wide %>%
  dplyr::mutate(range = max - min
                , any_range = range > 0
                ) %>%
  dplyr::group_by(any_range, variable) %>%
  tally()

## (based on the above, appears that avg metric would likely be the preferred)

## plot values for each of the data points gathered
ggplot(data = tumor_wide, aes(x = ParticipantBarcode)) + 
  geom_point(aes(y = avg), size = 1) + 
  geom_segment(aes(y = min, yend = max, xend = ParticipantBarcode), alpha = 0.2) + 
  facet_wrap(~variable, scale = 'free_y')

## how would plot look different if data points were ordered by % normal?
participants <- 
  tumor_wide %>% 
  dplyr::filter(variable == 'percent_normal_cells') %>% 
  dplyr::rename(percent_normal = avg) %>% 
  dplyr::select(ParticipantBarcode, percent_normal) %>% 
  dplyr::distinct(ParticipantBarcode) %>%
  dplyr::arrange(percent_normal) %>%
  dplyr::mutate(participant_id = seq_len(n()))

tumor_wide2 <- tumor_wide %>%
  dplyr::left_join(participants
    , by = 'ParticipantBarcode') %>%
  dplyr::arrange(percent_normal) 

## plot values for each of the data points gathered
ggplot(data = tumor_wide2, aes(x = participant_id)) + 
  geom_point(aes(y = avg), size = 1) + 
  geom_segment(aes(y = min, yend = max, xend = participant_id), alpha = 0.2) + 
  facet_wrap(~variable, scale = 'free_y')

## which of these percentages add to 100?
tumor_samples %>% 
  dplyr::mutate(tot_cells = max_percent_normal_cells + max_percent_stromal_cells + max_percent_tumor_cells) %>% 
  dplyr::select(tot_cells) %>% 
  unlist() %>% 
  summary()

## look at samples where tot_cells < 100
tot_cells_lessthan_100 <- tumor_samples %>% 
  dplyr::mutate(tot_cells = avg_percent_normal_cells + avg_percent_stromal_cells + avg_percent_tumor_cells) %>% 
  dplyr::mutate(tot_cells_lessthan_100 = ifelse(tot_cells < 100, 'less than 100', '>= 100')) %>% 
  dplyr::select(ParticipantBarcode, tot_cells_lessthan_100, tot_cells) %>%
  dplyr::distinct(ParticipantBarcode)

tumor_wide3 <- tumor_wide2 %>%
  dplyr::left_join(tot_cells_lessthan_100, by = 'ParticipantBarcode')

## plot values for each of the data points gathered
ggplot(data = tumor_wide3, aes(x = participant_id, colour = tot_cells_lessthan_100)) + 
  geom_point(aes(y = avg), size = 1) + 
  geom_segment(aes(y = min, yend = max, xend = participant_id), alpha = 0.2) + 
  facet_wrap(~variable, scale = 'free_y')

## add percent_necrosis to total
tot_cells_lessthan_100_2 <- tumor_samples %>% 
  dplyr::mutate(tot_cells = avg_percent_normal_cells + avg_percent_stromal_cells + avg_percent_tumor_cells + avg_percent_necrosis) %>% 
  dplyr::mutate(tot_cells_lessthan_100 = ifelse(tot_cells < 100, 'less than 100', '>= 100')) %>% 
  dplyr::select(ParticipantBarcode, tot_cells_lessthan_100, tot_cells) %>%
  dplyr::distinct(ParticipantBarcode)

summary(tot_cells_lessthan_100_2$tot_cells)

## ---- explore relationship among tumor-derived biospec data ----

tumor_samples %>%
  dplyr::select(starts_with('avg'), days_to_sample_procurement) %>%
  pairs()

## from this, looks like some correlation between percent necrosis & percent tumor cells. 
## IE one is to the exclusion of the other at the high end of ~100% tumor cells.   

## nevertheless, the most "independent" data points are:
   # 1. monocyte/lymphocyte infiltrates (need to handle outliers)
   # 2. percent_necrosis
   # 3. percent tumor cells
   # 4. percent tumor nuclei

## ---- explore data for normal samples ----

normal_samples <-
  bio_d %>%
  dplyr::filter(SampleType != 'Primary solid Tumor')

## how many duplicates do we have?
table(duplicated(normal_samples$ParticipantBarcode))


## ---- prep data ----

d <- clin_d %>%
  mutate(mortality = ifelse(vital_status == 'Dead', 1, 0),
         surv_time = ifelse(days_to_last_known_alive <= 1, 1, days_to_last_known_alive),
         t_status = ifelse(grepl(pathologic_T, pattern = '^T2'), 2,
                           ifelse(grepl(pathologic_T, pattern = '^T3'), 3,
                                  ifelse(grepl(pathologic_T, pattern = '^T4'), 4,
                                         ifelse(pathologic_T %in% c('T0','T1'), 1,
                                                4)))),
         current_smoker = ifelse(tobacco_smoking_history == 'Current smoker', 1, 0),
         number_pack_years = ifelse(is.na(number_pack_years_smoked), 0, number_pack_years_smoked),
         reformed_smoker = ifelse(grepl(tobacco_smoking_history, pattern = 'Current reformed smoker', ignore.case = T), 1, 0)
  ) %>%
  dplyr::left_join(tumor_samples %>% 
                     dplyr::filter(days_to_sample_procurement <= 100) %>%
                     dplyr::select(ParticipantBarcode, starts_with('avg'))
                   , by = c('ParticipantBarcode')
                   )

## ---- look at lymphocyte infiltration & survival ---- 

## according to http://www.nature.com/modpathol/journal/v22/n2s/pdf/modpathol20091a.pdf
## lymphocyte invasion is prognostic indicator for overall survival
## in PT0 cancers. However, we don't have very many of these PT0 samples

## let's see if this is a useful indicator in our sample overall
xph <- coxph(Surv(surv_time, mortality) ~
               I(avg_percent_lymphocyte_infiltration>0),
             data=d
)
summary(xph)

## how is this indicator related to PT status?
ggplot(d, aes(x = t_status, y = log1p(avg_percent_lymphocyte_infiltration))) + geom_jitter()

## and, to N status?
ggplot(d, aes(x = pathologic_N, y = avg_percent_lymphocyte_infiltration)) + geom_jitter()

## and, to M status?
ggplot(d, aes(x = pathologic_M, y = avg_percent_lymphocyte_infiltration)) + geom_jitter()

## and, to smoking status?
ggplot(d, aes(x = current_smoker, y = avg_percent_lymphocyte_infiltration)) + geom_jitter()


## How does this fare in an analysis adjusted for other confounders?
xph <- coxph(Surv(surv_time, mortality) ~
               age_at_initial_pathologic_diagnosis + gender + 
               pathologic_M + t_status + pathologic_N + height + 
               current_smoker + reformed_smoker + 
               I(avg_percent_lymphocyte_infiltration>0),
             data=d
)
summary(xph)

## Most papers cite lymphocyte infiltration as a negative. Why would lymphocyte infiltration be a bad thing?
# According to : http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4846726/, 
#    "Lymphatics serve as a major pathway for metastatic spread in many types of cancers. "

## ---- look at the tumor necrosis - a known prognostic indicator  ---- 

## let's look at tumor necrosis (>10%) vis-a-vis LVI
xph <- coxph(Surv(surv_time, mortality) ~
               I(avg_percent_necrosis>10),
             data=d
)
summary(xph)

xph <- coxph(Surv(surv_time, mortality) ~
               I(avg_percent_necrosis>10) + 
               I(avg_percent_lymphocyte_infiltration>0),
             data=d
)
summary(xph)

xph <- coxph(Surv(surv_time, mortality) ~
               age_at_initial_pathologic_diagnosis + gender + 
               pathologic_M + t_status + pathologic_N + height + 
               current_smoker + reformed_smoker + 
               I(avg_percent_lymphocyte_infiltration>0) + 
               I(avg_percent_necrosis>10),
             data=d
)
summary(xph)



