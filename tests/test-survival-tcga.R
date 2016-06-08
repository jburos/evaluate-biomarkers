install.packages("bigrquery")
library(bigrquery)

## --- get data ---- 
project <- "pici-1286"
sql <- 'SELECT * FROM [isb-cgc:tcga_201510_alpha.Clinical_data] where Study = "BLCA" LIMIT 1000'
d <- query_exec(sql, project = project)

## --- prep data ---
d <- d %>%
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
         )


## -- standard survival model -- 
sr <- survreg(Surv(surv_time, mortality) ~
                age_at_initial_pathologic_diagnosis + gender + 
                pathologic_M + t_status + pathologic_N + height + 
                current_smoker + reformed_smoker,
              data=d,
              dist="weibull"
              )
summary(sr)


## -- cox-ph model -- 
xph <- coxph(Surv(surv_time, mortality) ~
                age_at_initial_pathologic_diagnosis + gender + 
                pathologic_M + t_status + pathologic_N + height + 
                current_smoker + reformed_smoker,
              data=d
            )
summary(xph)



