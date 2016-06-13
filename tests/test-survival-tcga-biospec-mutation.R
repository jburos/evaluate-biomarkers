if (!require('bigrquery'))
  install.packages("bigrquery")
library(bigrquery)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
### don't load this since it messes with the namespace
### it is, however, required below - will be referenced directly.
### so you will want to install it if you haven't already
if (TRUE == FALSE) {
  if (!required('arm')) {
    install.packages('arm')
  }
}
arm_rescale <- arm::rescale
library(dendextend)

# ---- get data ----

project <- "pici-1286"
biospec_sql <- 'SELECT * FROM [isb-cgc:tcga_201510_alpha.Biospecimen_data] where Study = "BLCA" LIMIT 1000'
bio_d <- query_exec(biospec_sql, project = project)

clin_sql <- 'SELECT * FROM [isb-cgc:tcga_201510_alpha.Clinical_data] where Study = "BLCA" LIMIT 1000'
clin_d <- query_exec(clin_sql, project = project)

somatic_summary_sql <- 
  '
  select
    ParticipantBarcode,
    Study, 
    Variant_Classification,
    count(*) as num_mutations
  FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]
  WHERE 
    Study = "BLCA" 
  GROUP BY 
    1,2,3
  '
somatic_d <- query_exec(somatic_summary_sql, project = project)

## save this for future research!!
somatic_detail_sql <- 
  '
  SELECT 
    ParticipantBarcode,
    Study, 
    Annotation_Transcript, 
    COSMIC_Total_Alterations_In_Gene, 
    Chromosome, 
    DNARepairGenes_Role, 
    End_Position, 
    Entrez_Gene_Id, 
    GC_Content, 
    GENCODE_Transcript_Name, 
    GENCODE_Transcript_Type, 
    GO_Biological_Process, 
    GO_Cellular_Component, 
    GO_Molecular_Function, 
    Variant_Classification 
  FROM [isb-cgc:tcga_201510_alpha.Somatic_Mutation_calls]
  WHERE 
    Study = "BLCA" 
    and Variant_Classification != "Silent"
  '

## ---- explore somatic mutations data ----


somatic_wide <- somatic_d %>%
  tidyr::spread(Variant_Classification, num_mutations, fill = 0)

#somatic_wide %>%
#  dplyr::select(-ParticipantBarcode, -Study) %>%
#  dplyr::select(-Missense_Mutation, -Nonsense_Mutation) %>%
#  pairs()

somatic_trans <- somatic_d %>%
  dplyr::select(-Study) %>%
  dplyr::mutate(log_mutations = log1p(num_mutations)) %>%
  dplyr::select(-num_mutations) %>%
  tidyr::spread(ParticipantBarcode, log_mutations, fill = log1p(0))

hcdata <- somatic_trans %>%
  dplyr::select(-Variant_Classification)
rownames(hcdata) <- somatic_trans$Variant_Classification
hcres <- hclust(dist(hcdata), "ave")
plot(hcres)


## nice-looking plot of 

hcdend <- as.dendrogram(hcres) %>%
  hang.dendrogram(hang_height=0.1) %>%
  sort(type = "nodes") %>%
  color_branches(k = 4)

par(mar = c(8,3,3,3))
plot(hcdend)
legend("topright", legend = paste('Cluster',cutree(hcdend, k = 4)[order.dendrogram(hcdend)] %>% unique()), fill = rainbow_hcl(4))

## add group number to the original data frame
hc_grps <- cutree(hcdend, k = 4)
hc_result <- data.frame(Variant_Cluster_ID = hc_grps
                        , Variant_Classification = names(hc_grps)
                        )

## add cluster ids back into original data
somatic_d2 <- 
  somatic_d %>%
  left_join(hc_result, by = 'Variant_Classification')

## print out group #s to confirm they are correct
somatic_d2 %>% 
  dplyr::distinct(Variant_Cluster_ID, Variant_Classification) %>% 
  dplyr::select(Variant_Classification, Variant_Cluster_ID) %>% 
  dplyr::arrange(Variant_Cluster_ID) %>%
  dplyr::distinct(Variant_Cluster_ID)


### ---- analysis plan ---- 

## let's assume these obervations are related to 4 separate processes - with varying levels of impact on the protein

## 1. underlying mutation rate (cluster 4)
## 2. regulatory / transcription modification
## 3. alternate splicing / fusion 
## 4. modification of protein function

## each of these occurs at a different rate. 

## one method of analysis would be to "sum" across each of these & analyze the sum of each type as a "biomarker". 
## an alternate method will treat each of these as a "type" & propose a latent term for each. 

## we will estimate each and compare the approaches. 

### ---- analysis of mutation rate using aggregate of measures ---- 

## first we aggregate each of these 4 types
somatic_d_agg <- 
  somatic_d2 %>%
  group_by(ParticipantBarcode, Variant_Cluster_ID) %>%
  summarise(mutations = sum(num_mutations)) %>%
  ungroup() 

## the resulting data look like this, with number of mutations per Participant * Cluster_ID : 
somatic_d_agg %>% head()

## next we reshape to 'wide' format, resulting in one obs per Participant
## named 'resc_' because we plan to rescale the values
somatic_d_agg_wide_resc <- 
  somatic_d_agg %>%
  dplyr::mutate(Variant_Cluster_ID = paste("resc_cluster_",Variant_Cluster_ID, sep = '')) %>%
  tidyr::spread(Variant_Cluster_ID, mutations, fill = 0)

## we now have the count of each "cluster" of mutations per participant
somatic_d_agg_wide_resc %>% head()

## finally, we rescale each of the four mutation count types 
somatic_d_agg_wide_resc <- 
  somatic_d_agg_wide_resc %>%
  mutate_each(funs = funs(arm_rescale), starts_with('resc_cluster_'))

## let's keep the original data as well 
somatic_d_agg_wide <- 
  somatic_d_agg %>%
  dplyr::mutate(Variant_Cluster_ID = paste("cluster_",Variant_Cluster_ID, sep = '')) %>%
  tidyr::spread(Variant_Cluster_ID, mutations, fill = 0) %>%
  dplyr::inner_join(somatic_d_agg_wide_resc, by = 'ParticipantBarcode')



## notice that these are still highly correlated measures : 
pairs(somatic_d_agg_wide %>% dplyr::select(-ParticipantBarcode))


## ---- merge with biospec data ----

## note - exploration of biospec data done elsewhere
## here we will simply prep the data 
tumor_samples <- 
  bio_d %>%
  dplyr::filter(SampleType == 'Primary solid Tumor' & is_ffpe == 'NO')

stopifnot(!duplicated(tumor_samples$patid))

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
                   ) %>%
  dplyr::left_join(somatic_d_agg_wide, by = 'ParticipantBarcode')


## ---- look at cluster mutation counts <> survival ---- 

## let's see if this is a useful indicator in our sample overall
xph1 <- coxph(Surv(surv_time, mortality) ~
               cluster_1,
             data=d
)
summary(xph)

## let's compare the coefficient estimates (correlation with survival) between these mutation event
## types to see if the impact of each mutation is similar across these types 

est_cox <- function(coefname, data = d) {
  coxph(as.formula(paste("Surv(surv_time, mortality) ~ ",coefname,sep='')), data = data)
}
## estimate model for each coefficient 
model_results <- c(1:4) %>%
  purrr::map_chr(~paste('cluster_',.,sep='')) %>%
  purrr::map(est_cox) 
## extract coefficient estimates
model_results %>%
  purrr::map_dbl(coef)

## same using rescaled variables
## estimate model for each coefficient 
model_results_rescaled <- c(1:4) %>%
  purrr::map_chr(~paste('resc_cluster_',.,sep='')) %>%
  purrr::map(est_cox) 
## extract coefficient estimates
model_results_rescaled %>%
  purrr::map_dbl(coef)






