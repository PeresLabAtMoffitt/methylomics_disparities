# Import library
library(tidyverse)
library(mice)
library(cluster)
library(factoextra)


######################################################### I ### Load data
remp_res_Alu <- read_rds(paste0(here::here(), "/remp_res_Alu_annotation.rds"))
remp_res_ERV <- read_rds(paste0(here::here(), "/remp_res_ERV_annotation.rds"))
remp_res_L1 <- read_rds(paste0(here::here(), "/remp_res_L1_annotation.rds"))


######################################################### II ### Prep data for clustering
# Get Beta
L1_beta_results <- rempB(remp_res_L1)
L1_beta_results <- L1_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_L1@rowRanges@elementMetadata$RE.Index, .)

L1_df <- L1_beta_results %>%
  column_to_rownames("Index")
# a <- t(map_df1) %>% as_tibble() %>% mutate(rw = row.names(t(map_df1)), .before = 1)
a <- t(L1_df)
L1_df <- scale(a) %>% as_tibble() %>% mutate(patient_id = row.names(t(L1_df)), .before = 1)


Alu_beta_results <- rempB(remp_res_Alu)
Alu_beta_results <- Alu_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_Alu@rowRanges@elementMetadata$RE.Index, .)

Alu_df <- Alu_beta_results %>%
  column_to_rownames("Index")
# a <- t(map_df1) %>% as_tibble() %>% mutate(rw = row.names(t(map_df1)), .before = 1)
a <- t(Alu_df)
Alu_df <- scale(a) %>% as_tibble() %>% mutate(patient_id = row.names(t(Alu_df)), .before = 1)


ERV_beta_results <- rempB(remp_res_ERV)
ERV_beta_results <- ERV_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_ERV@rowRanges@elementMetadata$RE.Index, .)

ERV_df <- ERV_beta_results %>%
  column_to_rownames("Index")
# a <- t(map_df1) %>% as_tibble() %>% mutate(rw = row.names(t(map_df1)), .before = 1)
a <- t(ERV_df)
ERV_df <- scale(a) %>% as_tibble() %>% mutate(patient_id = row.names(t(ERV_df)), .before = 1)


RE_data <- full_join(Alu_df, L1_df, by = "patient_id") %>% 
  full_join(., ERV_df, by = "patient_id")

write_rds(RE_data, "RE_data.rds")
RE_data <- read_rds(paste0(here::here(), "/RE_data.rds"))

RE_data_imputed <- mice(RE_data, m=5, maxit=50, 
             seed=123, printFlag=F)
write_rds(RE_data_imputed, "RE_data_imputed.rds")
RE_data_imputed <- read_rds(paste0(here::here(), "/RE_data_imputed.rds"))
# check1 <- with(data = RE_data_imputed,
#                exp = coxph(
#                  Surv(pfs_time, pfs_event) ~ 
#                    raceeth1 + age_at_diagnosis + tnm_cs_mixed_group_stage + treatment_type + bmi
#                ))
# multiple_imputations_number <- round(max(pool(check1)$pooled$fmi),2)*100
# RE_data_imputed <- mice(RE_data, m=multiple_imputations_number, maxit=50, 
#                         seed=123, printFlag=F)



df <- USArrests

#remove rows with missing values
RE_data1 <- na.omit(RE_data)

#scale each variable to have a mean of 0 and sd of 1
df <- scale(df)

#view first six rows of dataset
head(df)




######################################################### III ### Clustering


#create plot of number of clusters vs total within sum of squares
fviz_nbclust(RE_data_imputed, kmeans, method = "wss", verbose = TRUE) # "wss" (for total within sum of square) 
RE_data_imputed %>% select(start)
fviz_nbclust(df, kmeans, method = "wss")
fviz_nbclust(df, kmeans, method = "wss")
fviz_nbclust(df, kmeans, method = "wss")
jpeg("elbow plot.jpeg")
fviz_nbclust(df, kmeans, method = "wss")
dev.off()

#make this example reproducible
set.seed(123)

#perform k-means clustering with k = 4 clusters
km <- kmeans(df, centers = 4, nstart = 25)

#view results
km

#add cluster assigment to original data
final_data <- cbind(data, cluster = km$cluster)

#view final data
head(final_data)



