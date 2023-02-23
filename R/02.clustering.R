# Import library
library(tidyverse)
library(cluster)
library(factoextra)


######################################################### I ### Load data
remp_res_Alu <- read_rds(paste0(here::here(), "/remp_res_Alu.rds"))
remp_res_ERV <- read_rds(paste0(here::here(), "/remp_res_ERV.rds"))
remp_res_L1 <- read_rds(paste0(here::here(), "/remp_res_L1.rds"))


######################################################### II ### Prep data for clustering
# Use the non trimmed data to create cluster because need to not have any NAs
# remp_res_L1 <- rempTrim(remp_res_L1,threshold=3,missingRate=0.2)
# # Also proceed with the aggregation step to eliminate RE repeat
# remp_res_L1 <- rempAggregate(remp_res_L1, NCpG = 2, ncore = 4)

# L1
# Get Beta
L1_beta_results <- rempB(remp_res_L1)
L1_beta_results <- L1_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_L1@rowRanges@elementMetadata$RE.Index, .)

L1_df <- L1_beta_results %>% 
  group_by(Index) %>% 
  mutate(n = n(), .before = 2) %>% 
  # Summarize RE when at least 2 predicted CpG sites ss rempAggregate() does
  filter(n > 1) %>% select(-n) %>% 
  # As rempAggregate() use the mean
  summarise(across(everything(), mean))

L1_df <- L1_df %>%
  column_to_rownames("Index")
# a <- t(map_df1) %>% as_tibble() %>% mutate(rw = row.names(t(map_df1)), .before = 1)
a <- t(L1_df)
# Remove rows with missing values
# a1 <- na.omit(L1_df)

L1_df <- scale(a) #%>% as_tibble(rownames = "patient_id") %>% 
  # column_to_rownames("patient_id")

# Alu
# Get Beta
Alu_beta_results <- rempB(remp_res_Alu)
Alu_beta_results <- Alu_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_Alu@rowRanges@elementMetadata$RE.Index, .)

Alu_df <- Alu_beta_results %>% 
  group_by(Index) %>% 
  mutate(n = n(), .before = 2) %>% 
  # Summarize RE when at least 2 predicted CpG sites ss rempAggregate() does
  filter(n > 1) %>% select(-n) %>% 
  # As rempAggregate() use the mean
  summarise(across(everything(), mean))

Alu_df <- Alu_df %>%
  column_to_rownames("Index")
a <- t(Alu_df)
Alu_df <- scale(a)

# ERV
# Get Beta
ERV_beta_results <- rempB(remp_res_ERV)
ERV_beta_results <- ERV_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_ERV@rowRanges@elementMetadata$RE.Index, .)

ERV_df <- ERV_beta_results %>% 
  group_by(Index) %>% 
  mutate(n = n(), .before = 2) %>% 
  # Summarize RE when at least 2 predicted CpG sites ss rempAggregate() does
  filter(n > 1) %>% select(-n) %>% 
  # As rempAggregate() use the mean
  summarise(across(everything(), mean))

ERV_df <- ERV_df %>%
  column_to_rownames("Index")
a <- t(ERV_df)
ERV_df <- scale(a)


######################################################### III ### Clustering
# Create plot of number of clusters vs total within sum of squares
fviz_nbclust(L1_df, kmeans, method = "wss")
jpeg("L1 elbow plot.jpeg")
fviz_nbclust(L1_df, kmeans, method = "wss")
dev.off()
fviz_nbclust(Alu_df, kmeans, method = "wss")
jpeg("Alu elbow plot.jpeg")
fviz_nbclust(Alu_df, kmeans, method = "wss")
dev.off()
fviz_nbclust(ERV_df, kmeans, method = "wss")
jpeg("ERV elbow plot.jpeg")
fviz_nbclust(ERV_df, kmeans, method = "wss")
dev.off()

# Perform k-means clustering with k clusters
set.seed(123)
L1_km <- kmeans(L1_df, centers = 2, nstart = 25)
Alu_km <- kmeans(Alu_df, centers = 2, nstart = 25)
ERV_km <- kmeans(ERV_df, centers = 2, nstart = 25)

# Bind cluster
Alu_cluster <- as_tibble(Alu_km$cluster) %>% `colnames<-`(c("Alu_cluster"))
L1_cluster <- as_tibble(L1_km$cluster) %>% `colnames<-`(c("L1_cluster"))
ERV_cluster <- as_tibble(ERV_km$cluster) %>% `colnames<-`(c("ERV_cluster"))

final_RE <- bind_cols(patient_id = colnames(Alu_beta_results)[2:ncol(Alu_beta_results)],
                      Alu_cluster, L1_cluster, ERV_cluster)

write_rds(final_RE, "final_RE.rds")

# END

