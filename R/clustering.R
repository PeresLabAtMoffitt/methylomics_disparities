# Import library
library(tidyverse)
library(cluster)
library(factoextra)


######################################################### I ### Load data
remp_res_Alu <- read_rds(paste0(here::here(), "/intermediary data/remp_results_Alu.rds"))
remp_res_ERV <- read_rds(paste0(here::here(), "/intermediary data/remp_results_ERV.rds"))
remp_res_L1 <- read_rds(paste0(here::here(), "/intermediary data/remp_results_L1.rds"))


######################################################### II ### Prep data for clustering
# Use the non trimmed data to create cluster because need to not have any NAs
# Aggregation step is then done manually to eliminate RE repeat

# L1
# Get Beta
L1_beta_results <- rempB(remp_res_L1)
L1_beta_results <- L1_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_L1@rowRanges@elementMetadata$RE.Index, .) %>% 
  select(-"X203717920012_R01C01")

L1_df <- L1_beta_results %>% 
  group_by(Index) %>% 
  mutate(n = n(), .before = 2) %>% 
  # Summarize RE when at least 2 predicted CpG sites
  # As rempAggregate() use the mean
  filter(n > 1) %>% select(-n) %>% 
  summarise(across(everything(), mean))

L1_df <- L1_df %>%
  column_to_rownames("Index")
a <- t(L1_df)
# Remove rows with missing values # Not needed now as we use the non trimmed data
# a1 <- na.omit(L1_df)
scaled_L1_df <- scale(a)


# Alu
# Get Beta
Alu_beta_results <- rempB(remp_res_Alu)
Alu_beta_results <- Alu_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_Alu@rowRanges@elementMetadata$RE.Index, .) %>% 
  select(-"X203717920012_R01C01")

Alu_df <- Alu_beta_results %>% 
  group_by(Index) %>% 
  mutate(n = n(), .before = 2) %>% 
  # Summarize RE when at least 2 predicted CpG sites
  # As rempAggregate() use the mean
  filter(n > 1) %>% select(-n) %>% 
  summarise(across(everything(), mean))

Alu_df <- Alu_df %>%
  column_to_rownames("Index")
a <- t(Alu_df)
scaled_Alu_df <- scale(a)


# ERV
# Get Beta
ERV_beta_results <- rempB(remp_res_ERV)
ERV_beta_results <- ERV_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_ERV@rowRanges@elementMetadata$RE.Index, .) %>% 
  select(-"X203717920012_R01C01")

ERV_df <- ERV_beta_results %>% 
  group_by(Index) %>% 
  mutate(n = n(), .before = 2) %>% 
  # Summarize RE when at least 2 predicted CpG sites
  # As rempAggregate() use the mean
  filter(n > 1) %>% select(-n) %>% 
  summarise(across(everything(), mean))

ERV_df <- ERV_df %>%
  column_to_rownames("Index")
a <- t(ERV_df)
scaled_ERV_df <- scale(a)

# Overall
full_beta <- bind_rows(Alu_df, L1_df, ERV_df)
a <- t(full_beta)
scaled_full_beta <- scale(a)
write_rds(scaled_full_beta, "scaled_full_beta.rds")


######################################################### III ### Clustering
# Create plot of number of clusters vs total within sum of squares
fviz_nbclust(scaled_L1_df, kmeans, method = "wss")
jpeg("L1 elbow plot.jpeg")
fviz_nbclust(scaled_L1_df, kmeans, method = "wss")
dev.off()
fviz_nbclust(scaled_Alu_df, kmeans, method = "wss")
jpeg("Alu elbow plot.jpeg")
fviz_nbclust(scaled_Alu_df, kmeans, method = "wss")
dev.off()
fviz_nbclust(scaled_ERV_df, kmeans, method = "wss")
jpeg("ERV elbow plot.jpeg")
fviz_nbclust(scaled_ERV_df, kmeans, method = "wss")
dev.off()
fviz_nbclust(scaled_full_beta, kmeans, method = "wss")
jpeg("full RE beta data elbow plot.jpeg")
fviz_nbclust(scaled_full_beta, kmeans, method = "wss")
dev.off()

# Perform k-means clustering with k clusters
set.seed(123)
L1_km <- kmeans(scaled_L1_df, centers = 2, nstart = 25)
Alu_km <- kmeans(scaled_Alu_df, centers = 2, nstart = 25)
ERV_km <- kmeans(scaled_ERV_df, centers = 2, nstart = 25)
RE_km_2 <- kmeans(scaled_full_beta, centers = 2, nstart = 25)
RE_km_4 <- kmeans(scaled_full_beta, centers = 4, nstart = 25)
RE_km_5 <- kmeans(scaled_full_beta, centers = 5, nstart = 25)

# Bind cluster
Alu_cluster <- as_tibble(Alu_km$cluster) %>% `colnames<-`(c("cluster_Alu"))
L1_cluster <- as_tibble(L1_km$cluster) %>% `colnames<-`(c("cluster_L1"))
ERV_cluster <- as_tibble(ERV_km$cluster) %>% `colnames<-`(c("cluster_ERV"))
RE_cluster_2 <- as_tibble(RE_km_2$cluster) %>% `colnames<-`(c("cluster_RE_2"))
RE_cluster_4 <- as_tibble(RE_km_4$cluster) %>% `colnames<-`(c("cluster_RE_4"))
RE_cluster_5 <- as_tibble(RE_km_5$cluster) %>% `colnames<-`(c("cluster_RE_5"))

final_RE <- bind_cols(patient_id = colnames(Alu_beta_results)[2:ncol(Alu_beta_results)],
                      Alu_cluster, L1_cluster, ERV_cluster, 
                      RE_cluster_2, RE_cluster_4, RE_cluster_5)

write_rds(final_RE, "final_RE.rds")

######################################################### IV ### Extract non correlated RE
library(caret)
Alu_vector <- findCorrelation(
  Alu_df,
  cutoff = 5,
  verbose = TRUE,
  names = TRUE,
  exact = FALSE
)
Alu_vector
write_rds(Alu_vector, paste0(here::here(), "/intermediary data/Alu_vector.rds"))
L1_vector <- findCorrelation(
  L1_df,
  cutoff = 5,
  verbose = TRUE,
  names = TRUE,
  exact = FALSE
)
L1_vector
write_rds(L1_vector, paste0(here::here(), "/intermediary data/L1_vector.rds"))
ERV_vector <- findCorrelation(
  ERV_df,
  cutoff = 5,
  verbose = TRUE,
  names = TRUE,
  exact = FALSE
)
ERV_vector
write_rds(ERV_vector, paste0(here::here(), "/intermediary data/ERV_vector.rds"))
RE_vector <- findCorrelation(
  scaled_full_beta,
  cutoff = 5,
  verbose = TRUE,
  names = TRUE,
  exact = FALSE
)
write_rds(RE_vector, paste0(here::here(), "/intermediary data/RE_vector.rds"))




# END Clustering

