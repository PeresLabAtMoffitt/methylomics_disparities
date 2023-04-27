########################################################## 
################## Consensus Clustering ##################
##########################################################

# Import library
library(tidyverse)
library(ConsensusClusterPlus)


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
  # Remove this sample as lots of NAs
  select(-"X203717920012_R01C01")

# Summarize to remove duplicate measurement
L1_df <- L1_beta_results %>%
  group_by(Index) %>%
  mutate(n = n(), .before = 2) %>%
  # Summarize RE when at least 2 predicted CpG sites
  # As rempAggregate() use the mean
  filter(n > 1) %>% select(-n) %>%
  summarise(across(everything(), mean))

# The input data format is a matrix where columns are samples (items), rows are features
L1_df_ <- L1_df %>%
  column_to_rownames("Index") %>% 
  as.matrix()

# Select the most informative featues
# Reduce the data set to the top 5,000 most variable genes, measured by median absolute deviation (MAD).
mads <- apply(L1_df_, 1, mad) 
L1_df <- L1_df_[rev(order(mads))[1:5000],]

# Transform the data with agglomerative hierarchical clustering algorithm using Pearson correlation distance
ready_data <- sweep(L1_df, 1, apply(L1_df, 1, median,na.rm=T))

# Run ConsensusClusterPlus
title <- paste0(getwd(), "/consensus plots")
results <- ConsensusClusterPlus(
  ready_data, # data to use
  maxK = 20, # a maximum evalulated k of 20 so that cluster counts of 2,3,4,...,20 are evaluated
  reps = 1000, # 1000 resamplings
  pItem = 0.8, # 80% item resampling
  pFeature = 1, # 80% gene resampling
  title = title, # output
  clusterAlg = "hc", # hierarchical clustering
  distance = "pearson",
  seed = 123, plot = "pngBMP"
  )

# Explore results
# Results for 3 k cluster
results[[3]][["consensusMatrix"]][1:5,1:5]

#consensusClass-thesampleclassifications 
results[[3]][["consensusClass"]][1:5]


# Calculate cluster consensus and item-consensus results
icl <- calcICL(results, title = title, plot = "pngBMP"
               )

icl[["clusterConsensus"]]

icl[["itemConsensus"]][1:10,]











