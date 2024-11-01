# Import library
library(tidyverse)
library(REMP)
# library(gtsummary)
library(patchwork)
# theme_gtsummary_compact()
theme_set(theme_classic())

# Figure S1 patients methylation pattern ----
# Load data
# LTR
remp_res_ERV <- read_rds(paste0(here::here(), "/intermediary data/remp_res_ERV_annotation_07262023.rds"))
remplot(remp_res_ERV, main = "ERV methylation", col = "blue")
annot_ERV_beta_results <- rempB(remp_res_ERV)
map_LTR <- annot_ERV_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_ERV@rowRanges@elementMetadata$RE.Index, .)

LTR <- map_LTR %>% 
  pivot_longer(cols = - Index) %>% 
  ggplot(aes(x= value, group= name))+
  geom_density(color = "black", linewidth = 0.1)+
  labs(x = "Methylation value (beta)")+
  xlim(0,1)+
  # scale_y_log10()+
  scale_y_continuous(transform = "log2"#,
                     # breaks = c(0.8, 1, 3, 10),
                     # labels = c(0.8, 1, 3, 10),
                     # limits = c(0.8,10),
                     # expand = c(-0.0001, 0.1)
  )+
  
  theme(legend.position = "none")

# L1
remp_res_L1 <- read_rds(paste0(here::here(), "/intermediary data/remp_res_L1_annotation_07262023.rds"))
remplot(remp_res_L1, main = "L1 methylation", col = "blue")
annot_L1_beta_results <- rempB(remp_res_L1)
map_L1 <- annot_L1_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_L1@rowRanges@elementMetadata$RE.Index, .)
library(ggbreak)
L1 <- map_L1 %>% 
  pivot_longer(cols = - Index) %>% 
  ggplot(aes(x= value, group= name))+
  geom_density(color = "black", linewidth = 0.1)+
  labs(x = "Methylation value (beta)")+
  xlim(0,1)+
  # scale_y_continuous(breaks = seq(0,40,10), limits = c(0,40)) +
  scale_y_break(c(6,15), #scales = "free",# space = 0.4,
                # ticklabels = c(seq(0,40,5)),
                # expand = c(0, 0)
                ) +
  theme(legend.position = "none")

# Alu
remp_res_Alu <- read_rds(paste0(here::here(), "/intermediary data/remp_res_Alu_annotation_07262023.rds"))
remplot(remp_res_Alu, main = "Alu methylation", col = "blue")
annot_Alu_beta_results <- rempB(remp_res_Alu)
map_ALU <- annot_Alu_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_Alu@rowRanges@elementMetadata$RE.Index, .)

ALU <- map_ALU %>% 
  pivot_longer(cols = - Index) %>% 
  ggplot(aes(x= value, group= name))+
  geom_density(color = "black", linewidth = 0.1)+
  labs(x = "Methylation value (beta)")+
  xlim(0,1)+
  theme(legend.position = "none")


L1 / ALU / LTR +
  plot_annotation(tag_levels = "A")+
  plot_layout(axes = "collect"#,
              # axis_titles = "collect"
              )

ggsave("Figure S1 patients methylation pattern.pdf",
       width = 4,
       height = 10, 
       dpi = 600)


# Batch scatter plot
library(tidyverse)
library("SQMtools")
theme_set(theme_classic())
met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
# cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))
# plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster.rds"))
# 
# alu_dat <- plot_cluster$Alu_dat # patient as column
# l1_dat <- plot_cluster$L1_dat
# erv_dat <- plot_cluster$ERV_dat

# Select most variable CPGs
dat <- bind_rows(map_ALU, map_L1, map_LTR) %>% 
  column_to_rownames("Index") %>% 
  t()

set.seed(1234)
dat1 <- mostVariable(dat, N=1000, bycol = TRUE) %>% 
  as_tibble(rownames = "patient_id") %>% 
  full_join(met_data %>% 
              select(patient_id, Batch) %>% 
              mutate(patient_id = paste0("X", patient_id)),
            .,
            by = "patient_id")

dat2 <- dat1 %>% select(-patient_id) %>% 
  group_by(Batch) %>% 
  summarize_all(.funs = ~ mean(.))

dat2 %>% 
  pivot_longer(cols = -c(Batch)) %>% 
  ggplot(aes(x= value, y= value, colour = Batch))+
  geom_point()



# PCA batch----
library(tidyverse)
library(RColorBrewer)
load(paste0(here::here(), "/cleaned_07082022.rda"))

phenoclean <- phenoclean %>% 
  mutate(Batch = case_when(
    Batch == "Peres"      ~ "Batch 1",
    Batch == "Peres2"     ~ "Batch 2",
    Batch == "Doherty"    ~ "Batch 3",
  ))
batch_col <- phenoclean$Batch # %>% select(Complete.Barcode, Batch)

betas_clean6[1:10, 1:5]
beta <- betas_clean6
# beta

pal <- brewer.pal(3,"Dark2")

# plotMDS(beta, top=2)
library(limma)
plotMDS(beta, top=1000, gene.selection="common", 
        pch = c(0, 15, 5),
        col=pal[factor(batch_col)])
# ggsave("Figure S1 batch-limna.pdf",
#        width = 4,
#        height = 10, 
#        dpi = 600)
library(minfi)
mdsPlot(beta, numPositions = 1000, 
        main = "", 
        sampGroups = batch_col, 
        pch = 1,
        pal = c(pal), 
        # xlim = c(-8, 8),
        # ylim = c(-7, 9),
        legendPos = "topright", 
        legendNCol = 1)
title(xlab = "Dimension 1", ylab = "Dimension 2")
# ggsave("Figure S1 batch-minfi with legend.pdf", # Doesn't work, need to extract h = 6, w = 7
#        width = 4,
#        height = 10,
#        dpi = 600)




# # Heatmap ComplexHeatmap NOT RUN ----
# library(tidyverse)
# library(ComplexHeatmap)
# theme_set(theme_classic())
# met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
# cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))
# plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster.rds"))
# 
# alu_dat <- plot_cluster$Alu_dat # patient as column
# 
# l1_dat <- plot_cluster$L1_dat
# 
# erv_dat <- plot_cluster$ERV_dat
# 
# heatmap_df <- t(alu_dat) %>% as_tibble(rownames = "patient_id")
# 
# annCol <- met_data %>% select(patient_id, 
#                               BRCA1_carrier, BRCA2_carrier
# ) %>% 
#   mutate(patient_id = paste0("X", patient_id)) %>% 
#   mutate(across(everything(), ~ as.character(.))) %>% 
#   mutate(across(everything(), ~ replace_na(., "Unknown")))
# coca_cluster_results <- cluster_res_list1$COCA$clust.res %>%
#   remove_rownames()
# 
# coca_cluster_results <- coca_cluster_results %>% 
#   full_join(., annCol, 
#             by = c("samID" = "patient_id")) %>% 
#   filter(BRCA1_carrier != "Unknown" |
#            BRCA2_carrier != "Unknown")
# 
# heatmap_df1 <- heatmap_df %>%
#   left_join(., coca_cluster_results, 
#             by = c("patient_id" = "samID")) %>% 
#   arrange(clust)
# heatmap_beta <- heatmap_df1 %>% 
#   unite(patient_id, c(patient_id, clust, BRCA1_carrier, BRCA2_carrier), sep = "_", remove = TRUE) %>%
#   column_to_rownames(var = "patient_id") %>% 
#   t()
# 
# Heatmap(heatmap_beta)
# 
# column_ho = HeatmapAnnotation(clust = c(heatmap_df1$clust),
#                               BRCA1_carrier = c(heatmap_df1$BRCA1_carrier), 
#                               col = list(clust = c("1" = "#932667FF", "2" = "#FDE725FF"),
#                                          BRCA1_carrier = c("Yes" = "red", "No"= "blue", "Unknown" = "black")
#                               ),
#                               na_col = "black")
# Heatmap(heatmap_beta, name = " ", 
#         # cluster_rows = FALSE, 
#         cluster_columns = FALSE,
#         top_annotation = column_ho)
# 





# Heatmap MOVICS----
# Each cluster separately
library(tidyverse)
# library(MOVICS)
library(ComplexHeatmap)
source("R/getMoHeatmap_updated_MOVICSfunc.R") # Need to chane the name of each cluster in the function
theme_set(theme_classic())
# met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))
plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster.rds"))
# plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster5000.rds"))

### NEMO
getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F), 
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$NEMO$clust.res, # cluster results
             clust.dend    = cluster_res_list1$NEMO$clust.dend, # show dendrogram for samples
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "heatmaps/Comprehensive heatmap of NEMO_10022024")



### COCA
getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$COCA$clust.res, 
             # clust.dend    = cluster_res_list1$COCA$clust.dend, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             fig.name      = "heatmaps/Comprehensive heatmap of COCA_10022024")



### ConsensusClustering
getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$ConsensusClustering$clust.res, 
             # clust.dend    = cluster_res_list1$ConsensusClustering$clust.dend, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             fig.name      = "heatmaps/Comprehensive heatmap of ConsensusClustering_10022024")



### IntNMF
getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$IntNMF$clust.res, 
             # clust.dend    = cluster_res_list1$IntNMF$clust.dend, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             fig.name      = "heatmaps/Comprehensive heatmap of IntNMF_10022024")



### MoCluster
getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$MoCluster$clust.res, 
             # clust.dend    = cluster_res_list1$MoCluster$clust.dend, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             fig.name      = "heatmaps/Comprehensive heatmap of MoCluster_10022024")



### iClusterBayes
getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$iClusterBayes$clust.res, 
             # clust.dend    = cluster_res_list1$iClusterBayes$clust.dend, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             fig.name      = "heatmaps/Comprehensive heatmap of iClusterBayes_10022024")


### BCC

cluster_res_list1$BCC$clust.res <- cluster_res_list1$BCC$clust.res %>%
  mutate(clust = case_when(
    clust == 1 ~ 2,
    clust == 2 ~ 1
  ))

getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","ERV"),
             is.binary     = c(F,F,F), 
             legend.name   = c("L1","Alu","ERV"),
             clust.res     = cluster_res_list1$BCC$clust.res, # cluster results
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "heatmaps/Heatmap of BCC_10022024")


# Heatmap BRCA 1 and 2----
library(tidyverse)
# library(MOVICS)
library(ComplexHeatmap)
source("R/getMoHeatmap_updated_MOVICSfunc.R") # Need to change the name of each cluster in the function
theme_set(theme_classic())
met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))
plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster.rds"))
# plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster5000.rds"))

annCol <- met_data %>% select(patient_id, 
                              BRCA1_carrier, BRCA2_carrier
) %>% 
  mutate(patient_id = paste0("X", patient_id)) %>% 
  mutate(across(everything(), ~ as.character(.))) %>% 
  mutate(across(everything(), ~ replace_na(., "Unknown")))
coca_cluster_results <- cluster_res_list1$COCA$clust.res %>%
  remove_rownames()

coca_cluster_results <- coca_cluster_results %>% 
  full_join(., annCol, 
            by = c("samID" = "patient_id")) %>% 
  filter(BRCA1_carrier != "Unknown" |
           BRCA2_carrier != "Unknown")

cluster_res_list1$COCA$clust.res <- cluster_res_list1$COCA$clust.res %>% 
  filter(str_detect(samID, paste0(coca_cluster_results$samID, collapse = "|")
                    )) %>% 
  mutate(clust = case_when(
    clust == 1 ~ "Active",
    clust == 2 ~ "Repressed"
  ))

annCol <- annCol %>% 
  column_to_rownames("patient_id") %>% 
  filter(BRCA1_carrier != "Unknown" |
           BRCA2_carrier != "Unknown") %>% 
  mutate(BRCA1_carrier = factor(BRCA1_carrier, levels = c("No", "Yes", "Unknown"))) %>% 
  mutate(BRCA2_carrier = factor(BRCA2_carrier, levels = c("No", "Yes", "Unknown")))

annColors <- list(BRCA1_carrier  = c("Yes" = "purple",
                                               "No"   = "lightgrey",
                                               "Unknown"   = "black"),
                  BRCA2_carrier  = c("Yes" = "purple",
                                            "No"   = "lightgrey",
                                            "Unknown"   = "black")
)

getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","LTR"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","LTR"),
             clust.res     = cluster_res_list1$COCA$clust.res,
             clust.dend    = NULL, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             fig.name      = "heatmaps/Heatmap of COCA_with_limited_BRCA_10022024"
             )

# # BRCA 1 and 2 - overall/germline/tumor
# annCol <- met_data %>% select(patient_id, 
#                               BRCA1_carrier, BRCA2_carrier,
#                               germline_mutation_BRCA1, tumor_mutation_BRCA1, 
#                               germline_mutation_BRCA2, tumor_mutation_BRCA2
# ) %>% 
#   mutate(patient_id = paste0("X", patient_id)) %>% 
#   mutate(across(everything(), ~ as.character(.))) %>% 
#   mutate(across(everything(), ~ replace_na(., "Unknown")))
# coca_cluster_results <- cluster_res_list1$COCA$clust.res %>%
#   remove_rownames()
# 
# coca_cluster_results <- coca_cluster_results %>% 
#   full_join(., annCol, 
#             by = c("samID" = "patient_id")) %>% 
#   filter(BRCA1_carrier != "Unknown" |
#            BRCA2_carrier != "Unknown")
# 
# cluster_res_list1$COCA$clust.res <- cluster_res_list1$COCA$clust.res %>% 
#   filter(str_detect(samID, paste0(coca_cluster_results$samID, collapse = "|")
#   ))
# 
# annCol <- annCol %>% 
#   column_to_rownames("patient_id") %>% 
#   filter(BRCA1_carrier != "Unknown" |
#            BRCA2_carrier != "Unknown")
# 
# annColors <- list(BRCA1_carrier  = c("Yes" = "red",
#                                      "No"   = "blue",
#                                      "Unknown"   = "lightgrey"),
#                   BRCA2_carrier  = c("Yes" = "red",
#                                      "No"   = "blue",
#                                      "Unknown"   = "lightgrey"),
#                   germline_mutation_BRCA1  = c("Yes" = "red",
#                                                "No"   = "blue",
#                                                "Unknown"   = "lightgrey"),
#                   tumor_mutation_BRCA1  = c("Yes" = "red",
#                                             "No"   = "blue",
#                                             "Unknown"   = "lightgrey"),
#                   germline_mutation_BRCA2  = c("Yes" = "red",
#                                                "No"   = "blue",
#                                                "Unknown"   = "lightgrey"),
#                   tumor_mutation_BRCA2  = c("Yes" = "red",
#                                             "No"   = "blue",
#                                             "Unknown"   = "lightgrey")
# )
# 
# getMoHeatmap(data          = plot_cluster,
#              row.title     = c("L1","Alu","LTR"),
#              is.binary     = c(F,F,F),
#              legend.name   = c("L1","Alu","LTR"),
#              clust.res     = cluster_res_list1$COCA$clust.res,
#              clust.dend    = NULL, 
#              color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
#                                   c("#0000FF", "white"  , "#FF0000"),
#                                   c("#0000FF", "white"  , "#FF0000")),
#              width         = 10, 
#              height        = 5, 
#              annCol        = annCol, # annotation for samples
#              annColors     = annColors, # annotation color
#              fig.name      = "heatmaps/Heatmap of COCA_with_limited_BRCA_6bars_05102024"
# )
# 
# # BRCA 1 and 2 - overall/germline/tumor
# annCol <- met_data %>% select(patient_id, 
#                               germline_mutation_BRCA1, tumor_mutation_BRCA1, 
#                               germline_mutation_BRCA2, tumor_mutation_BRCA2
# ) %>% 
#   mutate(patient_id = paste0("X", patient_id)) %>% 
#   mutate(across(everything(), ~ as.character(.))) %>% 
#   mutate(across(everything(), ~ replace_na(., "Unknown")))
# coca_cluster_results <- cluster_res_list1$COCA$clust.res %>%
#   remove_rownames()
# 
# coca_cluster_results <- coca_cluster_results %>% 
#   full_join(., annCol, 
#             by = c("samID" = "patient_id")) %>% 
#   filter(BRCA1_carrier != "Unknown" |
#            BRCA2_carrier != "Unknown")
# 
# cluster_res_list1$COCA$clust.res <- cluster_res_list1$COCA$clust.res %>% 
#   filter(str_detect(samID, paste0(coca_cluster_results$samID, collapse = "|")
#   ))
# 
# annCol <- annCol %>% 
#   column_to_rownames("patient_id") %>% 
#   filter(BRCA1_carrier != "Unknown" |
#            BRCA2_carrier != "Unknown")
# 
# annColors <- list(germline_mutation_BRCA1  = c("Yes" = "red",
#                                                "No"   = "blue",
#                                                "Unknown"   = "lightgrey"),
#                   tumor_mutation_BRCA1  = c("Yes" = "red",
#                                             "No"   = "blue",
#                                             "Unknown"   = "lightgrey"),
#                   germline_mutation_BRCA2  = c("Yes" = "red",
#                                                "No"   = "blue",
#                                                "Unknown"   = "lightgrey"),
#                   tumor_mutation_BRCA2  = c("Yes" = "red",
#                                             "No"   = "blue",
#                                             "Unknown"   = "lightgrey")
# )
# 
# getMoHeatmap(data          = plot_cluster,
#              row.title     = c("L1","Alu","LTR"),
#              is.binary     = c(F,F,F),
#              legend.name   = c("L1","Alu","LTR"),
#              clust.res     = cluster_res_list1$COCA$clust.res,
#              clust.dend    = NULL, 
#              color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
#                                   c("#0000FF", "white"  , "#FF0000"),
#                                   c("#0000FF", "white"  , "#FF0000")),
#              width         = 10, 
#              height        = 5, 
#              annCol        = annCol, # annotation for samples
#              annColors     = annColors, # annotation color
#              fig.name      = "heatmaps/Heatmap of COCA_with_limited_BRCA_4bars_05102024"
# )

# 1 heatmap for all cluster----
library(tidyverse)
library(ComplexHeatmap)
source("R/allclusters_getMoHeatmap_updated_MOVICSfunc.R")
# library(MOVICS)
theme_set(theme_classic())
met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))
plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster.rds"))
# plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster5000.rds"))

# BRCA 1 and 2
annCol <- met_data %>% select(patient_id, 
                              BRCA1_carrier, BRCA2_carrier
) %>% 
  mutate(patient_id = paste0("X", patient_id)) %>% 
  mutate(across(everything(), ~ as.character(.))) %>% 
  mutate(across(everything(), ~ replace_na(., "Unknown")))


coca_cluster_results <- cluster_res_list1$COCA$clust.res %>%
  rename(COCA = clust) %>% 
  full_join(., cluster_res_list1$NEMO$clust.res %>%
              rename(NEMO = clust),
            by = "samID"
            ) %>% 
  full_join(., cluster_res_list1$ConsensusClustering$clust.res %>%
              rename(ConsensusClustering = clust),
            by = "samID"
  ) %>% 
  full_join(., cluster_res_list1$IntNMF$clust.res %>%
              rename(IntNMF = clust),
            by = "samID"
  ) %>% 
  full_join(., cluster_res_list1$MoCluster$clust.res %>%
              rename(MoCluster = clust),
            by = "samID"
  ) %>% 
  full_join(., cluster_res_list1$iClusterBayes$clust.res %>%
              rename(iClusterBayes = clust),
            by = "samID"
  ) %>% 
  full_join(., cluster_res_list1$BCC$clust.res %>%
              mutate(clust = case_when(
                clust == 1 ~ 2,
                clust == 2 ~ 1
              )) %>% 
              rename(BCC = clust),
            by = "samID"
  ) %>% 
  remove_rownames()

coca_cluster_results <- coca_cluster_results %>% 
  full_join(., annCol, 
            by = c("samID" = "patient_id")) #%>% 
  # filter(BRCA1_carrier != "Unknown" |
  #          BRCA2_carrier != "Unknown")

# cluster_res_list1$COCA$clust.res <- cluster_res_list1$COCA$clust.res %>% 
#   filter(str_detect(samID, paste0(coca_cluster_results$samID, collapse = "|")
#   ))

annCol <- coca_cluster_results %>% 
  column_to_rownames("samID") %>% 
  # filter(BRCA1_carrier != "Unknown" |
  #          BRCA2_carrier != "Unknown") %>% 
  select(iClusterBayes,
         MoCluster, 
         ConsensusClustering, 
         IntNMF,
         NEMO, 
         BCC#, 
    # COCA
         )

annColors <- list(
                  iClusterBayes  = c("1" = "red",
                                     "2"   = "blue",
                                     "Unknown"   = "lightgrey"),
                               
                  MoCluster  = c("1" = "red",
                                         "2"   = "blue",
                                         "Unknown"   = "lightgrey"),
                  ConsensusClustering  = c("1" = "red",
                                         "2"   = "blue",
                                         "Unknown"   = "lightgrey"),
                  IntNMF  = c("1" = "red",
                                         "2"   = "blue",
                                         "Unknown"   = "lightgrey"),
                  NEMO  = c("1" = "red",
                                         "2"   = "blue",
                                         "Unknown"   = "lightgrey"),
                  BCC  = c("1" = "red",
                                     "2"   = "blue",
                                     "Unknown"   = "lightgrey"),
                  COCA  = c("1" = "red",
                                         "2"   = "blue",
                                         "Unknown"   = "lightgrey")
                               )


# cluster_res_list1$COCA$clust.res <- cluster_res_list1$COCA$clust.res %>% 
#   mutate(clust = case_when(
#     clust == 1 ~ "Active",
#     clust == 2 ~ "Repressed"
#   ))

getMoHeatmap(data          = plot_cluster,
             row.title     = c("L1","Alu","LTR"),
             is.binary     = c(F,F,F),
             legend.name   = c("L1","Alu","LTR"),
             clust.res     = cluster_res_list1$COCA$clust.res,
             clust.dend    = NULL, 
             color         = list(c("#0000FF", "white"  , "#FF3C38"), # col.list
                                  c("#0000FF", "white"  , "#FF0000"),
                                  c("#0000FF", "white"  , "#FF0000")),
             width         = 10, 
             height        = 5, 
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             fig.name      = "heatmaps/1 Heatmap for all clusters_10022024_b"
)

# Survival----
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
theme_set(theme_classic())
met_data <- read_rds(paste0(here::here(), "/met_data.rds"))



ggsave_workaround <- function(g){
  survminer:::.build_ggsurvplot(x = g,
                                surv.plot.height = NULL,
                                risk.table.height = NULL,
                                ncensor.plot.height = NULL)}

overall <- ggsurvplot(survfit(Surv(os_time_5year, os_event_5year) ~ coca_RE_cluster,
                   data=met_data),
           # title = "OS Analysis",
           font.main = c(16, "bold", "black"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.legend = c(12, "black"),
           font.tickslab = c(12, "bold", "black"),
           size = 1,
           
           xlab = "Time (months)",
           ylab = "OS (probability)",
           legend = "top",
           legend.title = "COCA cluster", palette = c("red", "blue"),
           legend.labs = c("Active", "Repressed"),
           pval = TRUE,
           conf.int = FALSE,
           # Censor
           censor = TRUE
) %++% guides(colour = guide_legend(ncol = 1))

fig_surv_overall <- ggsave_workaround(overall)

distant_data <- met_data %>%
  filter(stage_cat == "Late")

late_stage <- ggsurvplot(survfit(Surv(os_time_5year, os_event_5year) ~ coca_RE_cluster,
                   data=distant_data),
           # title = "OS Analysis",
           font.main = c(20, "bold", "black"),
           font.x = c(18, "bold", "black"),
           font.y = c(18, "bold", "black"),
           font.legend = c(16, "black"),
           font.tickslab = c(16, "bold", "black"),
           size = 1.5,
           
           xlab = "Time (months)",
           ylab = "OS (probability)",
           legend = "top",
           legend.title = "COCA cluster",
           legend.labs = c("Active", "Repressed"),
           pval = TRUE,
           conf.int = FALSE,
           # Censor
           censor = TRUE
) %++% guides(colour = guide_legend(ncol = 1))
fig_late_stage <- ggsave_workaround(late_stage)


# (overall) + (late_stage) +
#   plot_annotation(tag_levels = "A")

# ggsave("Survival plot.pdf", 
#               width = 7,
#               height = 5, 
#        dpi = 600)


ggarrange(fig_surv_overall,
          fig_late_stage,
          ncol = 2, nrow = 1, 
          labels="AUTO")

library(patchwork)
guide_area() /
  (overall$plot + late_stage$plot)  + 
  plot_layout( guides = 'collect')


ggsave("Survival panel 06192024.pdf", device = cairo_pdf,
       path = here::here(),
       width = 10, height = 5,
       units = c("in"),
       dpi=600,
       bg="white")


# Volcano plot ----
library(tidyverse)
library(REMP)
cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))

## ERV----
remp_res_ERV <- read_rds(paste0(here::here(), "/intermediary data/remp_res_ERV_annotation_07262023.rds"))
remp_res_ERV
annot_ERV_beta_results <- rempB(remp_res_ERV)
annot_ERV_beta_results <- annot_ERV_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_ERV@rowRanges@elementMetadata$RE.Index, .)

ERV_beta_results <-
  annot_ERV_beta_results %>%
  # Merge with annotation, get #chr, start, end
  left_join(., as_tibble(rempAnnot(remp_res_ERV)), by = "Index")

ERV_beta_results <- ERV_beta_results %>% 
  mutate(symbol = coalesce(InNM.symbol, InNR.symbol, InTSS.symbol, 
                           In5UTR.symbol, InCDS.symbol, InExon.symbol, 
                           In3UTR.symbol)) 
ERV_beta_results %>% 
  select(seqnames, strand, repFamily, symbol) %>% 
  tbl_summary(sort = everything() ~ "frequency")


ERV_beta_results_1 <- ERV_beta_results %>% 
  unite(Index, c(Index, seqnames:repName, symbol), sep = "; ") %>% 
  column_to_rownames("Index") %>% 
  select(starts_with("X")) %>% 
  t() %>% as_tibble(rownames = "patient_id") %>% 
  full_join(., cluster_res_list1[["COCA"]][["clust.res"]], 
            by = c("patient_id" = "samID")) %>% 
  select(patient_id, COCA = clust, everything())

# pvalue <- ERV_beta_results_1 %>% 
#   select(-patient_id) %>% 
#   # column_to_rownames("COCA") %>% 
#   t() %>% as_tibble(rownames = "Index") %>% 
#   rename(`1` = "COCA_1", `2` = "COCA_2")
#   summarize(across(everything(), ~ wilcox.test(x, alternative = "two.sided", na.rm = TRUE)))

ERV_beta_results_2 <- ERV_beta_results_1 %>% 
  group_by(COCA) %>% 
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% 
  column_to_rownames("COCA") %>% 
  select(-patient_id) %>% 
  t() %>% as_tibble(rownames = "Index") %>% 
  rename(`1` = "COCA_1", `2` = "COCA_2")

ERV_beta_results_3 <- ERV_beta_results_2 %>% 
  mutate(fold_change = COCA_2 - COCA_1 * 100) %>% 
  arrange(fold_change) %>% 
  dplyr::slice(1 : 20) %>% 
  ggplot(aes(x = fold_change, y = fold_change))+
  geom_point()







# coca_pval <- data.frame(matrix(nrow=0, ncol=1))
# 
# for(i in colnames(ERV_beta_results_1 %>% select(starts_with("ERV")))
#     ) {
#   
#   # name <- i
#   df <- ERV_beta_results_1 %>% select(COCA, all_of(i))
#   # print(df[1:2,])
#   COCA_1 <- df %>% filter(COCA == 1) %>% `colnames<-`(c("COCA", "ERV"))
#   COCA_2 <- df %>% filter(COCA == 2) %>% `colnames<-`(c("COCA", "ERV"))
#   # print(COCA_1)
#   # print(COCA_2)
#   
#   result <- wilcox.test(c(COCA_1$ERV), 
#                         COCA_2$ERV)
#   # print(i)
#   # print(result$p.value)
#   coca_pval[i,"pvalue"] <- result$p.value
#   # print(coca_pval)
# 
# }
# coca_pval1 <- coca_pval %>% 
#   rownames_to_column("Index") %>% 
#   select(Index, pvalue)
# write_rds(coca_pval1, "intermediary data/calculated pval for ERVs.rds")
remp_res_ERV <- read_rds(paste0(here::here(), "/intermediary data/calculated pval for ERVs.rds"))

a <- ERV_beta_results_2 %>% 
  full_join(., coca_pval1, by = "Index") %>% 
  mutate(fold_change = COCA_1 - COCA_2) %>% 
  mutate(fold_change_log2 = log2(fold_change)) %>% 
  mutate(symbol = str_extract(Index, '\\w+$')
  ) %>% 
  filter(symbol != "NA") 
write_rds(a %>% filter(pvalue <= 0.05), "intermediary data/calculated pval <= 0.05 for ERVs with gene.rds")

library(ggrepel)
library(plotly)
# a %>% 
#   # arrange(fold_change) %>% 
#   # dplyr::slice(1 : 500) %>% 
#   mutate(text = paste("Gene: ", symbol, sep="")) %>% 
#   ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
#   geom_point()+
#   geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(0.05), 
#              col = "gray", linetype = 'dashed')+
#   geom_text_repel(max.overlaps = Inf)

a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')+
  coord_cartesian(ylim = c(1.3, 3), xlim = c(0.01, 0.1))+
  geom_text_repel(max.overlaps = Inf)

a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  # filter(fold_change < 0.1 & pvalue < 0.01) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')+
  coord_cartesian(ylim = c(2.5, 10), xlim = c(-0.3, -0.15)
                  )+
  geom_text_repel(max.overlaps = Inf)

plot <- a %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), text=text))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')
plot <- ggplotly(plot, tooltip="text")
plot

## L1----
remp_res_L1 <- read_rds(paste0(here::here(), "/intermediary data/remp_res_L1_annotation_07262023.rds"))
remp_res_L1
annot_L1_beta_results <- rempB(remp_res_L1)
annot_L1_beta_results <- annot_L1_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_L1@rowRanges@elementMetadata$RE.Index, .)

L1_beta_results <-
  annot_L1_beta_results %>%
  # Merge with annotation, get #chr, start, end
  left_join(., as_tibble(rempAnnot(remp_res_L1)), by = "Index")

L1_beta_results <- L1_beta_results %>% 
  mutate(symbol = coalesce(InNM.symbol, InNR.symbol, InTSS.symbol, 
                           In5UTR.symbol, InCDS.symbol, InExon.symbol, 
                           In3UTR.symbol)) 
L1_beta_results %>% 
  select(seqnames, strand, repFamily, symbol) %>% 
  tbl_summary(sort = everything() ~ "frequency")


L1_beta_results_1 <- L1_beta_results %>% 
  unite(Index, c(Index, seqnames:repName, symbol), sep = "; ") %>% 
  column_to_rownames("Index") %>% 
  select(starts_with("X")) %>% 
  t() %>% as_tibble(rownames = "patient_id") %>% 
  full_join(., cluster_res_list1[["COCA"]][["clust.res"]], 
            by = c("patient_id" = "samID")) %>% 
  select(patient_id, COCA = clust, everything())

# pvalue <- L1_beta_results_1 %>% 
#   select(-patient_id) %>% 
#   # column_to_rownames("COCA") %>% 
#   t() %>% as_tibble(rownames = "Index") %>% 
#   rename(`1` = "COCA_1", `2` = "COCA_2")
#   summarize(across(everything(), ~ wilcox.test(x, alternative = "two.sided", na.rm = TRUE)))

L1_beta_results_2 <- L1_beta_results_1 %>% 
  group_by(COCA) %>% 
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% 
  column_to_rownames("COCA") %>% 
  select(-patient_id) %>% 
  t() %>% as_tibble(rownames = "Index") %>% 
  rename("COCA_1" = `1`,"COCA_2" = `2`)

L1_beta_results_3 <- L1_beta_results_2 %>% 
  mutate(fold_change = COCA_2 - COCA_1 * 100) %>% 
  arrange(fold_change) %>% 
  dplyr::slice(1 : 20) %>% 
  ggplot(aes(x = fold_change, y = fold_change))+
  geom_point()

coca_pval <- data.frame(matrix(nrow=0, ncol=1))

for(i in colnames(L1_beta_results_1 %>% select(starts_with("L1")))
) {
  
  # name <- i
  df <- L1_beta_results_1 %>% select(COCA, all_of(i))
  # print(df[1:2,])
  COCA_1 <- df %>% filter(COCA == 1) %>% `colnames<-`(c("COCA", "L1"))
  COCA_2 <- df %>% filter(COCA == 2) %>% `colnames<-`(c("COCA", "L1"))
  # print(COCA_1)
  # print(COCA_2)
  
  result <- wilcox.test(c(COCA_1$L1), 
                        COCA_2$L1)
  # print(i)
  # print(result$p.value)
  coca_pval[i,"pvalue"] <- result$p.value
  # print(coca_pval)
  
}
coca_pval1 <- coca_pval %>% 
  rownames_to_column("Index") %>% 
  select(Index, pvalue)
write_rds(coca_pval1, "intermediary data/calculated pval for L1s.rds")
coca_pval1 <- read_rds(paste0(here::here(), "/intermediary data/calculated pval for L1s.rds"))

a <- L1_beta_results_2 %>% 
  full_join(., coca_pval1, by = "Index") %>% 
  mutate(fold_change = COCA_1 - COCA_2) %>% 
  mutate(fold_change_log2 = log2(fold_change)) %>% 
  mutate(symbol = str_extract(Index, '\\w+$')
  ) %>% 
  filter(symbol != "NA") 
# write_rds(a %>% filter(pvalue <= 0.05), "intermediary data/calculated pval <= 0.05 for L1s with gene.rds")

library(ggrepel)
library(plotly)
# a %>% 
#   # arrange(fold_change) %>% 
#   # dplyr::slice(1 : 500) %>% 
#   mutate(text = paste("Gene: ", symbol, sep="")) %>% 
#   ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
#   geom_point()+
#   geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(0.05), 
#              col = "gray", linetype = 'dashed')+
#   geom_text_repel(max.overlaps = Inf)

a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')+
  coord_cartesian(ylim = c(1.3, 3), xlim = c(0.01, 0.1))+
  geom_text_repel(max.overlaps = Inf)

a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  # filter(fold_change < 0.1 & pvalue < 0.01) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')+
  coord_cartesian(ylim = c(2.5, 10), xlim = c(-0.3, -0.15)
  )+
  geom_text_repel(max.overlaps = Inf)

plot <- a %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), text=text))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')
plot <- ggplotly(plot, tooltip="text")
plot

## Alu----
remp_res_Alu <- read_rds(paste0(here::here(), "/intermediary data/remp_res_Alu_annotation_07262023.rds"))
remp_res_Alu
annot_Alu_beta_results <- rempB(remp_res_Alu)
annot_Alu_beta_results <- annot_Alu_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_Alu@rowRanges@elementMetadata$RE.Index, .)

Alu_beta_results <-
  annot_Alu_beta_results %>%
  # Merge with annotation, get #chr, start, end
  left_join(., as_tibble(rempAnnot(remp_res_Alu)), by = "Index")

Alu_beta_results <- Alu_beta_results %>% 
  mutate(symbol = coalesce(InNM.symbol, InNR.symbol, InTSS.symbol, 
                           In5UTR.symbol, InCDS.symbol, InExon.symbol, 
                           In3UTR.symbol)) 
Alu_beta_results %>% 
  select(seqnames, strand, repFamily, symbol) %>% 
  tbl_summary(sort = everything() ~ "frequency")


Alu_beta_results_1 <- Alu_beta_results %>% 
  unite(Index, c(Index, seqnames:repName, symbol), sep = "; ") %>% 
  column_to_rownames("Index") %>% 
  select(starts_with("X")) %>% 
  t() %>% as_tibble(rownames = "patient_id") %>% 
  full_join(., cluster_res_list1[["COCA"]][["clust.res"]], 
            by = c("patient_id" = "samID")) %>% 
  select(patient_id, COCA = clust, everything())

# pvalue <- Alu_beta_results_1 %>% 
#   select(-patient_id) %>% 
#   # column_to_rownames("COCA") %>% 
#   t() %>% as_tibble(rownames = "Index") %>% 
#   rename(`1` = "COCA_1", `2` = "COCA_2")
#   summarize(across(everything(), ~ wilcox.test(x, alternative = "two.sided", na.rm = TRUE)))

Alu_beta_results_2 <- Alu_beta_results_1 %>% 
  group_by(COCA) %>% 
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% 
  column_to_rownames("COCA") %>% 
  select(-patient_id) %>% 
  t() %>% as_tibble(rownames = "Index") %>% 
  rename(`1` = "COCA_1", `2` = "COCA_2")

Alu_beta_results_3 <- Alu_beta_results_2 %>% 
  mutate(fold_change = COCA_2 - COCA_1 * 100) %>% 
  arrange(fold_change) %>% 
  dplyr::slice(1 : 20) %>% 
  ggplot(aes(x = fold_change, y = fold_change))+
  geom_point()

# coca_pval <- data.frame(matrix(nrow=0, ncol=1))
# 
# for(i in colnames(Alu_beta_results_1 %>% select(starts_with("Alu")))
# ) {
#   
#   # name <- i
#   df <- Alu_beta_results_1 %>% select(COCA, all_of(i))
#   # print(df[1:2,])
#   COCA_1 <- df %>% filter(COCA == 1) %>% `colnames<-`(c("COCA", "Alu"))
#   COCA_2 <- df %>% filter(COCA == 2) %>% `colnames<-`(c("COCA", "Alu"))
#   # print(COCA_1)
#   # print(COCA_2)
#   
#   result <- wilcox.test(c(COCA_1$Alu), 
#                         COCA_2$Alu)
#   # print(i)
#   # print(result$p.value)
#   coca_pval[i,"pvalue"] <- result$p.value
#   # print(coca_pval)
#   
# }
# coca_pval1 <- coca_pval %>% 
#   rownames_to_column("Index") %>% 
#   select(Index, pvalue)
# write_rds(coca_pval1, "intermediary data/calculated pval for Alus.rds")
coca_pval1 <- read_rds(paste0(here::here(), "/intermediary data/calculated pval for Alus.rds"))

a <- Alu_beta_results_2 %>% 
  full_join(., coca_pval1, by = "Index") %>% 
  mutate(fold_change = COCA_1 - COCA_2) %>% 
  mutate(fold_change_log2 = log2(fold_change)) %>% 
  mutate(symbol = str_extract(Index, '\\w+$')
  ) %>% 
  filter(symbol != "NA") 

library(ggrepel)
library(plotly)
a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')#+
  # geom_text_repel(max.overlaps = Inf)

a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')+
  coord_cartesian(ylim = c(1.3, 3), xlim = c(0.01, 0.1))+
  geom_text_repel(max.overlaps = Inf)

a %>% 
  # arrange(fold_change) %>% 
  # dplyr::slice(1 : 500) %>% 
  # filter(fold_change < 0.1 & pvalue < 0.01) %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), label = symbol))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')+
  coord_cartesian(ylim = c(2.5, 10), xlim = c(-0.3, -0.15)
  )+
  geom_text_repel(max.overlaps = Inf)

plot <- a %>% 
  mutate(text = paste("Gene: ", symbol, sep="")) %>% 
  ggplot(aes(x = (fold_change), y = -log10(pvalue), text=text))+
  geom_point()+
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), 
             col = "gray", linetype = 'dashed')
plot <- ggplotly(plot, tooltip="text")
plot
# END volcano plot----

# Forest plot----
library(tidyverse)
library(survival)
theme_set(theme_classic())
met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
# forest_p <- coxph(Surv(time = met_data$os_time_5year,
#                        event = met_data$os_event_5year) ~ coca_RE_cluster + refage + stage_cat + 
#                     BRCA1_carrier + BRCA2_carrier,
#                   data = met_data)  %>%
#   tbl_regression(exponentiate = TRUE) %>%
#   as_tibble() %>% 
#   `colnames<-`(c("Characteristic", "HR", "CI95", "pvalue")) %>% 
#   fill(HR, CI95, pvalue, .direction = "up") %>% 
#   separate(CI95, into = c("lower", "upper"), sep = ", ") %>% 
#   filter(Characteristic %in% c("COCA cluster",
#                                "Age at diagnosis",
#                                "Stage",
#                                "BRCA1_carrier",
#                                "BRCA2_carrier") 
#   ) %>% 
#   mutate(across(everything(), ~ str_replace(., "_", " "))) %>% 
#   mutate(across(c(HR, lower, upper, pvalue), ~ as.numeric(.)))

forest_p <- tibble::tribble(
  ~Characteristic, ~`N.(Deaths)`,      ~`HR.(95%.CI)`, ~`p-value`, 
  "Model 1",   "175 (105)", "2.01 (1.01, 3.98)",      0.045,    
  "Model 2",   "125 (75)", "2.41 (1.04, 5.59)",       0.04,
  "Model 4",   "135 (83)", "2.91 (1.17, 7.20)",      0.021,
  "Model 3",   "101 (69)", "2.87 (1.15, 7.16)",      0.024 # restricted to early stage
  ) %>% 
  separate_wider_delim(cols = `HR.(95%.CI)`, delim = " (",
                       names = c("HR", "lower"), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = FALSE) %>% 
  separate_wider_delim(cols = lower, delim = ", ",
                       names = c("lower", "upper"), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  mutate(upper = str_remove(upper, "\\)")) %>% 
  mutate(Characteristic = factor(Characteristic, levels = c("Model 4",
                                                            "Model 3",
                                                            "Model 2",
                                                            "Model 1"
                                                            ))) %>% 
  mutate_at(c("HR", "lower", "upper"), ~ as.numeric(.)) %>% 
  mutate(text = paste0(`HR.(95%.CI)`, ", p-value = ", `p-value`))


forest_p <- ggplot(data = forest_p, 
       aes(y = Characteristic, x = HR, xmin = lower, xmax = upper), log10="y") +
  geom_point() +
  geom_errorbarh(height = .1) +
  geom_vline(xintercept = 1, color='black', linetype='dashed', alpha=.5
  )+ 
  geom_text(aes(y = Characteristic, label = text, x = after_stat(xmax)-1), vjust = -1.5)+
  theme_classic()+
  scale_x_continuous(transform = "log10",
                     breaks = c(0.8, 1, 3, 10),
                     labels = c(0.8, 1, 3, 10),
                     limits = c(0.8,10),
                     expand = c(-0.0001, 0.1)
                     )+
  labs(x = "HR (95% CI)", y = ""#, caption = "Model 3 = Adjusting for BRCA1 and BRCA2 mutations"
       )
forest_p
ggsave("Figure forest plot_10032024.pdf",
       width = 7,
       height = 4, 
       dpi = 600)

# forest_p <-
#   tbl_all_mutations %>%
#   as_tibble() %>% 
#   `colnames<-`(c("Group", "Characteristic", "remove_n1", "remove_event1", 
#                  "HR_tumor", "CI95_tumor", "pvalue_tumor", "remove_n2", "remove_event1", 
#                  "HR_germline", "CI95_germline", "pvalue_germline")) %>% 
#   select(-contains("remove_")) %>% 
#   fill(Group, .direction = "down") %>% 
#   fill(HR_tumor, CI95_tumor, pvalue_tumor, HR_germline, CI95_germline, pvalue_germline, .direction = "up") %>% 
#   mutate(across(everything(), ~ str_remove_all(., "__"))) %>% 
#   filter(Characteristic %in% c("coca_RE_cluster",
#                  # "refage",
#                  # "stage_cat",
#                  "tumor_mutation_BRCA1",
#                  "tumor_mutation_BRCA2",
#                  "tumor_mutation_TP53") 
#          ) %>% 
#   mutate(name= case_when(
#     str_detect(Characteristic, "tumor_mutation")      ~ "mutation",
#     TRUE                                                ~ Characteristic
#   )) %>% 
#   # unite(Characteristic, c(Group, Characteristic), sep = "_", remove = TRUE) %>% 
#   select(-name) %>% 
#   pivot_longer(cols = -c(Group, Characteristic),
#                names_pattern = "(.*)_(.*)$", 
#                names_to = c("estimate", "name"))

# ggplot(data=forest_p, aes(y=Characteristic, x=OR, xmin=lower, xmax=upper)) +
#   geom_point() + 
#   geom_errorbarh(height=.1) +
#   geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
#   theme_classic()+
#   scale_x_log10(breaks = c(0, 0.25, 0.5, 1, 2.5, 5, 10, 40))+
#   labs(y="")


# Figure 2 KM + Forest ----

library(patchwork)

(free(overall$plot) + free(forest_p)) + 
  plot_annotation(tag_levels = "A")+
  plot_layout(
              # heights = c(3, 6),
              widths = c(3, 3.5)
              )

ggsave("Figure 2 KM forest plot_10032024.pdf", device = cairo_pdf,
       path = here::here(),
       width = 11, height = 5,
       units = c("in"),
       dpi=600,
       bg="white")


