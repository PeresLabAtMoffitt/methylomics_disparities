# Import library
library(tidyverse)
library(REMP)
library(gtsummary)
library(patchwork)
theme_gtsummary_compact()
theme_set(theme_classic())

# Figure S1 patients methylation pattern ----
# Load data
# LTR
remp_res_ERV <- read_rds(paste0(here::here(), "/intermediary data/remp_res_ERV_annotation_07262023.rds"))
remplot(remp_res_ERV, main = "ERV methylation", col = "blue")
annot_ERV_beta_results <- rempB(remp_res_ERV)
map_df <- annot_ERV_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_ERV@rowRanges@elementMetadata$RE.Index, .)

LTR <- map_df %>% 
  pivot_longer(cols = - Index) %>% 
  ggplot(aes(x= value, group= name))+
  geom_density(color = "black", linewidth = 0.1)+
  labs(x = "Methylation value (beta)")+
  xlim(0,1)+
  theme(legend.position = "none")

# L1
remp_res_L1 <- read_rds(paste0(here::here(), "/intermediary data/remp_res_L1_annotation_07262023.rds"))
remplot(remp_res_L1, main = "L1 methylation", col = "blue")
annot_L1_beta_results <- rempB(remp_res_L1)
map_df <- annot_L1_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_L1@rowRanges@elementMetadata$RE.Index, .)

L1 <- map_df %>% 
  pivot_longer(cols = - Index) %>% 
  ggplot(aes(x= value, group= name))+
  geom_density(color = "black", linewidth = 0.1)+
  labs(x = "Methylation value (beta)")+
  xlim(0,1)+
  theme(legend.position = "none")

# Alu
remp_res_Alu <- read_rds(paste0(here::here(), "/intermediary data/remp_res_Alu_annotation_07262023.rds"))
remplot(remp_res_Alu, main = "Alu methylation", col = "blue")
annot_Alu_beta_results <- rempB(remp_res_Alu)
map_df <- annot_Alu_beta_results %>% as_tibble() %>%
  cbind(Index = remp_res_Alu@rowRanges@elementMetadata$RE.Index, .)

ALU <- map_df %>% 
  pivot_longer(cols = - Index) %>% 
  ggplot(aes(x= value, group= name))+
  geom_density(color = "black", linewidth = 0.1)+
  labs(x = "Methylation value (beta)")+
  xlim(0,1)+
  theme(legend.position = "none")


LTR / L1 / ALU +
  plot_annotation(tag_levels = "A")+
  plot_layout(axes = "collect"#,
              # axis_titles = "collect"
              )

ggsave("Figure S1 patients methylation pattern.pdf",
       width = 4,
       height = 10, 
       dpi = 600)


# Heatmap----
library(tidyverse)
library(MOVICS)
theme_set(theme_classic())
met_data <- read_rds(paste0(here::here(), "/met_data.rds"))
cluster_res_list1 <- read_rds(paste0(here::here(), "/clustering/MOVICS/movics_cluster_res_list1_04112024.rds"))
plot_cluster <- read_rds(paste0(here::here(), "/plot_cluster.rds"))

# BRCA 1 and 2
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
                    ))

annCol <- annCol %>% 
  column_to_rownames("patient_id") %>% 
  filter(BRCA1_carrier != "Unknown" |
           BRCA2_carrier != "Unknown")

annColors <- list(BRCA1_carrier  = c("Yes" = "red",
                                               "No"   = "blue",
                                               "Unknown"   = "lightgrey"),
                  BRCA2_carrier  = c("Yes" = "red",
                                            "No"   = "blue",
                                            "Unknown"   = "lightgrey")
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
             annColors     = annColors#, # annotation color
             # fig.name      = "Heatmap of COCA_with_limited_BRCA_05102024"
             )

# BRCA 1 and 2 - overall/germline/tumor
annCol <- met_data %>% select(patient_id, 
                              BRCA1_carrier, BRCA2_carrier,
                              germline_mutation_BRCA1, tumor_mutation_BRCA1, 
                              germline_mutation_BRCA2, tumor_mutation_BRCA2
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
  ))

annCol <- annCol %>% 
  column_to_rownames("patient_id") %>% 
  filter(BRCA1_carrier != "Unknown" |
           BRCA2_carrier != "Unknown")

annColors <- list(BRCA1_carrier  = c("Yes" = "red",
                                     "No"   = "blue",
                                     "Unknown"   = "lightgrey"),
                  BRCA2_carrier  = c("Yes" = "red",
                                     "No"   = "blue",
                                     "Unknown"   = "lightgrey"),
                  germline_mutation_BRCA1  = c("Yes" = "red",
                                               "No"   = "blue",
                                               "Unknown"   = "lightgrey"),
                  tumor_mutation_BRCA1  = c("Yes" = "red",
                                            "No"   = "blue",
                                            "Unknown"   = "lightgrey"),
                  germline_mutation_BRCA2  = c("Yes" = "red",
                                               "No"   = "blue",
                                               "Unknown"   = "lightgrey"),
                  tumor_mutation_BRCA2  = c("Yes" = "red",
                                            "No"   = "blue",
                                            "Unknown"   = "lightgrey")
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
             fig.name      = "Heatmap of COCA_with_limited_BRCA_6bars_05102024"
)

# BRCA 1 and 2 - overall/germline/tumor
annCol <- met_data %>% select(patient_id, 
                              germline_mutation_BRCA1, tumor_mutation_BRCA1, 
                              germline_mutation_BRCA2, tumor_mutation_BRCA2
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
  ))

annCol <- annCol %>% 
  column_to_rownames("patient_id") %>% 
  filter(BRCA1_carrier != "Unknown" |
           BRCA2_carrier != "Unknown")

annColors <- list(germline_mutation_BRCA1  = c("Yes" = "red",
                                               "No"   = "blue",
                                               "Unknown"   = "lightgrey"),
                  tumor_mutation_BRCA1  = c("Yes" = "red",
                                            "No"   = "blue",
                                            "Unknown"   = "lightgrey"),
                  germline_mutation_BRCA2  = c("Yes" = "red",
                                               "No"   = "blue",
                                               "Unknown"   = "lightgrey"),
                  tumor_mutation_BRCA2  = c("Yes" = "red",
                                            "No"   = "blue",
                                            "Unknown"   = "lightgrey")
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
             fig.name      = "Heatmap of COCA_with_limited_BRCA_4bars_05102024"
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
           font.main = c(20, "bold", "black"),
           font.x = c(18, "bold", "black"),
           font.y = c(18, "bold", "black"),
           font.legend = c(16, "black"),
           font.tickslab = c(16, "bold", "black"),
           size = 1.5,
           
           xlab = "Time (days)",
           ylab = "OS (probability)",
           legend = "top",
           legend.title = "COCA cluster",
           legend.labs = c("1", "2"),
           pval = TRUE,
           conf.int = FALSE,
           # Censor
           censor = TRUE
) + guides(colour = guide_legend(ncol = 1))

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
           
           xlab = "Time (days)",
           ylab = "OS (probability)",
           legend = "top",
           legend.title = "COCA cluster",
           legend.labs = c("1", "2"),
           pval = TRUE,
           conf.int = FALSE,
           # Censor
           censor = TRUE
) + guides(colour = guide_legend(ncol = 1))
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


