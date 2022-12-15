# Import packages
library(tidyverse)
# library(minfi)
library(REMP)
# HELP chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://bioconductor.org/packages/release/bioc/vignettes/REMP/inst/doc/REMP.pdf

############################################################# Load data

# load("/Volumes/Peres_Research/Methylomics R01/Data from Lucas/Unfiltered betas.RData")
# betas_sesame R01C01 column and cg ids rowname

load("/Volumes/Peres_Research/Methylomics R01/Data from Lucas/cleaned_07082022.rda")
# beta_clean6 R01C01 column and cg ids rowname
# and phenoclean

beta_data <- betas_clean6 # Looks already normalized
summary(beta_data)
# 1 Groommethylation data

grooMethy(beta_data)


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)
genomic_location_matrix <- getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)

groomed_data <- grooMethy(beta_data, Seq.GR = genomic_location_matrix)

# 2 Prepare annotation data

remparcel <- initREMP(arrayType = "Sequencing", Seq.GR = genomic_location_matrix,
                      REtype = "Alu", 
                      annotation.source = "AH", genome = "hg19", 
                      ncore = 1)

remparcel
saveParcel(remparcel)

# 3 Runprediction

remp.res <- remp(groomed_data, REtype = 'Alu',
                 parcel = remparcel, ncore = 1, seed = 777)
remp.res
# Prediction results can be obtained by accessors: 
#Predicted RE-CpG methylation value (Beta value) 
rempB(remp.res)
#Predicted RE-CpG methylation value (M value) 
rempM(remp.res)
  
  
  
  
  
  
  
  









