######## REMP dataset and tests

# Import packages
library(tidyverse)
library(REMP)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
############################################################# Load data


REMP::Alu.hg38.demo
Alu.hg38.demo <- Alu.hg38.demo

# 1. Groomm methylation data
# Get GM12878 methylation data (450k array)
GM12878_450k <- getGM12878('450k')
GM12878_450k@assays@data@listData
GM12878_450k <- grooMethy(GM12878_450k)
GM12878_450k

# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# 2. Prepare annotation data
data(Alu.hg19.demo)
remparcel <- initREMP(arrayType = "450k", 
                      REtype = "Alu",
                      annotation.source = "AH", 
                      genome = "hg19", 
                      RE = Alu.hg19.demo,
                      ncore = 1)
remparcel
saveParcel(remparcel)

# 3. Run prediction
remp.res <- remp(GM12878_450k,
                 REtype = 'Alu',
                 parcel = remparcel, ncore = 1, seed = 777)
remp.res
# Display more detailed information 
details(remp.res)
# Prediction results can be obtained by accessors: 
#Predicted RE-CpG methylation value (Beta value) 
rempB(remp.res)
#Predicted RE-CpG methylation value (M value) 
rempM(remp.res)

#GenomiclocationinformationofthepredictedRE-CpG 
#Functioninheritfromclass'RangedSummarizedExperiment' 
rowRanges(remp.res)
#Standarderror-scaledpermutationimportanceofpredictors 
rempImp(remp.res)
#Retriveseednumberusedforthereesults 
metadata(remp.res)$Seed

# Trimofflessreliablepredictedresults: 
#AnypredictedCpGvalueswithqualityscorelessthan 
#threshold(default=1.7)willbereplacedwithNA. 
#CpGscontainmorethanmissingRate*100%(default=20%) 
#missingrateacrosssampleswillbediscarded. 
remp.res <- rempTrim(remp.res,threshold=1.7,missingRate=0.2) 
details(remp.res)

# (Optional) Aggregate the predicted methylation of CpGs in RE by averaging them to obtain the RE-specific methylation level: 
remp.res <- rempAggregate(remp.res, NCpG = 2)
details(remp.res)

# Aggregating CpGs in the same RE for RE-level methylation data is beneficial because 1) it greatly reduces the 
# data dimension for downstream analysis and 2) it may produce more robust RE methylation estimation. Note that 
# by default, RE with 2 or more predicted CpG sites will be aggregated. Therefore, the downside of doing this is the 
# reduced coverage of RE. The assumption of doing this is the CpG methylation level within each RE are similar.
# To add genomic regions annotation of the predicted REs: 
# By default gene symbol annotation will be added 
remp.res<-decodeAnnot(remp.res) 
rempAnnot(remp.res)

# 4. Plot prediction
remplot(remp.res, main = "Alu methylation (GM12878)", col = "blue")








