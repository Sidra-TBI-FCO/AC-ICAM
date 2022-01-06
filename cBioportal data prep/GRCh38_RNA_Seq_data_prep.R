
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))

library(stringr)

# Set parameters
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.Complete.dataset.EDAseq.QN.HPC.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# subset
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, -which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% 
                                           excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)])]
dim(RNASeq.QN.LOG2)

colnames(RNASeq.QN.LOG2) = paste("SER-SILU-CC-P0", colnames(RNASeq.QN.LOG2), sep = "")
colnames(RNASeq.QN.LOG2) = gsub("T_P", "-PT-01-A-01", colnames(RNASeq.QN.LOG2))
colnames(RNASeq.QN.LOG2) = gsub("LM1", "-LM-01-A-01", colnames(RNASeq.QN.LOG2))
colnames(RNASeq.QN.LOG2) = gsub("LM2", "-LM-02-A-01", colnames(RNASeq.QN.LOG2))
colnames(RNASeq.QN.LOG2) = gsub("T_B", "-PT-02-A-01", colnames(RNASeq.QN.LOG2))
colnames(RNASeq.QN.LOG2) = gsub("T_C", "-PT-03-A-01", colnames(RNASeq.QN.LOG2))

# Remove liver metastasis samples
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, -grep("LM-", colnames(RNASeq.QN.LOG2))]
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, -grep("-PT-02", colnames(RNASeq.QN.LOG2))]
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, -grep("-PT-03", colnames(RNASeq.QN.LOG2))]

dim(RNASeq.QN.LOG2)

write.csv(RNASeq.QN.LOG2, file = "./Processed_Data/Shared_Data/For cbioportal/2405_2021_SER-SILU-RNASeq_QN_LOG2_for_cbioportal.csv") 

# z score will be calculated by cbioportal importer tool!
# z score matrix
#RNASeq.QN.LOG2_z = RNASeq.QN.LOG2 
#for(j in 1: nrow(RNASeq.QN.LOG2_z))  {
# RNASeq.QN.LOG2_z[j,] = (RNASeq.QN.LOG2[j,]-mean(RNASeq.QN.LOG2[j,]))/sd(RNASeq.QN.LOG2[j,]) # z-score the enrichment matrix
#}
#RNASeq.QN.LOG2_z[which(is.na(RNASeq.QN.LOG2_z))] = 0
#min(RNASeq.QN.LOG2_z)
#max(RNASeq.QN.LOG2_z)

#write.csv(RNASeq.QN.LOG2_z, file = "./Processed_Data/Shared_Data/For cbioportal/SER-SILU-RNASeq_Z_SCORE_matrix_GRCh37_for_cbioportal.csv") 
