
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Set parameters
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_nonsynonymous_filter_all_samples.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")
excluded_df = excluded_df[which(excluded_df$Reason.excluded %in% exclude),]

# Analysis
MAF_df = finalMafFiltered
MAF_df = MAF_df[-which(MAF_df$Patient_ID %in% excluded_df$Patient_ID),]

MAF_df$Patient_ID = paste("SER-SILU-CC-P0", MAF_df$Patient_ID, sep = "")
MAF_df$Sample_ID = paste("SER-SILU-CC-P0", MAF_df$Sample_ID, sep = "")
MAF_df$Sample_ID = gsub("T", "-PT-01-B-02", MAF_df$Sample_ID)

# Little trick to relabel liver metastasis consistently
MAF_df$Sample_ID = gsub("LM2", "x", MAF_df$Sample_ID)
MAF_df$Sample_ID = gsub("LM", "LM1", MAF_df$Sample_ID)
MAF_df$Sample_ID = gsub("LM11", "LM1", MAF_df$Sample_ID)
MAF_df$Sample_ID = gsub("LM1", "-LM-01-B-02", MAF_df$Sample_ID)
MAF_df$Sample_ID = gsub("x", "-LM-02-B-02", MAF_df$Sample_ID)

MAF_df = MAF_df[-grep("LM-", MAF_df$Sample_ID),]

unique(MAF_df$Sample_ID)

MAF_df$Sample_ID = substring(MAF_df$Sample_ID, 1, 23)

write.table(MAF_df, file = "./Processed_Data/Shared_Data/For cbioportal/2405_2021_SER-SILU-WES_Mutect_Strelka_for_cbioportal.tsv",
          row.names = FALSE, sep = "\t")
