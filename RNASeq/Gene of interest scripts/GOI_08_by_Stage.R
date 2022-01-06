
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")

# Prepare data
plot_df = data.frame(Sample_ID = colnames(RNASeq.QN.LOG2), Gene_expression = RNASeq.QN.LOG2[which(rownames(RNASeq.QN.LOG2) == "DUX4"),],
                     Stage = NA)

