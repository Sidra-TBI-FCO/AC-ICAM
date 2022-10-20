
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# BiocManager::install("survcomp")
ipak(c("ggplot2", "survival", "survivalAnalysis", "data.table", 
       "corrplot", "creditmodel", "dplyr", "glmnet", "doParallel", "extrafont",
       "extrafontdb", "stringr", "stringi", "factoextra", "ggfortify", "survminer",
       "gbm", "survcomp"))

source("./R code/Microbiome/Risk_model_glmnet/get_functions.R")

# Load data
load("./Processed_Data/Microbiome/All_Input_Data_for_Risk_Model/full_genus_matrix_not_filtered_42_tumor_samples_only.Rdata")
load("./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_AC_ICAM_246_for_glmnet.Rdata")
load("../NGS_data/Processed_Data/Survival Data/clinical_data_348_patients.Rdata")

#Perform validation on 42 additional samples
subset_42_Genus_full_abundance <- Genus_full_abundance
subset_42_Genus_subset_abundance <- subset_42_Genus_full_abundance[rownames(revised_Genus_full_abundance),]
revised_colnames <- str_replace(colnames(subset_42_Genus_subset_abundance),pattern = "T","")
colnames(subset_42_Genus_subset_abundance) <- revised_colnames
subset_42_clinical_data <- clinical_data[clinical_data$Patient_ID %in% colnames(subset_42_Genus_subset_abundance),]

#Build the survival dataframe
subset_42_survival_df <- as.data.frame(cbind(subset_42_clinical_data[,c("Patient_ID","OS.Status","OS.Time")], t(subset_42_Genus_subset_abundance)))

#Make the survival df with event == 1 (death) and alive == 0 or censored
rev_subset_42_survival_df <- subset_42_survival_df
colnames(rev_subset_42_survival_df)[c(1:3)] <- c("samples","status","time")
rev_subset_42_survival_df$time <- as.numeric(as.vector(rev_subset_42_survival_df$time))
rev_subset_42_survival_df[rev_subset_42_survival_df$status=="Alive",]$status <- 0
rev_subset_42_survival_df[rev_subset_42_survival_df$status=="Dead",]$status <- 1
rev_subset_42_survival_df$status <- as.numeric(as.vector(rev_subset_42_survival_df$status))
for (i in 2:ncol(rev_subset_42_survival_df)){
  rev_subset_42_survival_df[,i] <- as.numeric(as.vector(rev_subset_42_survival_df[,i]))
}
rev_subset_42_survival_df$time <- rev_subset_42_survival_df$time/365

#Normalize the 42 sample data
# mean_vector and sd_vector from the 246 matrix
for (i in c(1:(ncol(rev_subset_42_survival_df)-3))){
  rev_subset_42_survival_df[,i+3] <- (rev_subset_42_survival_df[,i+3]-mean_vector[i])/sd_vector[i]
}

save(rev_subset_42_survival_df, colnames_of_interest, y_cols,
     file = "./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_AC_ICAM_42_for_glmnet.Rdata")



