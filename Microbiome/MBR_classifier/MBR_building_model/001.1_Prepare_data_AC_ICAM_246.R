
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

loadfonts()
registerDoParallel(cores=10)
source("./R code/Microbiome/Risk_model_glmnet/get_functions.R")

#Get the data for AC-ICAM
load("./Processed_Data/Microbiome/All_Input_Data_for_Risk_Model/AC-ICAM_filtered_microbiome_10_percent_abundance_yes_0.01_T_samples.Rdata")
dataset <- "AC-ICAM"

#Perform preprocessing to only have the tumor samples
normal_ids <- grep("N",colnames(Genus_full_abundance))
if (length(normal_ids)>0)
{
  revised_Genus_full_abundance <- Genus_full_abundance[,-normal_ids]
  colnames(revised_Genus_full_abundance) <- str_replace(colnames(revised_Genus_full_abundance),pattern = "T","")
}else{
  revised_Genus_full_abundance <- Genus_full_abundance
  colnames(revised_Genus_full_abundance) <- str_replace(colnames(revised_Genus_full_abundance),pattern = "T","")
}

#Get the correlations between the different genera
corr_matrix <- cor(t(revised_Genus_full_abundance))

#The genera are not heavily correlated (majority of them have |cor| < 0.2)
plot(density(corr_matrix))
hist(corr_matrix, breaks = 20)

#Get the clinical information
load("../NGS_data/Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
subset_clinical_data <- clinical_data[clinical_data$Patient_ID %in% colnames(revised_Genus_full_abundance),]
final_Genus_full_abundance <- revised_Genus_full_abundance[,colnames(revised_Genus_full_abundance) %in% subset_clinical_data$Patient_ID]

colnames(final_Genus_full_abundance) == subset_clinical_data$Patient_ID # check if all is TRUE before cbind in next step

#Build the survival dataframe
survival_df <- as.data.frame(cbind(subset_clinical_data[,c("Patient_ID","OS.Status","OS.Time")], t(final_Genus_full_abundance)))

#Make the survival df with event == 1 (death) and alive == 0 or censored
rev_survival_df <- survival_df
colnames(rev_survival_df)[c(1:3)] <- c("samples","status","time")
rev_survival_df$time <- as.numeric(as.vector(rev_survival_df$time))
rev_survival_df[rev_survival_df$status=="Alive",]$status <- 0
rev_survival_df[rev_survival_df$status=="Dead",]$status <- 1
rev_survival_df$status <- as.numeric(as.vector(rev_survival_df$status))

for (i in 2:ncol(rev_survival_df)){
  rev_survival_df[,i] <- as.numeric(as.vector(rev_survival_df[,i]))
}

rev_survival_df$time <- rev_survival_df$time/365

#Normalize the dataset
mean_vector <- NULL
sd_vector <- NULL

for (i in 4:ncol(rev_survival_df)){
  mean_val <- mean(rev_survival_df[,i])
  sd_val <- sd(rev_survival_df[,i])
  rev_survival_df[,i] <- (rev_survival_df[,i]-mean_val)/sd_val
  mean_vector <- c(mean_vector, mean_val)
  sd_vector <- c(sd_vector,sd_val)
}

#Removing patient id, status and time
colnames_of_interest <- colnames(rev_survival_df)[-c(1,2,3)]
y_cols <- c("time","status")

dir.create("./Analysis/Microbiome/glmnet_model_microbiome", showWarnings = FALSE)
save(rev_survival_df, colnames_of_interest, y_cols, revised_Genus_full_abundance, mean_vector, sd_vector,
     file = "./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_AC_ICAM_246_for_glmnet.Rdata")

