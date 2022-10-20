
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
load("./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_AC_ICAM_246_for_glmnet.Rdata")

# Get the abundance information
tcga_bacteria_df <- fread("./Processed_Data/Microbiome/All_Input_Data_for_Risk_Model/bacteria.sample.relabund.genus.txt")
tcga_bacteria_df <- as.data.frame(tcga_bacteria_df)
bacteria_names <- tcga_bacteria_df$name

# Remove samples from biological replicates (duplicated patients)
duplicated_patients = substring(colnames(tcga_bacteria_df), 1, 12)[duplicated(substring(colnames(tcga_bacteria_df), 1, 12))]
duplicated_samples = colnames(tcga_bacteria_df)[which(substring(colnames(tcga_bacteria_df), 1, 12) %in% duplicated_patients)]

duplicated_samples = duplicated_samples[order(duplicated_samples)]
to_remove = duplicated_samples[duplicated(substring(duplicated_samples, 1, 12))]

tcga_bacteria_df = tcga_bacteria_df[,-which(colnames(tcga_bacteria_df) %in% to_remove)]

save(tcga_bacteria_df, file = "../NGS_Data_TCGA_COAD_Jessica/Processed_Data/Microbiome/Dohlman_TCGA_COAD_no_duplicate_patients_microbiome.Rdata")

tcga_bacteria_abundance_df <- tcga_bacteria_df[,c(2:ncol(tcga_bacteria_df))]

#Get TCGA clinical information
tcga_clinical_df <- fread("./Processed_Data/Microbiome/All_Input_Data_for_Risk_Model/TCGA_CLINICAL_DATA_CELL_2018_S1.csv")
tcga_clinical_df <- as.data.frame(tcga_clinical_df)
colnames_tcga_df <- str_replace(str_replace(str_replace(colnames(tcga_bacteria_abundance_df),"-01A",""),"-01B",""),"-01C","")
subset_tcga_clinical_df <- tcga_clinical_df[tcga_clinical_df$bcr_patient_barcode %in% colnames_tcga_df,]
colnames(tcga_bacteria_abundance_df) <- colnames_tcga_df
unique_colnames_tcga <- subset_tcga_clinical_df$bcr_patient_barcode

#Get the matching abundance information between clinical data and abundance info and if more than one patient exist then sum up the abundance levels
revised_tcga_bacteria_abundance_df <- NULL

i=9
for (i in 1:length(unique_colnames_tcga))
{
  name <- unique_colnames_tcga[i]
  ids <- which(colnames_tcga_df==name)
  if (length(ids)>1)
  {
    print("multiple samples per patients present")
    temp <- tcga_bacteria_abundance_df[,ids[1]]
    #temp <- rowSums(tcga_bacteria_abundance_df[,ids])                       #If more than one patient sample match then sum up abundance levels
    #temp <- rowMeans(tcga_bacteria_abundance_df[,ids])
  }else{
    print("sample is unique")
    temp <- tcga_bacteria_abundance_df[,ids]
  }
  revised_tcga_bacteria_abundance_df <- cbind(revised_tcga_bacteria_abundance_df, temp)
}
revised_tcga_bacteria_abundance_df <- as.data.frame(revised_tcga_bacteria_abundance_df)
colnames(revised_tcga_bacteria_abundance_df) <- unique_colnames_tcga
rownames(revised_tcga_bacteria_abundance_df) <- bacteria_names

#Make the mapping from the current bacteria names to tcga bacteria names and ids
all_bacterias_in_dataset <- colnames(rev_survival_df)[c(4:ncol(rev_survival_df))]
mapping_df <- NULL
for (i in 1:length(all_bacterias_in_dataset))
{
  Split <- strsplit(all_bacterias_in_dataset[i],split="__")
  covariate <- tail(Split[[1]], 1)
  mapped_bacteria <- NULL
  mapped_id <- NULL
  for (j in 1:length(bacteria_names))
  {
    ids <- grep(bacteria_names[j], covariate)
    if (length(ids)>0)
    {
      mapped_bacteria <- c(mapped_bacteria,bacteria_names[j])
      mapped_id <- c(mapped_id, j)
    }
  }
  if (is.null(mapped_bacteria)){
    temp <- cbind(all_bacterias_in_dataset[i],covariate,NA, NA)
  }else{
    temp <- cbind(all_bacterias_in_dataset[i],covariate, mapped_bacteria, mapped_id)
  }
  mapping_df <- rbind(mapping_df, temp)
}
mapping_df <- as.data.frame(mapping_df)
mapping_df$mapped_id <- as.numeric(as.vector(mapping_df$mapped_id))

#Build the TCGA Abundance matrix in the same format as our training dataset as then we can pass TCGA matrix as input to prediction model
tcga_abundance_df <- NULL
for (i in 1:nrow(mapping_df))
{
  name <- mapping_df$V1[i]
  if (name!="D_0__Bacteria D_1__Proteobacteria D_2__Gammaproteobacteria D_3__Enterobacteriales D_4__Enterobacteriaceae D_5__Escherichia-Shigella")
  {
    id <- mapping_df$mapped_id[i]
    if (!is.na(id))
    {
      temp <- revised_tcga_bacteria_abundance_df[id,]
    }else{
      temp <- rep(0, 117)
    }
    tcga_abundance_df <- rbind(tcga_abundance_df, temp)
  }
}
tcga_abundance_df <- rbind(tcga_abundance_df,colSums(revised_tcga_bacteria_abundance_df[mapping_df[mapping_df$V1=="D_0__Bacteria D_1__Proteobacteria D_2__Gammaproteobacteria D_3__Enterobacteriales D_4__Enterobacteriaceae D_5__Escherichia-Shigella",]$mapped_id,]))
tcga_abundance_df <- as.data.frame(tcga_abundance_df)
rownames(tcga_abundance_df) <- c(setdiff(mapping_df$V1,"D_0__Bacteria D_1__Proteobacteria D_2__Gammaproteobacteria D_3__Enterobacteriales D_4__Enterobacteriaceae D_5__Escherichia-Shigella"),"D_0__Bacteria D_1__Proteobacteria D_2__Gammaproteobacteria D_3__Enterobacteriales D_4__Enterobacteriaceae D_5__Escherichia-Shigella")
tcga_abundance_df <- as.data.frame(t(tcga_abundance_df))
N_col <- ncol(tcga_abundance_df)

# put Escherichia-Shigella in correct location in the matrix (6th before last position)
tcga_abundance_df <- tcga_abundance_df[,c(1:c(N_col-6),(N_col),(N_col-5),(N_col-4),(N_col-3),(N_col-2),(N_col-1))]

# JR QC-check:
test = data.frame(Microbiome = colnames(tcga_abundance_df), colMeans(tcga_abundance_df))
# confirmed with original dataframe subsetted to the same 117 samples !

#Build the validation set of TCGA survival and abundance
tcga_survival_df <- as.data.frame(cbind(subset_tcga_clinical_df[,c("bcr_patient_barcode","OS","OS.time")],tcga_abundance_df))
colnames(tcga_survival_df)[c(1:3)] <- c("samples","status","time")
tcga_survival_df$time <- as.numeric(as.vector(tcga_survival_df$time))
tcga_survival_df$status <- as.numeric(as.vector(tcga_survival_df$status))
tcga_survival_df$time <- tcga_survival_df$time/365

#Normalize the TCGA data
for (i in c(1:(ncol(tcga_survival_df)-3))){
  tcga_survival_df[,i+3] <- (tcga_survival_df[,i+3]-mean_vector[i])/sd_vector[i]
}

save(tcga_survival_df, colnames_of_interest, y_cols,
     file = "./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_TCGA_COAD_117_for_glmnet.Rdata")


