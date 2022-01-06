# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "openxlsx")
ipak(required.packages)

patient.samples = read.csv("./Processed_Data/Shared_Data/For cbioportal/New_to_be_uploaded/check_EA/SER-SILU-Sample_to_patient_association_v3.csv")
clinical = read.csv("./Processed_Data/Shared_Data/For cbioportal/New_to_be_uploaded/SER-SILU-CC-clinical_data_for_cbioportal_2405_2021.csv")

patient.samples = patient.samples[which(patient.samples$Patient_ID %in% clinical$Patient_ID),]
patient.samples = patient.samples[,c(1,2,6)]
patient.samples = patient.samples[-grep("AN", patient.samples$Sample),]
patient.samples = patient.samples[-grep("LM", patient.samples$Sample),]
patient.samples = patient.samples[-grep("-PT-02", patient.samples$Sample),]
patient.samples = patient.samples[-grep("-PT-03", patient.samples$Sample),]

patient.samples$Patient_ID %in% clinical$Patient_ID

write.csv(patient.samples, file = "./Processed_Data/Shared_Data/For cbioportal/New_to_be_uploaded/Patient_sample_association_files.csv")
