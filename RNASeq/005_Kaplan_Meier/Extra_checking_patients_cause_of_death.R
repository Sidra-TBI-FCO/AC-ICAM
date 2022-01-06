

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 10                                        # SET cut-off
Group.of.interest = "ICR_HML"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = ""
CMS = ""
ICR_cluster_with_metastasis = ""

# Read in the clinical data file
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")

# Add ICR as a variable and assign ICR cluster according to table cluster assignment

if(ICR_cluster_with_metastasis == "ICR_cluster_with_metastasis"){
  load(paste0("./Analysis/Trimmed_p/ICR Consensus Clustering/With_metastasis/JSREP_ICR_cluster_assignment_k2-6.Rdata"))
  table_cluster_assignment = table_cluster_assignment[-grep("LM", rownames(table_cluster_assignment)),]
  table_cluster_assignment = table_cluster_assignment[-grep("B", rownames(table_cluster_assignment)),]
  table_cluster_assignment = table_cluster_assignment[-grep("C", rownames(table_cluster_assignment)),]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(colnames(RNASeq.QN.LOG2) %in% rownames(table_cluster_assignment))]
}else{
  load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
}

load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("ICR High", "ICR Medium", "ICR Low"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

# Exclude adjuvant treated
if(exclude_adjuvant == "adjuvant_treated_excluded"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$Adjuvant_treatment == "No"),]
}

# Exclude Stage IV
if(exclude_stage_IV == "exclude_stage_IV" ){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$ajcc_pathologic_tumor_stage == "4"),]
}

if(CMS == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$CMS == CMS),]
}

Merged_dataset$FU_in_years = as.numeric(Merged_dataset$last_contact_days_to)/365
Merged_dataset$Final_age = round((as.numeric(Merged_dataset$age_at_initial_pathologic_diagnosis) + Merged_dataset$FU_in_years), 1)
Merged_dataset = Merged_dataset[order(Merged_dataset$ICR_HML),]

dir.create("./Analysis/Trimmed_p/005_Extra_checking_patients", showWarnings = FALSE)
write.csv(Merged_dataset, file = "./Analysis/Trimmed_p/005_Extra_checking_patients/005_Extra_checking_patients.csv",
          row.names = FALSE)
