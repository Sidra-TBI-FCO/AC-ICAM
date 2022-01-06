
# Overview Clinical Characteristics of Patient subgroup for TCR Sequencing

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Set parameters


# Load data
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
clinical_data = clinical_data[which(clinical_data$Patient_ID %in% TCR_Overview$Patient_ID),]

# Median follow-up
median(as.numeric(as.character(clinical_data$last_contact_days_to)))/365

table(clinical_data$gender)

# Age
clinical_data$age_at_initial_pathologic_diagnosis <- as.numeric(as.character(clinical_data$age_at_initial_pathologic_diagnosis))
mean(clinical_data$age_at_initial_pathologic_diagnosis)

# Categories age in <65 years, 65-74 years, >=75 years.
age_category1 <- sum(clinical_data$age_at_initial_pathologic_diagnosis<50)
age_category2 <- sum(clinical_data$age_at_initial_pathologic_diagnosis>= 50 & clinical_data$age_at_initial_pathologic_diagnosis < 65)
age_category3 <- sum(clinical_data$age_at_initial_pathologic_diagnosis>= 65 & clinical_data$age_at_initial_pathologic_diagnosis < 75)
age_category4 <- sum(clinical_data$age_at_initial_pathologic_diagnosis>= 75)
age_category1
age_category2
age_category3
age_category4
age_category1 + age_category2 + age_category3 + age_category4

# TNM Stage
table(clinical_data$ajcc_tumor_pathologic_pt)
table(clinical_data$ajcc_nodes_pathologic_pn)
table(clinical_data$ajcc_metastasis_pathologic_pm)
table(clinical_data$ajcc_pathologic_tumor_stage)

# Tumour anatomic location
table(clinical_data$tumour_anatomic_site)


# Tumour morphology
table(clinical_data$Tumor_morphology)

# Adjuvant treatment
table(clinical_data$Adjuvant_treatment)

# Year of diagnosis
table(clinical_data$year_of_initial_diagnosis)

# Recurrences
table(clinical_data$locoregional_recurrence_status)
table(clinical_data$distant_recurrence_status)

# History of cancer
table(clinical_data$history_other_malignancy)

