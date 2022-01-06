
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

Files = list.files("./Analysis/TCR/7_Overlapping_sequences_analysis/7.2_Per_patient_TN")
Patients = substring(Files, 1, 3)

results = data.frame(Patient = Patients, Total = NA, Number_Tumor_Enriched = NA, Number_Normal_Enriched = NA, Number_not_Enriched = NA)

i=3
for (i in 1:10){
  Patient = Patients[i]
  load(paste0("./Analysis/TCR/7_Overlapping_sequences_analysis/7.2_Per_patient_TN/", Patient, "_TN_clones.Rdata"))
  subset$T_N_ratio = subset[, paste0(Patient, "T")]/subset[, paste0(Patient, "N")]
  subset$cat_beausang = "Not enriched"
  subset$cat_beausang[which(subset$T_N_ratio > 32 & subset$cat == "Overlapping" & subset[,paste0(Patient, "T")] > 0.10)] = "Tumor enriched"
  subset$cat_beausang[which(subset[,paste0(Patient, "T")] > 0.10 & subset$cat == "Tumor restricted")] = "Tumor enriched"
  subset$cat_beausang[which(subset$T_N_ratio < 1/32 & subset$cat == "Overlapping" & subset[,paste0(Patient, "N")] > 0.10)] = "Normal enriched"
  subset$cat_beausang[which(subset[,paste0(Patient, "N")] > 0.10 & subset$cat == "Normal restricted")] = "Normal enriched"
  subset$cat_beausang = factor(subset$cat_beausang, levels = c("Tumor enriched", "Not enriched", "Normal enriched"))
  subset$abundance = "non abundant sequence"
  subset$abundance[which(subset[,paste0(Patient, "N")] > 0.10 | subset[,paste0(Patient, "T")] > 0.10)] = "abundant sequence"
  dir.create("./Analysis/TCR/7_Overlapping_sequences_analysis/8_Per_patient_TN/", showWarnings = FALSE)
  save(subset, file = paste0("./Analysis/TCR/7_Overlapping_sequences_analysis/8_Per_patient_TN/", Patient, "_TN_clones_as_beausang.Rdata"))
  subset_v2 = subset[which(subset$abundance == "abundant sequence"),]
  tbl = table(subset_v2$cat_beausang)
  results$Total[which(results$Patient == Patient)] = nrow(subset_v2)
  results$Number_Tumor_Enriched[which(results$Patient == Patient)] = tbl["Tumor enriched"]
  results$Number_Normal_Enriched[which(results$Patient == Patient)] = tbl["Normal enriched"]
  results$Number_not_Enriched[which(results$Patient == Patient)] = tbl["Not enriched"]
}

results$Fraction_Tumor_Enriched = results$Number_Tumor_Enriched/results$Total
  
dir.create("./Analysis/TCR/8_As_Beausang_TN_Enriched/", showWarnings = FALSE)
save(results, file = "./Analysis/TCR/8_As_Beausang_TN_Enriched/8_Tumor_enriched_sequences.Rdata")
