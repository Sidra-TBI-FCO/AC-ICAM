
### QC report of TCR Sequencing

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Load data
kit1_QC = read.table("./TCR_Experiments/Analysis/Experiments/QC/qcReport.2020-08-04_08-02-21.tsv", sep = "\t", header = TRUE, row.names = NULL)
kit2_QC = read.table("./TCR_Experiments/Analysis/Experiments/QC/qcReport.2020-08-04_07-58-50.tsv", sep = "\t", header = TRUE, row.names = NULL)
kit3_QC = read.table("./TCR_Experiments/Analysis/Experiments/QC/qcReport.2020-08-04_08-05-08.tsv", sep = "\t", header = TRUE, row.names = NULL)

kit_data = read.csv("./Processed_Data/TCR/TCR SampleOverview/SampleOverview_kit_1_2_3.csv", stringsAsFactors = FALSE)
kit_data$X = NULL

# Clean up files
colnames(kit1_QC) = colnames(kit1_QC)[2:10]
kit1_QC[,10] = NULL

colnames(kit2_QC) = colnames(kit2_QC)[2:10]
kit2_QC[,10] = NULL

colnames(kit3_QC) = colnames(kit3_QC)[2:10]
kit3_QC[,10] = NULL

all = rbind(kit1_QC, kit2_QC, kit3_QC)

kit_data$Coverage = all$Coverage[match(kit_data$sample_name, all$Sample.Name)]
kit_data$Adaptive_kit_run_ID = all$Run.ID[match(kit_data$sample_name, all$Sample.Name)]
kit_data$sample_tags = NULL
kit_data$sku = NULL
kit_data$test_name = NULL

kit_data = kit_data[, c(1, 11, 10, 2:9)]

dir.create("./Analysis/TCR/QC_001", showWarnings = FALSE)
write.csv(kit_data, file = "./Analysis/TCR/QC_001/overview_results_TCR.csv", row.names = FALSE)

