
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

# Load data
load("./Analysis/TCR/8_As_Beausang_TN_Enriched/8.2_Tumor_enriched_sequences.Rdata")

Patients = unique(plot_df$Patient)

results = data.frame(Patient = Patients, number_of_unique_enriched_sequences = NA, productive_frequency_sum = NA)

i=1
for (i in 1:10){
  Patient = Patients[i]
  plot_df_sub = plot_df[which(plot_df$Patient == Patient),]
  results$number_of_unique_enriched_sequences[which(results$Patient == Patient)] = nrow(plot_df_sub)
  results$productive_frequency_sum[which(results$Patient == Patient)] = sum(plot_df_sub$Tumor)
}

save(results, file = "./Analysis/TCR/8_As_Beausang_TN_Enriched/8.3_Tumor_enriched_sequences_stats.Rdata")
