
### T cell metrics

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Load data
MiXCR_orig = read.csv("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Diversity.csv", 
                 stringsAsFactors = FALSE)
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")

# Prepare data
MiXCR = MiXCR_orig
colnames(MiXCR)[1] = "Sample_ID"
MiXCR$Patient_ID = substring(MiXCR$Sample_ID, 1, 3)
MiXCR$Tissue = substring(MiXCR$Sample_ID, 4, 6)

# Filter only tumor samples
MiXCR = MiXCR[which(MiXCR$Tissue == "T-P"),]
MiXCR = MiXCR[which(substring(MiXCR$Sample_ID, 7, 10) == "_TRB"),]

MiXCR = MiXCR[which(MiXCR$Patient_ID %in% substring(rownames(immune_sig_df), 1, 3)),]

missing = rownames(immune_sig_df)[-which(substring(rownames(immune_sig_df), 1, 3) %in% MiXCR$Patient_ID)]

MiXCR = MiXCR[-which(is.na(MiXCR$Pielou)),]
MiXCR$TCR_Clonality = 1 - MiXCR$Pielou

MiXCR = MiXCR[-which(MiXCR$TCR_Clonality == -Inf),]

ES = t(immune_sig_df)
ES = as.data.frame(ES)

ES = as.matrix(ES)

variables = rownames(ES)
N.variables = length(variables)

results_clonality = data.frame(Variable = variables, cor_coef = NA, p_val = NA)

i = 1
for (i in 1:N.variables){
  variable = variables[i]
  MiXCR[,variable] = ES[variable, ][match(substring(MiXCR$Sample_ID, 1, 4), substring(colnames(ES), 1, 4))]
  MiXCR[,variable] = as.numeric(MiXCR[,variable])
  
  cor_test = cor.test(x = MiXCR[,variable], y = MiXCR$TCR_Clonality, method = "pearson")
  correlation = cor_test$estimate
  p_val = cor_test$p.value
  results_clonality$cor_coef[which(results_clonality$Variable == variable)] = correlation
  results_clonality$p_val[which(results_clonality$Variable == variable)] = p_val
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_MiXCR_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  #png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_v3_Scatterplot_", variable, "_TCR_cell_clonality.png"), res = 600,
  #   width = 3.2, height = 3, units = "in")
  #plot(plot4)
  #dev.off()
}

results_clonality = results_clonality[order(results_clonality$cor_coef, decreasing = TRUE),]
results_clonality$cor_coef = round(results_clonality$cor_coef, 3)
results_clonality$p_val = signif(results_clonality$p_val, 3)
results_clonality$FDR = p.adjust(results_clonality$p_val, method = "BH")

dir.create("./Analysis/TCR/006_MiXCR_T_cell_infiltration_cor", showWarnings = FALSE)
write.csv(results_clonality, file = "./Analysis/TCR/006_MiXCR_T_cell_infiltration_cor/April_2021_006_MiXCR_1-Pierlou_Evenness.csv", row.names = FALSE)