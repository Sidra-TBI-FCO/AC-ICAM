
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/heatmap.3.R"))

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper", "heatmap3", "plyr", "spatstat",
                      "heatmap3", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
gene_of_interest = "AGR2" # MUC2 AGR2
morphology_subset = "" # mucineus adenocarcinoom (dutch version)
hypermut_subset = "" # "hypermutated" or "nonhypermutated"

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")

frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"
frequency_df$Patient_ID = as.character(frequency_df$Patient_ID)

if(morphology_subset == ""){}else{
  included_patients = clinical_data$Patient_ID[which(clinical_data$Tumor_morphology == morphology_subset)]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% included_patients)]
}

if(hypermut_subset == ""){}else{
  included_patients = frequency_df$Patient_ID[which(frequency_df$Mutation_cat == hypermut_subset)]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% included_patients)]
}

plot_df = data.frame(Sample = colnames(RNASeq.QN.LOG2), Expression = RNASeq.QN.LOG2[gene_of_interest,], ICR = NA)

plot_df$ICR = table_cluster_assignment$ICRscore[match(plot_df$Sample, rownames(table_cluster_assignment))]

dir.create("./Figures/Trimmed_p/Gene_of_interest_plots", showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/02_Gene_by_ICR_cluster_boxplot", showWarnings = FALSE)
plot = ggplot(plot_df, aes(x = ICR, y = Expression)) +
  geom_point(size = 0.4) +
  stat_cor(method = "pearson", size = 5) +
  geom_smooth(method="lm") +
  xlab("ICR score") +
  ylab(gene_of_interest) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 17, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        axis.text.x = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 17, colour = "black"))

dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/11_Gene_by_ICRscore_scatter", showWarnings = FALSE)
png(filename = paste0("./Figures/Trimmed_p/Gene_of_interest_plots/11_Gene_by_ICRscore_scatter/", gene_of_interest, "_ICR_scatterplot",
                      morphology_subset, hypermut_subset, ".png"), 
    height = 4, width = 4, units = "in", res = 600)
plot(plot)
dev.off()
