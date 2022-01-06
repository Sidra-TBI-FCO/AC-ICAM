
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
mean_of = "" #c("HLA-A", "HLA-B", "HLA-C") or if no mean ""
morphology_subset = "non-mucinous" # mucineus adenocarcinoom (dutch version)
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

matrix_all = RNASeq.QN.LOG2

if(morphology_subset == ""){}

if(morphology_subset == "non-mucinous"){
  included_patients = clinical_data$Patient_ID[-which(clinical_data$Tumor_morphology == "mucineus adenocarcinoom")]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% included_patients)]
  }else{
  included_patients = clinical_data$Patient_ID[which(clinical_data$Tumor_morphology == morphology_subset)]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% included_patients)]
}

if(hypermut_subset == ""){}else{
  included_patients = frequency_df$Patient_ID[which(frequency_df$Mutation_cat == hypermut_subset)]
  #included_patients = included_patients[which(included_patients %in% substring(colnames(RNASeq.QN.LOG2), 1, 3))]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% included_patients)]
}


if(mean_of > 1){
  calc_matrix = RNASeq.QN.LOG2[mean_of,]
  calc_df = as.data.frame(t(calc_matrix))
  calc_df$Expression = rowMeans(calc_df)
  plot_df = data.frame(Sample = rownames(calc_df), Expression = calc_df$Expression, ICR = NA)
}else{
  plot_df = data.frame(Sample = colnames(RNASeq.QN.LOG2), Expression = RNASeq.QN.LOG2[gene_of_interest,], ICR = NA)
}

plot_df$ICR = table_cluster_assignment$ICR_HML[match(plot_df$Sample, rownames(table_cluster_assignment))]
plot_df$ICR = factor(plot_df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))


symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

dir.create("./Figures/Trimmed_p/Gene_of_interest_plots", showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/02_Gene_by_ICR_cluster_boxplot", showWarnings = FALSE)

png(filename = paste0("./Figures/Trimmed_p/Gene_of_interest_plots/02_Gene_by_ICR_cluster_boxplot/", gene_of_interest, "_ICR_box_plot",
                      morphology_subset, hypermut_subset, ".png"), 
    height = 3, width = 2, units = "in", res = 600)
plot = ggplot(plot_df, aes(x = ICR, y = Expression, fill = ICR)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("blue", "green", "red")) +
  ylim(c(min(matrix_all[gene_of_interest,]) * 0.8, max(matrix_all[gene_of_interest,]) *1.3)) +
  xlab("") +
  ylab(paste0(gene_of_interest, " expression")) +
  stat_compare_means(comparisons = list(c("ICR Low", "ICR Medium"),
                                        c("ICR Medium", "ICR High"),
                                        c("ICR Low", "ICR High")), method = "t.test", label = "p.signif") +
  #ggtitle(gene_of_interest) +
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 0.95, vjust = 0.5),
        legend.position = "none",
        axis.text.y = element_text(color = "black"))
plot(plot)
dev.off()

levels(plot_df$ICR) = c(1, 2, 3)
plot_df$ICR = as.numeric(plot_df$ICR)
cor_test = cor.test(plot_df$ICR, plot_df$Expression, method = "spearman")
cor_test$p.value
cor_test$estimate
