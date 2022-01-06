
# Boxplot of somatic mutation in gene(s) of interest versus ICR score

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))


required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
selection = "cancer_driver_QGP" # "cancer_driver_501" or "cancer_driver_QGP" or ""
genes = c("B2M")
subset = "all_patients"  # "all_patients" or "nonhypermutated" or "hypermutated"
#genes = c("CXCL10", "CXCL6", "CXCL14", "CXCL13",
          #"CXCL2", "CXCL12", "CXCL11", "CXCL16", "CXCL5", "CXCL17", "CXCL3")
#genes = c("TTN", "APC", "TP53", "SYNE1", "MUC16", "KRAS", "FAT4", "MUC5B", "OBSCN", 
#          "DNAH5", "LRP1B", "KMT2C", "RYR2", "PCLO", "BRAF", "FAT2", "CSMD1", "RYR1", "USH2A", 
#          "FAT3", "FLG", "PDPR", "CSMD3", "PLEC", "SPTA1", "ZFHX4", "CROCC", "DNAH11", "MUC4", "RYR3")

dir.create("./Figures/WES/003_Boxplot_ICRscore_MUT_vs_WT", showWarnings = FALSE)
dir.create("./Analysis/WES/003_ICRscore_MUT_vs_WT", showWarnings = FALSE)
           
# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(paste0("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata"))
load("./Processed_Data/External/Gene_collections/501_genes_Michele.Rdata")
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")  # QGPC genes are 1226 genes in total
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")

finalMafFiltered$Hugo_Symbol = as.character(finalMafFiltered$Hugo_Symbol)
all_patients = unique(finalMafFiltered$Patient_ID)
frequency_df$Patient_ID = as.character(frequency_df$Patient_ID)
nonhypermutated = frequency_df$Patient_ID[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)]
hypermutated = frequency_df$Patient_ID[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)]
table_cluster_assignment$Patient_ID = substring(rownames(table_cluster_assignment), 1, 3)

finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Patient_ID %in% get(subset)),]
ICR_sub = table_cluster_assignment[which(table_cluster_assignment$Patient_ID %in% get(subset)),]

all_genes = as.character(unique(finalMafFiltered$Hugo_Symbol))
results = data.frame(Gene = all_genes, p_value = NA, CI_lower = NA, CI_upper = NA, mean_WT = NA, mean_MUT = NA, n_WT = NA, n_MUT = NA)

j= all_genes[1]
for (j in all_genes){
  gene = j
  sub = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  mutations_freq_df = data.frame(table(sub$Patient_ID))
  if(nrow(mutations_freq_df) == 1){next}
  colnames(mutations_freq_df) = c("Patient_ID", "Number_of_mutations")
  plot_df = data.frame(Patient_ID = ICR_sub$Patient_ID, ICRscore = ICR_sub$ICRscore,
                       Mutation_status = NA)
  plot_df$Mutation_status[which(plot_df$Patient_ID %in% mutations_freq_df$Patient_ID)] = "MUT"
  plot_df$Mutation_status[-which(plot_df$Patient_ID %in% mutations_freq_df$Patient_ID)] = "WT"
  plot_df$Mutation_status = factor(plot_df$Mutation_status, levels = c("WT", "MUT"))
  t_test = t.test(plot_df$ICRscore[which(plot_df$Mutation_status == "WT")], plot_df$ICRscore[which(plot_df$Mutation_status == "MUT")])
  results$p_value[which(results$Gene == gene)] = t_test$p.value
  results$CI_lower[which(results$Gene == gene)] = t_test$conf.int[1]
  results$CI_upper[which(results$Gene == gene)] = t_test$conf.int[2]
  results$mean_WT[which(results$Gene == gene)] = t_test$estimate[1]
  results$mean_MUT[which(results$Gene == gene)] = t_test$estimate[2]
  results$n_WT[which(results$Gene == gene)] = table(plot_df$Mutation_status)["WT"]
  results$n_MUT[which(results$Gene == gene)] = table(plot_df$Mutation_status)["MUT"]
}

results$delta = results$mean_MUT - results$mean_WT
results = results[order(results$p_value),]
results = results[-which(is.na(results$mean_WT)),]
results_all = results

if(selection == "cancer_driver_501"){
  results = results[which(results$Gene %in% genes_Michele),]
}
if(selection == "cancer_driver_QGP"){
  results = results[which(results$Gene %in% QGPC_genes),]
}
if(selection == "All"){
  results = results
}

results$FDR = p.adjust(results$p_value, method = "fdr", n = nrow(results))
write.csv(results, file = paste0("./Analysis/WES/003_ICRscore_MUT_vs_WT/", subset, "_", selection ,"_table_MUT_WT_ICRscore",
                                 ".csv"),
          row.names = FALSE)

results = results[which(results$n_MUT >4),]
results$FDR = p.adjust(results$p_value, method = "fdr", n = nrow(results))
results_sub_FDR = results[which(results$FDR < 0.05),]
results_sub_sig = results[which(results$p_value < 0.05),]
results_sub_sig = results_sub_sig[order(results_sub_sig$delta),]

write.csv(results_sub_sig, file = paste0("./Analysis/WES/003_ICRscore_MUT_vs_WT/", subset, "_", selection ,"_Table_significant_difference_MUT_WT_ICRscore_min_5_MUT.csv"),
          row.names = FALSE)
genes = as.character(results_sub_sig$Gene)
i=1
for (i in genes){
  gene = i
  finalMafFiltered_sub = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  mutations_freq_df = data.frame(table(finalMafFiltered_sub$Patient_ID))
  colnames(mutations_freq_df) = c("Patient_ID", "Number_of_mutations")
  plot_df = data.frame(Patient_ID = ICR_sub$Patient_ID, ICRscore = ICR_sub$ICRscore,
                       Mutation_status = NA)
  plot_df$Mutation_status[which(plot_df$Patient_ID %in% mutations_freq_df$Patient_ID)] = "MUT"
  plot_df$Mutation_status[-which(plot_df$Patient_ID %in% mutations_freq_df$Patient_ID)] = "WT"
  plot_df$Mutation_status = factor(plot_df$Mutation_status, levels = c("WT", "MUT"))
  
  plot = ggplot(plot_df, aes(x = Mutation_status, y = ICRscore, fill = Mutation_status)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("WT" = "#0064FF", "MUT" = "#FF0166")) +
    geom_jitter(width = 0.2, size = 0.8) +
    theme_bw() +
    ggtitle(gene) +
    xlab(paste0("Mutation status")) +
    ylab("ICR score") +
    ylim(3.5, 11) +
    stat_compare_means(method = "t.test", comparisons = list(c("MUT", "WT"))) +
    theme(axis.title.x = element_text(size = 15, color = "black"),
          axis.title.y = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 15, color = "black"),
          title = element_text(size = 15, color = "black"),
          legend.position = "none")
  
  dir.create(paste0("./Figures/WES/003_Boxplot_ICRscore_MUT_vs_WT/", subset), showWarnings = FALSE)
  png(filename = paste0("./Figures/WES/003_Boxplot_ICRscore_MUT_vs_WT/", subset, "/Boxplot_Mut_status_", gene, "_vs_ICRscore.png"),
      res = 600, width = 2.8, height = 4, units = "in")
  plot(plot)
  dev.off()
}


