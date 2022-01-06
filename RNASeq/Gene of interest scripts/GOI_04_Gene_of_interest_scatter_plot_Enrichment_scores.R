
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ggplot2", "ggpubr", "graphics")
ipak(required.packages)

# Set parameters
Gene_of_interest = "CDH1" #"ITGAE"
Gene_set = "ConsensusTME_COAD"

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/", Gene_set, "_ES.Rdata"))

# Analysis
plot_df = data.frame(t(ES))
plot_df$Gene_of_interest = RNASeq.QN.LOG2[Gene_of_interest,][match(rownames(plot_df), colnames(RNASeq.QN.LOG2))]

ES_vars = colnames(plot_df)[-which(colnames(plot_df) == "Gene_of_interest")]

dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/04_scatterplot_ES_vs_expression", showWarnings = FALSE)
i = 1
for (i in 1:length(ES_vars)){
  var = ES_vars[i]
  cor_coef = cor(plot_df[, var], plot_df$Gene_of_interest, method = "pearson")
  plot = ggplot(plot_df, aes(x = Gene_of_interest, y = get(var))) +
    geom_point(size = 0.2) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", size = 3) +
    theme_bw() +
    xlab(Gene_of_interest) +
    ylab(gsub("\\.", " ", var)) +
    #geom_text(aes(x = max(Gene_of_interest) * 0.88, y = max(plot_df[, var]), 
     #             label = paste0("cor = ", round(cor_coef, 2)))) +
    theme(text = element_text(colour = "black", size = 8),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 12))
  
  png(paste0("./Figures/Trimmed_p/Gene_of_interest_plots/04_scatterplot_ES_vs_expression/", Gene_of_interest, "_", var,
             "_scatterplot.png"), height = 3, width = 3, units = "in", res = 600)
  plot(plot)
  dev.off()
}
