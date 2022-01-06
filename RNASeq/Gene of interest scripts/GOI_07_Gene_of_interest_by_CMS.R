
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
genes = c("IL12A")
#c("LRRC32","LTBP1","LTBP2","LTBP3","LTBP4", "MMP2", "MMP9",
 # "THBS1","ITGB6","ITGB3", "TGFB1", "TGFB2", "TGFB3")


mean_of = "" #c("HLA-A", "HLA-B", "HLA-C") or if no mean ""

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

for (i in 1:length(genes)){
  gene_of_interest = genes[i]
  
  if(mean_of > 1){
    calc_matrix = RNASeq.QN.LOG2[mean_of,]
    calc_df = as.data.frame(t(calc_matrix))
    calc_df$Expression = rowMeans(calc_df)
    plot_df = data.frame(Sample = rownames(calc_df), Expression = calc_df$Expression, CMS = NA)
  }else{
    plot_df = data.frame(Sample = colnames(RNASeq.QN.LOG2), Expression = RNASeq.QN.LOG2[gene_of_interest,], CMS = NA)
  }
  
  plot_df$CMS = Rfcms$RF.predictedCMS[match(plot_df$Sample, rownames(Rfcms))]
  plot_df$CMS = as.character(plot_df$CMS)
  plot_df$CMS[which(is.na(plot_df$CMS))] = "mixed"
  
  plot_df$CMS = factor(plot_df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "mixed"))
  
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
  
  
  dir.create("./Figures/Trimmed_p/Gene_of_interest_plots", showWarnings = FALSE)
  dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/07_Gene_by_CMS_boxplot", showWarnings = FALSE)
  
  png(filename = paste0("./Figures/Trimmed_p/Gene_of_interest_plots/07_Gene_by_CMS_boxplot/", gene_of_interest, "_CMS_box_plot.png"), 
      height = 4.5, width = 3.5, units = "in", res = 600)
  plot = ggplot(plot_df, aes(x = CMS, y = Expression, color = CMS)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.2) +
    theme_bw() +
    scale_color_manual(values = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                  "mixed" = "grey")) +
    #ylim(c(0, max(plot_df$Expression) *1.15)) +
    #ylim(8, 19) +
    xlab("") +
    ylab(paste0(gene_of_interest, " expression")) +
    stat_compare_means(comparisons = list(c("CMS3", "CMS4"), c("CMS2", "CMS4"), c("CMS1", "CMS4")), method = "t.test", label = "p.signif") +
    #ggtitle(gene_of_interest) +
    theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 0.95, vjust = 0.5, size = 15),
          legend.position = "none",
          axis.text.y = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 17))
  plot(plot)
  dev.off()
}





