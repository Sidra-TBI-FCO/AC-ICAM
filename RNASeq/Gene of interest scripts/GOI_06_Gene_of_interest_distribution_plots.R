
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ggplot2", "easyGgplot2")
ipak(required.packages)

# Set parameters
load("./Processed_Data/External/Gene_collections/DUX4_56_Gene_Symbol.Rdata")
load("./Processed_Data/External/Gene_collections/DUX4_Yao_genes.Rdata")
genes_of_interest =  Yao_DUX4_genes

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
genes_of_interest = genes_of_interest[which(genes_of_interest %in% rownames(RNASeq.QN.LOG2))]

i=1
for (i in 1:length(genes_of_interest)){
  gene_of_interest = genes_of_interest[i]
  if(gene_of_interest %in% rownames(RNASeq.QN.LOG2)){
    plot_df = data.frame(Sample = colnames(RNASeq.QN.LOG2), Expression = RNASeq.QN.LOG2[gene_of_interest,])
    
    plot = ggplot(plot_df, aes(Expression)) +
      geom_histogram(color = "black", fill = "lightblue") +
      theme_bw() +
      xlab(paste0(gene_of_interest, " expression value")) +
      ylab("number of patients")
    
    dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/06_Distribution_plots", showWarnings = FALSE)
    
    assign(paste0("p", i), plot)
    
    #png(paste0("./Figures/Trimmed_p/Gene_of_interest_plots/06_Distribution_plots/", gene_of_interest, "_histogram.png"),
     #   width = 3, height = 2, units = "in", res = 600)
    #plot(plot)
    #dev.off()
  }else{
    print(paste0(gene_of_interest, " not available"))
  }
}

plots = paste("p", 1:length(genes_of_interest), sep = "")
list_of_plots = mget(plots)

dir.create(paste0("./Figures/Trimmed_p/Gene_of_interest_plots/06_Distribution_plots/Multiplot"), showWarnings = FALSE)
png(paste0("./Figures/Trimmed_p/Gene_of_interest_plots/06_Distribution_plots/Multiplot/Genes.png"),
    width = 30, height = 10, units = "in", res = 600)

#width = 8.27, height = 10, units = "in", res = 600)
#png(paste0("./Figures/Multi_panel_boxplots/ggpaired_core_network_analysis_1_", signif.cutoff, "_multiplot_p", 
#number_start,"-p", number_end,".png"), width = 11, height = 11, units = "in", res = 600)
ggplot2.multiplot(plotlist = list_of_plots[1:length(genes_of_interest)], cols=13)
dev.off()
