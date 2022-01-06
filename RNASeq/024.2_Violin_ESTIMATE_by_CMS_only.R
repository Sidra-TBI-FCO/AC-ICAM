
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "easyGgplot2")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
geneset = "ESTIMATE"
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Analysis
table_cluster_assignment$Patient_ID = substring(rownames(table_cluster_assignment), 1, 3)
table_cluster_assignment$CMS = Rfcms$RF.predictedCMS[match(table_cluster_assignment$Patient_ID, substring(rownames(Rfcms), 1, 3))]
table_cluster_assignment = table_cluster_assignment[-which(table_cluster_assignment$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]

plot_df = data.frame(Patient_ID = table_cluster_assignment$Patient_ID, ES = NA, ICR_cluster = table_cluster_assignment$ICR_HML, 
                     CMS = table_cluster_assignment$CMS)
plot_df$ICR_cluster = factor(plot_df$ICR_cluster, levels =c("ICR Low", "ICR Medium", "ICR High"))
plot_df$CMS = factor(plot_df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))

dir.create("./Figures/Trimmed_p/024_violin_ESTIMATE", showWarnings = FALSE)

plot_df = plot_df[-which(is.na(plot_df$CMS)),]

i=1
for (i in 1:ncol(ESTIMATE)){
  var = colnames(ESTIMATE)[i]
  plot_df$ES = ESTIMATE[,var][match(plot_df$Patient_ID, substring(rownames(ESTIMATE), 1, 3))]
  plot = ggplot(plot_df, aes(x=CMS, y = ES, fill = CMS)) +
    scale_fill_manual(values = c("CMS1" = alpha("#FFD09E", 1),
                                 "CMS2" = alpha("#97BDD8", 1), 
                                 "CMS3" = alpha("#F4C0D5", 1), 
                                 "CMS4" = alpha("#99CFBE", 1))) +
    theme_bw() +
    geom_violin() +
    geom_boxplot(width=.1, outlier.shape = NA) +
    ylab(paste0(var)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 17),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 10),
          strip.background = element_blank(),
          strip.text = element_blank(),
          aspect.ratio = 1.3/1) +
   # stat_compare_means(method = "t.test", comparisons = list(c("CMS1", "CMS2"),
    #                                                         c("CMS1", "CMS3"),
     #                                                        c("CMS1", "CMS4"),
      #                                                       c("CMS2", "CMS3"),
       #                                                      c("CMS2", "CMS4"),
        #                                                     c("CMS3", "CMS4")),
         #              label = "p.signif")
  
  png(paste0("./Figures/Trimmed_p/024_violin_ESTIMATE/024.2_", var, "_violinplot.png"),
      width = 6, height = 2, units = "in", res = 600)
  plot(plot)
  dev.off()
}