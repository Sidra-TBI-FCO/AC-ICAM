
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
  plot = ggplot(plot_df, aes(x=ICR_cluster, y = ES, fill = ICR_cluster)) +
    facet_grid(.~CMS) +
    geom_rect(aes(fill = CMS), xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.3) +
    scale_fill_manual(values = c("ICR High" = alpha("red", 1), 
                                 "ICR Medium" = alpha("green", 1),
                                 "ICR Low" = alpha("blue", 1),
                                 "CMS1" = alpha("#FFD09E", 1),
                                 "CMS2" = alpha("#97BDD8", 1), 
                                 "CMS3" = alpha("#F4C0D5", 1), 
                                 "CMS4" = alpha("#99CFBE", 1))) +
    geom_violin() +
    geom_boxplot(width=.1, outlier.shape = NA) +
    #geom_jitter(width = 0.15, size = 0.2) +
    #theme_bw() +
    #stat_compare_means(method = "t.test", label = "p.signif",
    #comparisons = list(c("ICR High", "ICR Low"),
    #c("ICR Medium", "ICR Low"),
    #c("ICR High", "ICR Medium"))) +
    ylab(paste0(var)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 17),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 10),
          strip.background = element_blank(),
          strip.text = element_blank(),
          aspect.ratio = 1.3/1) 
  
  png(paste0("./Figures/Trimmed_p/024_violin_ESTIMATE/", var, "_violinplot.png"),
      width = 6, height = 2, units = "in", res = 600)
  plot(plot)
  dev.off()
  
  assign(paste0("p", i), plot)
}

plots = paste("p", 1:ncol(ESTIMATE), sep = "")
list_of_plots = mget(plots)

png(paste0("./Figures/Trimmed_p/024_violin_ESTIMATE/Multiplot.png"),
    width = 10, height = 2*length(To_plot), units = "in", res = 600)
ggplot2.multiplot(plotlist = list_of_plots[1:length(To_plot)], rows=length(To_plot),
                  cols = 1)
dev.off()


i=1
for (i in 1:nrow(ES)){
  Cell = rownames(ES)[i]
  plot_df$ES = ES[Cell,][match(plot_df$Patient_ID, substring(colnames(ES), 1, 3))]
  for (j in 1:5){
    CMS = unique(plot_df$CMS)[j]
    plot_df_CMS = plot_df[which(plot_df$CMS == CMS),]
    test1 = t.test(plot_df_CMS$ES[which(plot_df_CMS$ICR_cluster == "ICR High")],
                   plot_df_CMS$ES[which(plot_df_CMS$ICR_cluster == "ICR Low")])
    test2 = t.test(plot_df_CMS$ES[which(plot_df_CMS$ICR_cluster == "ICR Medium")],
                   plot_df_CMS$ES[which(plot_df_CMS$ICR_cluster == "ICR Low")])
    test3 = t.test(plot_df_CMS$ES[which(plot_df_CMS$ICR_cluster == "ICR High")],
                   plot_df_CMS$ES[which(plot_df_CMS$ICR_cluster == "ICR Medium")])
    
    results = data.frame(Comparison = c("ICR High- ICR Low", 
                                        "ICR Medium-ICR Low",
                                        "ICR High-ICR Medium"),
                         p_value = NA,
                         CI_lower = NA,
                         CI_upper = NA)
    results$p_value[which(results$Comparison == "ICR High- ICR Low")] = test1$p.value
    results$CI_lower[which(results$Comparison == "ICR High- ICR Low")] = test1$conf.int[1]
    results$CI_upper[which(results$Comparison == "ICR High- ICR Low")] = test1$conf.int[2]
    
    results$p_value[which(results$Comparison == "ICR Medium-ICR Low")] = test2$p.value
    results$CI_lower[which(results$Comparison == "ICR Medium-ICR Low")] = test2$conf.int[1]
    results$CI_upper[which(results$Comparison == "ICR Medium-ICR Low")] = test2$conf.int[2]
    
    results$p_value[which(results$Comparison == "ICR High-ICR Medium")] = test3$p.value
    results$CI_lower[which(results$Comparison == "ICR High-ICR Medium")] = test3$conf.int[1]
    results$CI_upper[which(results$Comparison == "ICR High-ICR Medium")] = test3$conf.int[2]
    
    assign(paste0(Cell, CMS), results)
    
  }
  
}


