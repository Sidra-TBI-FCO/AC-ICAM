
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("reshape2", "ggplot2", "plotrix", "ggpubr"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"


# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))

rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)
abundance_N = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "N")]
colnames(abundance_N) = substring(colnames(abundance_N), 1, 3)
abundance_N = t(abundance_N)

rownames(abundance_N) == rownames(abundance_T) # check if all is TRUE
colnames(abundance_N) == colnames(abundance_T) # check if all is TRUE
abundance = abundance_T - abundance_N

df = data.frame(Patient_ID = unique(substring(colnames(rank_abundancies), 1, 3)),
                Abundance = NA,
                ICR_cluster = NA,
                stringsAsFactors = FALSE)

df$ICR_cluster = table_cluster_assignment$ICR_HML[match(df$Patient_ID,substring(rownames(table_cluster_assignment), 1, 3))]
df$ICR_cluster = factor(df$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

dir.create("./Figures/Microbiome/007.3_Delta_TN_by_ICR_cluster", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/007.3_Delta_TN_by_ICR_cluster/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/007.3_Delta_TN_by_ICR_cluster/", Type, "/", Rank), showWarnings = FALSE)

results = data.frame(Name = rownames(rank_abundancies), p_val = NA, Spearman_cor = NA,
                     mean_ICR_Low = NA, mean_ICR_Medium = NA, mean_ICR_High = NA)

i = 14
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[,micr][match(df$Patient_ID, rownames(abundance))]
  
  df_calc = df
  levels(df_calc$ICR_cluster) = c("1", "2", "3")
  df_calc$ICR_cluster = as.numeric(df_calc$ICR_cluster)
  
  cor = cor.test(df_calc$Abundance, df_calc$ICR_cluster, method = "spearman")
  results$p_val[which(results$Name == micr)] = cor$p.value
  results$Spearman_cor[which(results$Name == micr)] = cor$estimate
  results$mean_ICR_Low[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$ICR_cluster == 1)])
  results$mean_ICR_Medium[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$ICR_cluster == 2)])
  results$mean_ICR_High[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$ICR_cluster == 3)])
  
  ylabel = gsub(".*\\D_5__", "", micr)
  
  if(ylabel %in% c("Fusobacterium", "Campylobacter", "Parabacteroides", "Alistipes") | make_plot == "make_plot"){
    plot = ggplot(df, aes(x = ICR_cluster, y = Abundance)) +
      geom_violin(trim = TRUE, aes(fill=ICR_cluster)) +
      stat_compare_means(method = "wilcox", comparisons = list(c("ICR High", "ICR Medium"),
                                                               c("ICR Medium", "ICR Low"),
                                                               c("ICR High", "ICR Low"))) +
      #geom_jitter(position = "dodge",size=0.5) +
      scale_fill_manual(values = c("ICR Low" = "blue", "ICR Medium" = "green", "ICR High" = "red")) +
      geom_boxplot(width=0.1, outlier.shape = NA) +
      ylab(paste0(ylabel, "\n", "Delta tumor-normal")) +
      theme_bw() +
      #stat_summary(fun.y = "median", colour = "grey", size = 1, geom = "point") +
      #scale_y_log10() +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"),
            strip.background = element_blank(),
            strip.text = element_text(size = 15, colour = "black"),
            legend.position = "none")
    
    png(paste0("./Figures/Microbiome/007.3_Delta_TN_by_ICR_cluster/", Type, "/", Rank, "/",
               micr, "_boxplot_Delta_T_N.png"), res = 600, units = "in", width = 2.8, height = 3.5)
    plot(plot)
    dev.off()
  }
}

results = results[order(results$p_val),]
results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

dir.create("./Analysis/Microbiome/007.3_Delta_TN_by_ICR_cluster", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/007.3_Delta_TN_by_ICR_cluster/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/007.3_Delta_TN_by_ICR_cluster/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/007.3_Delta_TN_by_ICR_cluster/",  Type, "/", Rank,
                            "/","Spearman_by_ICR_cluster_", Type, "_", Rank, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/007.3_Delta_TN_by_ICR_cluster/",  Type, "/", Rank,
                                 "/", "Spearman_by_ICR_cluster_", Type, "_", Rank, ".csv"),
          row.names = FALSE)
