
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full"    # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"
subset = "" # "nonhypermutated" or "hypermutated"

# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))
rank_abundancies = get(paste0(Rank, "_abundance"))

frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

dir.create("./Figures/Microbiome/008.2_ICR_wilcoxon_boxplot", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/008.2_ICR_wilcoxon_boxplot/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/008.2_ICR_wilcoxon_boxplot/", Type, "/", Rank), showWarnings = FALSE)

# Analysis
abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)

abundance_T = as.data.frame(abundance_T)
abundance_T$ICR_cluster = table_cluster_assignment$ICR_HML[match(rownames(abundance_T),
                                                                 substring(rownames(table_cluster_assignment), 1, 3))]
abundance_T$ICRscore = table_cluster_assignment$ICRscore[match(rownames(abundance_T),
                                                               substring(rownames(table_cluster_assignment), 1, 3))]


df_plot = data.frame(Patient_ID = rownames(abundance_T),
                     ICR_cluster = abundance_T$ICR_cluster,
                     ICRscore = abundance_T$ICRscore,
                     Mutation_cat = NA,
                     Abundance = NA)
df_plot$ICR_cluster = factor(df_plot$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
df_plot$Mutation_cat = frequency_df$Mutation_cat[match(df_plot$Patient_ID, frequency_df$Patient_ID)]

if(subset %in% c("nonhypermutated", "hypermutated")){
  df_plot = df_plot[which(df_plot$Mutation_cat == subset),]
}


results = data.frame(Name = rownames(rank_abundancies), p_val = NA,
                     mean_ICR_Low = NA,
                     mean_ICR_Medium = NA, mean_ICR_High = NA)
i= 1
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df_plot$Abundance = abundance_T[, micr][match(df_plot$Patient_ID, rownames(abundance_T))]
  
  df_calc = df_plot
  test = wilcox.test(df_calc$Abundance[which(df_plot$ICR_cluster == "ICR High")], 
                     df_calc$Abundance[which(df_plot$ICR_cluster == "ICR Low")])
  
  results$p_val[which(results$Name == micr)] = test$p.value
  results$mean_ICR_Low[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$ICR_cluster == "ICR Low")])
  results$mean_ICR_Medium[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$ICR_cluster == "ICR Medium")])
  results$mean_ICR_High[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$ICR_cluster == "ICR High")])
  
  if(Type == "Relative"){
    df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
  }
  
  if(is.na(test$p.value)){next}
  if(test$p.value < 0.053 | make_plot == "make_plot"){
    plot = ggplot(df_plot, aes(x = ICR_cluster, y = Abundance)) +
      stat_summary(fun.y = "median", colour = "grey", size = 1, geom = "point") +
      geom_violin(aes(fill=ICR_cluster)) +
      geom_boxplot(width=0.1, outlier.shape = NA) +
      scale_fill_manual(values = c("ICR Low" = "blue", "ICR Medium" = "green", "ICR High" = "red")) +
      #geom_jitter(width = 0.1, size = 0.8, aes(color = Mutation_cat)) +
      #geom_point(size = 0.3) +
      #scale_color_manual(values = c("nonhypermutated" = "darkgreen", "hypermutated" = "orchid")) +
      theme_bw() +
      #stat_summary(fun.y = "median", colour = "grey", size = 2, geom = "point") +
      scale_y_log10() +
      ylab(paste0(Type, " abundance", "\n", gsub(".*\\D5__", "", micr))) +
      xlab("") +
      #facet_grid(.~Mutation_cat) +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black",
                                       angle = 90, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"),
            strip.background = element_blank(),
            legend.position = "none",
            strip.text = element_text(size = 15, colour = "black"))
    
    png(paste0("./Figures/Microbiome/008.2_ICR_wilcoxon_boxplot/", Type, "/", Rank, "/v3_", micr, "_", subset,
               "_by_ICR_cluster_boxplot.png"), res = 600, width = 3, height = 4, units = "in")
    plot(plot)
    dev.off()
  }
  
}

results = results[order(results$p_val),]
results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

dir.create(paste0("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot"), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot/",  Type, "/", Rank,
                            "/", subset, "_by_ICR_cluster_", Type, "_", Rank, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot/",  Type, "/", Rank,
                                 "/", subset ,"_by_ICR_cluster_", Type, "_", Rank, ".csv"),
          row.names = FALSE)



