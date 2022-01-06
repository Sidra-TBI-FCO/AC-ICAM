
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "plotrix"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
Tissue = "T" # "T" or "N" or "Difference"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"


# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")

rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)
abundance_N = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "N")]
colnames(abundance_N) = substring(colnames(abundance_N), 1, 3)
abundance_N = t(abundance_N)

rownames(abundance_N) == rownames(abundance_T) # check if all is TRUE
colnames(abundance_N) == colnames(abundance_T) # check if all is TRUE
abundance_Difference = abundance_T - abundance_N
  
abundance = get(paste("abundance_", Tissue, sep = ""))

df = data.frame(Patient_ID = unique(substring(colnames(rank_abundancies), 1, 3)),
                Abundance = NA,
                Anatomical_location = NA,
                Primary_tumor_side = NA,
                stringsAsFactors = FALSE)

df$Anatomical_location = clinical_data$tumour_anatomic_site[match(df$Patient_ID, clinical_data$Patient_ID)]
df$Primary_tumor_side[which(df$Anatomical_location %in% c("ceceum", "colon ascendens", "flexura hepatica", 
                                                          "colon transversum"))] = "Right sided"
df$Primary_tumor_side[which(df$Anatomical_location %in% c("flexura lienalis", "colon descendens", "colon sigmoideum",
                                                          "rectosigmoideum"))] = "Left sided"
df$Primary_tumor_side = factor(df$Primary_tumor_side, levels = c("Right sided", "Left sided"))
table(df$Primary_tumor_side)

results = data.frame(Name = rownames(rank_abundancies), p_val = NA, wilcox_statistic = NA,
                     mean_right = NA, mean_left = NA)

dir.create("./Figures/Microbiome/014_Sideness_wilcoxon", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/014_Sideness_wilcoxon/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/014_Sideness_wilcoxon/", Type, "/", Rank), showWarnings = FALSE)

i = 12
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[,micr][match(df$Patient_ID, rownames(abundance))]
  
  test = wilcox.test(x = df$Abundance[which(df$Primary_tumor_side == "Right sided")], y = df$Abundance[which(df$Primary_tumor_side == "Left sided")], paired = FALSE,
                     alternative = "two.sided")
  results$p_val[which(results$Name == micr)] = test$p.value
  results$wilcox_statistic[which(results$Name == micr)] = test$statistic
  results$mean_right[which(results$Name == micr)] = mean(df$Abundance[which(df$Primary_tumor_side == "Right sided")])
  results$mean_left[which(results$Name == micr)] = mean(df$Abundance[which(df$Primary_tumor_side == "Left sided")])
  
  if(Type == "Relative" & Tissue %in% c("T", "N")){
    df$Abundance[which(df$Abundance == 0)] = 1e-4
  }
  
  melt = df
  
  if(is.na(test$p.value)){next}
  if(test$p.value < 0.05 | make_plot == "make_plot"){
    plot = ggplot(melt, aes(x = Primary_tumor_side, y = Abundance)) +
      geom_boxplot(outlier.shape = NA) +
      geom_line(aes(group = Patient_ID), color = alpha("lightgrey", 0.5)) +
      geom_point(size = 0.8) +
      ylab(paste0(micr, "\n", Type, " abudance")) +
      theme_bw() +
      scale_y_log10() +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black",
                                       angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"))
    
    png(paste0("./Figures/Microbiome/014_Sideness_wilcoxon/", Type, "/", Rank, "/", Tissue, "_",
               micr, "_boxplot_by_side.png"), res = 600, units = "in", width = 3, height = 3.5)
    plot(plot)
    dev.off()
  }
}

results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

results = results[order(results$p_val),]

results$Direction = NA
results$Direction[which(results$mean_right > results$mean_left)] = "Enriched in right"
results$Direction[which(results$mean_right <= results$mean_left)] = "Enriched in left"
results$ratio = results$mean_right / results$mean_left
results$log_ratio = log(results$ratio,10)

results$Proportion_difference = results$mean_right - results$mean_left

dir.create("./Analysis/Microbiome/014_Sideness_wilcoxon", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/014_Sideness_wilcoxon/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/014_Sideness_wilcoxon/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/014_Sideness_wilcoxon/",  Type, "/", Rank,
                            "/Wilcoxon_test_", Type, "_", Rank, "_", Tissue, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/014_Sideness_wilcoxon/",  Type, "/", Rank,
                                 "/Wilcoxon_test_", Type, "_", Rank, "_", Tissue, ".csv"),
          row.names = FALSE)

