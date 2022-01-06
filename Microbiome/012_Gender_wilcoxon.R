
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "plotrix"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Phylum" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
Tissue = "N"
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

abundance = get(paste("abundance_", Tissue, sep = ""))

df = data.frame(Patient_ID = unique(substring(colnames(rank_abundancies), 1, 3)),
                Abundance = NA,
                Gender = NA,
                stringsAsFactors = FALSE)

df$Gender = clinical_data$gender[match(df$Patient_ID, clinical_data$Patient_ID)]

results = data.frame(Name = rownames(rank_abundancies), p_val = NA, wilcox_statistic = NA,
                     mean_female = NA, mean_male = NA)

dir.create("./Figures/Microbiome/012_Wilcoxon_test_Gender", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/012_Wilcoxon_test_Gender/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/012_Wilcoxon_test_Gender/", Type, "/", Rank), showWarnings = FALSE)

i = 12
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[,micr][match(df$Patient_ID, rownames(abundance))]
  
  test = wilcox.test(x = df$Abundance[which(df$Gender == "FEMALE")], y = df$Abundance[which(df$Gender == "MALE")], paired = FALSE,
                     alternative = "two.sided")
  results$p_val[which(results$Name == micr)] = test$p.value
  results$wilcox_statistic[which(results$Name == micr)] = test$statistic
  results$mean_female[which(results$Name == micr)] = mean(df$Abundance[which(df$Gender == "FEMALE")])
  results$mean_male[which(results$Name == micr)] = mean(df$Abundance[which(df$Gender == "MALE")])
  
  if(Type == "Relative"){
    df$Abundance[which(df$Abundance == 0)] = 1e-4
  }
  
  melt = df
  melt$Gender = factor(melt$Gender, levels = c("FEMALE", "MALE"))
  
  if(is.na(test$p.value)){next}
  if(test$p.value < 0.05 | make_plot == "make_plot"){
    plot = ggplot(melt, aes(x = Gender, y = Abundance)) +
      geom_violin(trim = TRUE, aes(fill=Gender)) +
      #geom_jitter(position = "dodge",size=0.5) +
      #scale_fill_manual(values = c("ICR Low" = "blue", "ICR Medium" = "green", "ICR High" = "red")) +
      geom_boxplot(width=0.1, outlier.shape = NA) +
      ylab(paste0(micr, "\n", Type, " abundance")) +
      theme_bw() +
      stat_summary(fun.y = "median", colour = "grey", size = 1, geom = "point") +
      xlab("") +
      scale_y_log10() +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"),
            strip.background = element_blank(),
            strip.text = element_text(size = 15, colour = "black"))
    
    png(paste0("./Figures/Microbiome/012_Wilcoxon_test_Gender/", Type, "/", Rank, "/", Tissue, "_",
               micr, "_boxplot_by_gender.png"), res = 600, units = "in", width = 4, height = 3.5)
    plot(plot)
    dev.off()
  }
}

results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

results = results[order(results$p_val),]

results$Direction = NA
results$Direction[which(results$mean_female > results$mean_male)] = "Enriched in female"
results$Direction[which(results$mean_female <= results$mean_male)] = "Enriched in male"
results$ratio = results$mean_male / results$mean_female
results$log_ratio = log(results$ratio,10)

results$Proportion_difference = results$mean_male - results$mean_female

dir.create("./Analysis/Microbiome/012_Wilcoxon_test_Gender", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/012_Wilcoxon_test_Gender/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/012_Wilcoxon_test_Gender/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/012_Wilcoxon_test_Gender/",  Type, "/", Rank,
                            "/Wilcoxon_test_", Type, "_", Rank, "_", Tissue, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/012_Wilcoxon_test_Gender/",  Type, "/", Rank,
                                 "/Wilcoxon_test_", Type, "_", Rank, "_", Tissue, ".csv"),
          row.names = FALSE)
