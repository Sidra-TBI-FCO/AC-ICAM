
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Phylum" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
Tissue = "N"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"
option = "option 3" # "option 1" (for <40 compared to >= 40)

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
                Age = NA,
                Age_cat = NA,
                stringsAsFactors = FALSE)

df$Age = clinical_data$age_at_initial_pathologic_diagnosis[match(df$Patient_ID, clinical_data$Patient_ID)]
df$Age = as.numeric(df$Age)

if(option == "option 1"){
  df$Age_cat[which(df$Age < 40)] = "Group1"
  df$Age_cat[which(df$Age >=40)] = "Group2"
  table(df$Age_cat)
}

if(option == "option 2"){
  df$Age_cat[which(df$Age < 45)] = "Group1"
  df$Age_cat[which(df$Age >=45)] = "Group2"
  table(df$Age_cat)
}

if(option == "option 3"){
  df$Age_cat[which(df$Age < 45)] = "Group1"
  df$Age_cat[which(df$Age >= 85)] = "Group2"
  df = df[-which(is.na(df$Age_cat)),]
  table(df$Age_cat)
}

df$Age_cat = factor(df$Age_cat, levels = c("Group1", "Group2"))

dir.create("./Figures/Microbiome/013.2_Age", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/013.2_Age/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/013.2_Age/", Type, "/", Rank), showWarnings = FALSE)

results = data.frame(Name = rownames(rank_abundancies), p_val = NA,
                     mean_group1 = NA,
                     mean_group2 = NA)
i= 1
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[, micr][match(df$Patient_ID, rownames(abundance))]
  
  test = wilcox.test(x = df$Abundance[which(df$Age_cat == "Group1")], y = df$Abundance[which(df$Age_cat == "Group2")], paired = FALSE,
                     alternative = "two.sided")
  results$p_val[which(results$Name == micr)] = test$p.value
  results$mean_group1[which(results$Name == micr)] = mean(df$Abundance[which(df$Age_cat == "Group1")])
  results$mean_group2[which(results$Name == micr)] = mean(df$Abundance[which(df$Age_cat == "Group2")])
  
  if(Type == "Relative"){
    df$Abundance[which(df$Abundance == 0)] = 1e-4
  }
  
  if(is.na(test$p.value)){next}
  if(test$p.value < 0.05 | make_plot == "make_plot"){
    plot = ggplot(df, aes(x = Age_cat, y = Abundance)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_jitter(width = 0.1, size = 0.8, aes(color = Mutation_cat)) +
      geom_point(size = 0.3) +
      theme_bw() +
      scale_y_log10() +
      ylab(paste0(micr, "\n", Type, " abundance")) +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black",
                                       angle = 90, vjust = 0.5, hjust = 0.5),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"))
    
    png(paste0("./Figures/Microbiome/013.2_Age/", Type, "/", Rank, "/", Tissue, "_", micr,
               "_by_age_category_boxplot.png"), res = 600, width = 4, height = 4, units = "in")
    plot(plot)
    dev.off()
  }
}

results = results[order(results$p_val),]
results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

dir.create(paste0("./Analysis/Microbiome/013.2_Age"), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/013.2_Age/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/013.2_Age/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/013.2_Age/",  Type, "/", Rank,
                            "/", Tissue, "_Spearman_by_age_category_", Type, "_", Rank, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/013.2_Age/",  Type, "/", Rank,
                                 "/", Tissue, "_Spearman_by_age_category_", Type, "_", Rank, ".csv"),
          row.names = FALSE)
