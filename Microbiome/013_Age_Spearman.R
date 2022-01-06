
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Phylum" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
Tissue = "T"
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
                Age = NA,
                Age_cat = NA,
                stringsAsFactors = FALSE)

df$Age = clinical_data$age_at_initial_pathologic_diagnosis[match(df$Patient_ID, clinical_data$Patient_ID)]
df$Age = as.numeric(df$Age)
df$Age_cat[which(df$Age < 50)] = "<50"
df$Age_cat[which(df$Age >=50 & df$Age < 65)] = "50-65"
df$Age_cat[which(df$Age >=65 & df$Age < 75)] = "65-75"
df$Age_cat[which(df$Age >=75)] = ">=75"

table(df$Age_cat)

df$Age_cat = factor(df$Age_cat, levels = c("<50", "50-65", "65-75", ">=75"))

dir.create("./Figures/Microbiome/013_Age_Spearman", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/013_Age_Spearman/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/013_Age_Spearman/", Type, "/", Rank), showWarnings = FALSE)

results = data.frame(Name = rownames(rank_abundancies), p_val = NA, rho = NA,
                     mean_50 = NA,
                     mean_50_65 = NA, 
                     mean_65_75 = NA,
                     mean_75 = NA)
i= 1
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[, micr][match(df$Patient_ID, rownames(abundance))]
  
  df_calc = df
  levels(df_calc$Age_cat) = c("1", "2", "3", "4")
  df_calc$Age_cat = as.numeric(df_calc$Age_cat)
  cor = cor.test(df_calc$Age_cat, df_calc$Abundance, method = "spearman")
  results$p_val[which(results$Name == micr)] = cor$p.value
  results$rho[which(results$Name == micr)] = cor$estimate
  results$mean_50[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$Age_cat == "1")])
  results$mean_50_65[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$Age_cat == "2")])
  results$mean_65_75[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$Age_cat == "3")])
  results$mean_75[which(results$Name == micr)] = mean(df_calc$Abundance[which(df_calc$Age_cat == "4")])
  
  if(Type == "Relative"){
    df$Abundance[which(df$Abundance == 0)] = 1e-4
  }
  
  if(is.na(cor$p.value)){next}
  if(cor$p.value < 0.05 | make_plot == "make_plot"){
    plot = ggplot(df, aes(x = Age_cat, y = Abundance)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_jitter(width = 0.1, size = 0.8, aes(color = Mutation_cat)) +
      geom_point(size = 0.3) +
      theme_bw() +
      scale_y_log10() +
      ylab(paste0(micr, "\n", Type, " abudance")) +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black",
                                       angle = 90, vjust = 1, hjust = 0.5),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"))
    
    png(paste0("./Figures/Microbiome/013_Age_Spearman/", Type, "/", Rank, "/", Tissue, "_", micr,
               "_by_age_category_boxplot.png"), res = 600, width = 4, height = 4, units = "in")
    plot(plot)
    dev.off()
  }
}

results = results[order(results$p_val),]
results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

dir.create(paste0("./Analysis/Microbiome/013_Age_Spearman"), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/013_Age_Spearman/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/013_Age_Spearman/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/013_Age_Spearman/",  Type, "/", Rank,
                            "/", Tissue, "_Spearman_by_age_category_", Type, "_", Rank, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/013_Age_Spearman/",  Type, "/", Rank,
                                 "/", Tissue, "_Spearman_by_age_category_", Type, "_", Rank, ".csv"),
          row.names = FALSE)
