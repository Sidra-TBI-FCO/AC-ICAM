
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "plotrix"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
Tissue = "N" # "T" or "N" or "Difference"
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

df$Anatomical_location = factor(df$Anatomical_location, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                   "colon transversum", "flexura lienalis", "colon descendens", 
                                                                   "colon sigmoideum", "rectosigmoideum"))

table(df$Anatomical_location)
levels(df$Anatomical_location) = c("1", "2", "3", "4", "5", "6", "7", "8")
df$Anatomical_location = as.numeric(df$Anatomical_location)
table(df$Anatomical_location)

results = data.frame(Name = rownames(rank_abundancies), cor = NA, p_val = NA)

dir.create("./Figures/Microbiome/014.2_Spearman", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/014.2_Spearman/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/014.2_Spearman/", Type, "/", Rank), showWarnings = FALSE)

i = 12
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[,micr][match(df$Patient_ID, rownames(abundance))]
  
  cor = cor.test(df$Abundance, df$Anatomical_location, method = "spearman")
  results$p_val[which(results$Name == micr)] = cor$p.value
  results$cor[which(results$Name == micr)] = cor$estimate
  
  if(Type == "Relative" & Tissue %in% c("T", "N")){
    df$Abundance[which(df$Abundance == 0)] = 1e-4
  }
  
  melt = df
  
  if(is.na(cor$p.value)){next}
  if(cor$p.value < 0.05 | make_plot == "make_plot"){
    melt$Anatomical_location = factor(melt$Anatomical_location)
    levels(melt$Anatomical_location) = c("ceceum", "colon ascendens", "flexura hepatica", 
                                         "colon transversum", "flexura lienalis", "colon descendens", 
                                         "colon sigmoideum", "rectosigmoideum")
    label = gsub(".*\\D_5__", "", micr)
    
    plot = ggplot(melt, aes(x = Anatomical_location, y = Abundance)) +
      geom_boxplot(outlier.shape = NA) +
      geom_line(aes(group = Patient_ID), color = alpha("lightgrey", 0.5)) +
      geom_point(size = 0.8) +
      ylab(paste0(label)) +
      theme_bw() +
      scale_y_log10() +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black",
                                       angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"))
    
    png(paste0("./Figures/Microbiome/014.2_Spearman/", Type, "/", Rank, "/", Tissue, "_",
               micr, "_boxplot_anatomic location.png"), res = 600, units = "in", width = 3, height = 3.5)
    plot(plot)
    dev.off()
  }
}

results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

results = results[order(results$p_val),]

dir.create("./Analysis/Microbiome/014.2_Anatomic_location_Spearman", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/014.2_Anatomic_location_Spearman/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/014.2_Anatomic_location_Spearman/",  Type, "/", Rank), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Microbiome/014.2_Anatomic_location_Spearman/",  Type, "/", Rank,
                            "/Spearman_", Type, "_", Rank, "_", Tissue, ".Rdata"))
write.csv(results, file = paste0("./Analysis/Microbiome/014.2_Anatomic_location_Spearman/",  Type, "/", Rank,
                                 "/Spearman_", Type, "_", Rank, "_", Tissue, ".csv"),
          row.names = FALSE)

