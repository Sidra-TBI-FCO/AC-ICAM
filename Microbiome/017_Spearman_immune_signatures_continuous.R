

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
subset = "" # "nonhypermutated" or "hypermutated" or "Right sided" or "Left sided"

# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Clean_TRB_clonality_341_patients.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/Deconvolution_and_GSEA/Selected.pathways_ES.Rdata")
ES_path = ES

# Analysis
rank_abundancies = get(paste0(Rank, "_abundance"))
abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)

abundance_T = as.data.frame(abundance_T)
abundance_T$ICR_cluster = table_cluster_assignment$ICR_HML[match(rownames(abundance_T),
                                                                 substring(rownames(table_cluster_assignment), 1, 3))]
abundance_T$ICRscore = table_cluster_assignment$ICRscore[match(rownames(abundance_T),
                                                               substring(rownames(table_cluster_assignment), 1, 3))]

rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1, 4)

# Add TCR data
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
immune_sig_df$TCR_Clonality_ImmunoSeq = TCR_Overview$productive_clonality[match(substring(rownames(immune_sig_df), 1, 3), TCR_Overview$Patient_ID)]
MiXCR = MiXCR[which(MiXCR$Tissue == "T-P"),]
immune_sig_df$MiXCR_TRB_Clonality = MiXCR$MiXCR_Clonality[match(substring(rownames(immune_sig_df), 1, 3), MiXCR$Patient_ID)]
immune_sig_df$MiXCR_TRB_Clonality[which(immune_sig_df$MiXCR_TRB_Clonality == -Inf)] = NA

# Add mutational load
immune_sig_df$Mutational_load = frequency_df$Nonsilent_mutational_burden_per_Mb[match(substring(rownames(immune_sig_df), 1, 3),
                                                                                      frequency_df$Patient_ID)]

# Select right rows (patients)
immune_sig_df = immune_sig_df[which(substring(rownames(immune_sig_df), 1, 3) %in%
                                      rownames(abundance_T)),]
#immune_sig_df = immune_sig_df[, 1:3]

# Add mutation categories for subset
frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

# Add sideness categories for subset
clinical_data$Primary_tumor_side = NA
clinical_data$Primary_tumor_side[which(clinical_data$tumour_anatomic_site %in% c("ceceum", "colon ascendens", "flexura hepatica",
                                                                                 "colon transversum"))] = "Right sided"
clinical_data$Primary_tumor_side[which(clinical_data$tumour_anatomic_site %in% c("flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                                 "rectosigmoideum"))] = "Left sided"
clinical_data$Primary_tumor_side = factor(clinical_data$Primary_tumor_side, levels = c("Right sided", "Left sided"))

dir.create("./Figures/Microbiome/017_Spearman_continous", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/017_Spearman_continous/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/017_Spearman_continous/", Type, "/", Rank), showWarnings = FALSE)

# Analysis
df_plot = data.frame(Patient_ID = rownames(abundance_T),
                     ICR_cluster = abundance_T$ICR_cluster,
                     Enrichment = NA,
                     Mutation_cat = NA,
                     Side = NA,
                     Abundance = NA)
df_plot$ICR_cluster = factor(df_plot$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
df_plot$Mutation_cat = frequency_df$Mutation_cat[match(df_plot$Patient_ID, frequency_df$Patient_ID)]
df_plot$Side = clinical_data$Primary_tumor_side[match(df_plot$Patient_ID, clinical_data$Patient_ID)]

if(subset %in% c("nonhypermutated", "hypermutated")){
  df_plot = df_plot[which(df_plot$Mutation_cat == subset),]
}

if(subset %in% c("Left sided", "Right sided")){
  df_plot = df_plot[which(df_plot$Side == subset),]
}

results_all = data.frame(Name = NA, signature = NA, p_val = NA, rho = NA)

j = 8 # 22 (= TREM1 data) # 78 (= attractor G) # 8 (PDL1 data) # 63 (GRANS PCA)
for(j in 1:ncol(immune_sig_df)){
  sig = colnames(immune_sig_df)[j]
  df_plot$Enrichment = immune_sig_df[, sig][match(df_plot$Patient_ID, substring(rownames(immune_sig_df), 1, 3))]
  
  results = data.frame(Name = rownames(rank_abundancies),
                       signature = NA, p_val = NA, rho = NA)
  
  i= 371 # 371 (= Ruminococcus 1) # 374 (= subdo) 414 (=Selenomonas) 277 (Fusicatenibacter)
  for (i in 1:nrow(rank_abundancies)){
    micr = rownames(rank_abundancies)[i]
    df_plot$Abundance = abundance_T[, micr][match(df_plot$Patient_ID, rownames(abundance_T))]
    
    df_calc = df_plot
    df_plot_edit = df_plot
    df_calc$Enrichment = as.numeric(df_calc$Enrichment)
    cor = cor.test(df_calc$Enrichment, df_calc$Abundance, method = "spearman")
    results$p_val[which(results$Name == micr)] = cor$p.value
    results$rho[which(results$Name == micr)] = cor$estimate
    
    if(Type == "Relative"){
      df_plot_edit$Abundance[which(df_plot$Abundance == 0)] = 1e-4
    }
    
    if(is.na(cor$p.value)){next}
    if(make_plot == "make_plot" & cor$p.value < 0.001){
      plot = ggplot(df_plot_edit, aes(x = Enrichment, y = Abundance)) +
        geom_point(size = 0.3) +
        theme_bw() +
        ylab(paste0(gsub(".*\\D5__", "", micr))) +
        xlab(gsub(".*\\|  ", "", sig)) +
        #scale_y_log10() +
        scale_y_log10(limits = c(0.0001, 0.2)) +
        #ylim(0.0001, 0.075) +
        #facet_grid(.~Mutation_cat) +
        #stat_cor(method = "spearman", size = 5, position = "upperright") +
        geom_smooth(method="lm") +
        theme(axis.text.y = element_text(size = 20, colour = "black"),
              axis.text.x = element_text(size = 20, colour = "black"),
              axis.title.x = element_text(size = 20, colour = "black"),
              axis.title.y = element_text(size = 20, colour = "black"),
              strip.background = element_blank(),
              strip.text = element_text(size = 15, colour = "black"),
              aspect.ratio = 1/1)
      
      png(paste0("./Figures/Microbiome/017_Spearman_continous/", Type, "/", Rank, "/facet_", micr, "_", subset,
                 "_by_", sig, "_scatterplot.png"), res = 600, width = 4, height = 4, units = "in")
      plot(plot)
      dev.off()
    }
  }
  #results = results[-which(is.na(results$p_val)),]
  results$signature = sig
  results_all = rbind(results_all, results)
}

results_all = results_all[order(results_all$p_val),]

dir.create(paste0("./Analysis/Microbiome/017_Spearman_continous"), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/017_Spearman_continous/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/017_Spearman_continous/",  Type, "/", Rank), showWarnings = FALSE)

save(results_all, file = paste0("./Analysis/Microbiome/017_Spearman_continous/",  Type, "/", Rank,
                                "/v5_April_2021_", subset, "_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank, ".Rdata"))
write.csv(results_all, file = paste0("./Analysis/Microbiome/017_Spearman_continous/",  Type, "/", Rank,
                                     "/v5_April_2021_", subset ,"_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank, ".csv"),
          row.names = FALSE)
