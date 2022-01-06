
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "ComplexHeatmap"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Phylum" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)

# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Analysis
rank_abundancies = get(paste0(Rank, "_abundance"))

Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"
MANTIS = MANTIS[-which(MANTIS$Tissue == "LM"),]

annotation_df = data.frame(Sample_ID = colnames(rank_abundancies), ICR_cluster = NA, Tissue_type = NA, CMS = NA, MSI = NA)
annotation_df$Tissue_type = substring(annotation_df$Sample_ID, 4, 4)
annotation_df$ICR_cluster = table_cluster_assignment$ICR_HML[match(substring(annotation_df$Sample_ID, 1, 3), substring(rownames(table_cluster_assignment), 1, 3))]
annotation_df$CMS = Rfcms$RF.predictedCMS[match(substring(annotation_df$Sample_ID, 1, 3), substring(rownames(Rfcms), 1, 3))]
annotation_df$MSI = MANTIS$MSI[match(substring(annotation_df$Sample_ID, 1, 3), MANTIS$Patient_ID)]

annotation_df$ICR_cluster = factor(annotation_df$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

sample_order = annotation_df$Sample_ID

rank_abundancies = as.matrix(rank_abundancies)
mode(rank_abundancies)

rank_abundancies = rank_abundancies[, sample_order]
rank_abundancies = rank_abundancies 

# z-score abundance matrix
rank_abundancies.z = rank_abundancies +1
for(j in 1: nrow(rank_abundancies.z))  {
  rank_abundancies.z[j,] = (rank_abundancies[j,]-mean(rank_abundancies[j,]))/sd(rank_abundancies[j,]) # z-score the matrix
}

# Annotation for heatmap
ha = HeatmapAnnotation(Tissue_type = annotation_df$Tissue_type,
                       ICR = annotation_df$ICR_cluster,
                       CMS = annotation_df$CMS,
                       MSI = annotation_df$MSI,
                       col = list(ICR = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  MSI = c("MSI-H" = "purple", "MSS" = "white", "Not determined" = "lightgrey"),
                                  CMS = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                            "mixed" = "lightgrey"),
                                  Tissue_type = c("N" = "purple", "T" = "orange")))

dir.create("./Figures/Microbiome/011_Heatmap_Relative_abundance", showWarnings = FALSE)

# Plot Heatmap
png(paste0("./Figures/Microbiome/011_Heatmap_Relative_abundance/", Rank, "_", Type, "_Heatmap.png"), res = 600, 
    width = 8, height = 6, units = "in")
Heatmap(rank_abundancies.z, cluster_rows = FALSE, cluster_columns = TRUE,
        show_column_names = FALSE, top_annotation = ha, name = "Relative abundance\n z score"
)
dev.off()

