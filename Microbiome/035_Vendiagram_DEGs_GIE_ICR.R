
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "VennDiagram"))

# Load data
load("./Analysis/Microbiome/033_GIE_wilcoxon_boxplot/Relative/Genus_full/_by_ICR_cluster_Relative_Genus_full.Rdata")
results_GIE = results
DEGs_GIE = results_GIE[which(results_GIE$p_val < 0.05),]

load("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot/Relative/Genus_full/_by_ICR_cluster_Relative_Genus_full.Rdata")
results_ICR = results
DEGs_ICR = results_ICR[which(results_ICR$p_val < 0.05),]

dir.create("./Figures/Microbiome/035_Vendiagram_DEGs_GIE_ICR", showWarnings = FALSE)

DEGs_GIE = DEGs_GIE[-grepl("uncul", DEGs_GIE$Name),]

venn.diagram(
  x = list(DEGs_GIE$Name, DEGs_ICR$Name),
  category.names = c("GIE" , "ICR"),
  filename = "./Figures/Microbiome/035_Vendiagram_DEGs_GIE_ICR/Venn1.png",
  output=TRUE
)


# IES and ICR
load("./Analysis/Microbiome/034_IES_Spearman/Relative/Genus_full/_by_ICR_cluster_Relative_Genus_full.Rdata")
DEGs_IES = results[which(results$p_val < 0.05),]

venn.diagram(
  x = list(DEGs_GIE$Name, DEGs_ICR$Name, DEGs_IES$Name),
  category.names = c("GIE" , "ICR", "IES"),
  filename = "./Figures/Microbiome/035_Vendiagram_DEGs_GIE_ICR/Venn2.png",
  output=TRUE
)

DEGs_Name_IES = DEGs_IES$Name
DEGs_Name_IES = DEGs_Name_IES[-which(DEGs_Name_IES %in% DEGs_ICR$Name)]
DEGs_Name_IES = DEGs_Name_IES[-which(DEGs_Name_IES %in% DEGs_GIE$Name)]

DEGs_IES = DEGs_IES[which(DEGs_IES$Name %in% DEGs_Name_IES),]

write.csv(DEGs_IES, file = "./Analysis/Microbiome/035_Vendiagram_results/DEGs_IES_11_only.csv", row.names = FALSE)
