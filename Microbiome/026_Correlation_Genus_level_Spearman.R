
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("gplots", "ComplexHeatmap", "circlize") #"color"
ipak(required.packages)

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full" 
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
Tissue = "T"

# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))
rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)
abundance_N = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "N")]
colnames(abundance_N) = substring(colnames(abundance_N), 1, 3)
abundance_N = t(abundance_N)

abundance = get(paste("abundance_", Tissue, sep = ""))

abundance = abundance[,-which(colSums(abundance) == 0)]

# Calculate correlation
mat_cor = cor(abundance, method = "spearman")
mat_cor = as.matrix(mat_cor)

dir.create("./Figures/Microbiome/026_Spearman_correlation_plot", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/026_Spearman_correlation_plot/026_Spearman_correlation_heatmap.2_", 
           Type, "_", Rank, "_", Tissue, ".png"), res = 600, width = 13, height = 13, units = "in")
hm = heatmap.2(mat_cor,
               # dendrogram control
               Rowv = TRUE,
               distfun = dist,
               hclustfun = hclust,
               dendrogram = c("both"),
               col = "bluered",
               trace = "none",
               margins=c(10,10),
               key = "T", density = "none")
dev.off()


#
#column_ha = HeatmapAnnotation(Immune_trait_category = annotation$Immune_trait_category,
                  #            show_annotation_name = FALSE,
                   #           col = list(Immune_trait_category = c("Expression Signature " = "#b4e2c5", "Leukocyte Subset ES (ConsensusTME) " = "#eda8ff",
                    #                                               "Attractor Metagene " = "#ffbe63")),
                     #         show_legend = FALSE)


col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "#ff2626"))


png(paste0("./Figures/Microbiome/026_Spearman_correlation_plot/026_Spearman_correlation_ComplexHeatmap_", 
           Type, "_", Rank, "_", Tissue, ".png"), res = 600, width = 8, height = 8, units = "in")
HM = Heatmap(mat_cor,
             row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5),
             column_title_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
             cluster_rows = TRUE, cluster_columns = T,
             show_column_names = T, 
             #top_annotation = column_ha, 
             #col = col_fun,
             show_heatmap_legend = TRUE,
             #col = color,
             # heatmap_legend_param =list(title_gp=gpar(fontsize=10, fontface="bold"),legend_width=unit(8,"cm"),legend_position = "left"),
             #row_names_max_width = unit(5, "cm")
             # theme(legend.position = "none")
             row_names_max_width = unit(6, "in")
)
HM = draw(HM, heatmap_legend_side = "left")
dev.new()
dev.off()

order = rownames(mat_cor)[row_order(HM)]
dir.create("./Analysis/Microbiome/026_Correlation_plot_Order_Genus", showWarnings = FALSE)
save(order, file = "./Analysis/Microbiome/026_Correlation_plot_Order_Genus/026_order_genus.Rdata")
