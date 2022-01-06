
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("GenVisR", "dplyr", "ggplot2", "maftools", "ComplexHeatmap",
                      "data.table", "tidyr", "stringr", "cowplot", "tidyverse")
ipak(required.packages)

# Set parameters
selection = "cancer_driver_QGP" # "cancer_driver_501" or "cancer_driver_QGP" or "All" or "artifacts_removed"

# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")
load("./Processed_Data/External/Gene_collections/501_genes_Michele.Rdata")
driver_genes = data.frame(gene = genes_Michele)
driver_genes$gene = as.character(driver_genes$gene)
QGPC_genes = data.frame(gene = QGPC_genes)
QGPC_genes$gene = as.character(QGPC_genes$gene)

# Gene count ranking
finalMafFiltered$Tumor_Sample_Barcode = as.character(finalMafFiltered$Tumor_Sample_Barcode)
MAF_gene_Count= finalMafFiltered %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(count=n()) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%     
  summarise(count=n()) %>% group_by(Hugo_Symbol) %>% 
  summarise(count=n()) %>% arrange(desc(count)) 

MAF_gene_Count$Percentage = round((MAF_gene_Count$count / 281) * 100, 1)
write.csv(MAF_gene_Count, file = "./Analysis/WES/001_Oncoplot/Gene_Count_MAF.csv")

MAF_gene_Count = MAF_gene_Count %>% 
  filter(count >= 2)  %>% select(Hugo_Symbol)

MAF_gene_Count_GENE = MAF_gene_Count$Hugo_Symbol

# Intersection with cancer driver genes (501 from Michele)
genes = MAF_gene_Count
if(selection == "cancer_driver_501"){
  genes = genes %>% inner_join(driver_genes, by=c("Hugo_Symbol" = "gene"))
}
if(selection == "cancer_driver_QGP"){
  genes = genes %>% inner_join(QGPC_genes, by=c("Hugo_Symbol" = "gene"))
}

if(selection == "artifacts_removed"){
  artifacts = c("LOC","ENS","FAM","GOL","PRA","NBP","POT","DEF","MUC","KRT","WAS","ANK","TRI","FRG",paste0("OR",1:9))
  genes = genes[-which(substring(genes$Hugo_Symbol, 1, 3) %in% artifacts),]
  artifacts2 = c("PLIN","CELA","SRA1")
  genes = genes[-which(substring(genes$Hugo_Symbol, 1, 4) %in% artifacts2),]
  artifacts_genes = c("ATXN1","PBRM1","ZNF814","MSH3","TTN","USH2A")
  genes = genes[-which(genes$Hugo_Symbol %in% artifacts_genes),]
}
genes = genes$Hugo_Symbol

# Plotting
## plotting only for diffGene for only IBC sample:
silu= read.maf(maf = finalMafFiltered)

#One can use any colors, here in this example color palette from RColorBrewer package is used
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

dir.create("./Figures/WES/001_Waterfall_and_Oncoplot", showWarnings = FALSE)
png(paste0("./Figures/WES/001_Waterfall_and_Oncoplot/Oncoplot_primary_tumor_clean_", 
           selection,".png"), res = 600,
    height = 6, width = 8, units = "in")
oncoplot(maf=silu, genes = genes[1:50],
         #writeMatrix = TRUE,
         fontSize = 0.4, 
         gene_mar = 5,
         #colors = vc_cols,
         draw_titv = TRUE,
         removeNonMutated = FALSE)
dev.off()

finalMafFiltered = data.frame(finalMafFiltered)

png(paste0("./Figures/WES/001_Waterfall_and_Oncoplot/Waterfall_plot_primary_tumor_clean", 
           selection, ".png"), res = 600,
    height = 5, width = 10, units = "in")
waterfall(finalMafFiltered, 
          fileType = "MAF", 
          geneOrder = genes[1:25], 
          mainXlabel = FALSE, 
          #variant_class_order = most_deleterious,
          plotMutBurden = FALSE)
dev.off()

############
png("./Figures/WES/001_Waterfall_and_Oncoplot/Plot_mafSummary_primary_tumor_clean.png", res = 600,
    height = 6, width = 7, units = "in")
plotmafSummary(maf = silu, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

