#################################################################
###
### This script creates heatmaps for ICR genes using EDASeq 
### normalized RNASeq data. Samples are ordered by ICR 
### score and a heatmap is created for ICR genes including 
### ICR cluster allocation for k = 3, 4 and 6 and the proposed
### HML classfication.
### Heatmaps are saved as png files at location: 
### 
###
#################################################################

## Create Heatmap of TCGA RNASeq ICR genes
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
Log_file = paste0("./Logfiles/ICR Heatmaps/ICR_heatmaps_",gsub(":",".",gsub(" ","_",date())),".txt")
with_metastasis = ""
ha_genes = "ICR genes" #"BRCAness genes"

# Load data
load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))
BRCAness_df = read.csv("./Processed_Data/External/Gene_collections/Konstantinopoulus_BRCAness_Signature_Supplementary_Table_1.csv", stringsAsFactors = FALSE)

# Create folders and log file
dir.create("./Figures/Trimmed_p",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Heatmaps", showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Heatmaps/ICR Heatmaps", showWarnings = FALSE)
dir.create("./Logfiles/ICR Heatmaps", showWarnings = FALSE)

load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")

ICR_genes = c("IFNG", "TBX21", "CD8A", "CD8B", "IL12B", "STAT1", "IRF1",
              "CXCL9", "CXCL10", "CCL5",
              "GNLY", "PRF1", "GZMA", "GZMB", "GZMH",
              "CD274", "CTLA4", "FOXP3", "IDO1", "PDCD1")

if(with_metastasis == "with_metastasis"){
  load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms_with_metastasis.Rdata")
}else{
  load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
}

load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

if(with_metastasis == "with_metastasis"){
  load(paste0("./Analysis/Trimmed_p/ICR Consensus Clustering/With_metastasis/JSREP_ICR_cluster_assignment_k2-6.Rdata"))
  table_cluster_assignment = table_cluster_assignment[-grep("LM", rownames(table_cluster_assignment)),]
  table_cluster_assignment = table_cluster_assignment[-grep("B", rownames(table_cluster_assignment)),]
  table_cluster_assignment = table_cluster_assignment[-grep("C", rownames(table_cluster_assignment)),]
  RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, which(colnames(RNASeq.QN.LOG2) %in% rownames(table_cluster_assignment))]
}else{
  load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
}

#load("./Analysis/CMS/Rfcms.Rdata")

if(ha_genes == "ICR genes"){
  ICR_subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% ICR_genes, ])
  ICR_subset_RNAseq_log2 = ICR_subset_RNAseq_log2[,ICR_genes]
}

if(ha_genes == "BRCAness genes"){
  BRCAness_genes = BRCAness_df$Gene.Symbol
  BRCAness_genes = BRCAness_genes[which(BRCAness_genes %in% rownames(RNASeq.QN.LOG2))]
  ICR_subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% BRCAness_genes, ])
  ICR_subset_RNAseq_log2 = ICR_subset_RNAseq_log2[,BRCAness_genes]
}

table_cluster_assignment$ICR_HML = factor(table_cluster_assignment$ICR_HML, levels = c("ICR High", "ICR Medium", "ICR Low"))

#table_cluster_assignment = table_cluster_assignment[order(table_cluster_assignment$ICRscore),]
table_cluster_assignment = table_cluster_assignment[order(table_cluster_assignment$ICR_HML, decreasing = TRUE),]
table_cluster_assignment$Sample_ID = rownames(table_cluster_assignment)
table_cluster_assignment$CMS = Rfcms$RF.predictedCMS[match(table_cluster_assignment$Sample_ID, rownames(Rfcms))]
table_cluster_assignment$CMS[which(is.na(table_cluster_assignment$CMS))] = "Not predicted"
table_cluster_assignment$MSI = MANTIS$MSI[match(substring(rownames(table_cluster_assignment), 1, 3),
                                                MANTIS$Patient_ID)]
table_cluster_assignment$MSI[which(is.na(table_cluster_assignment$MSI))] = "Not determined"

table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(RNASeq.QN.LOG2)),]
sample_order = rownames(table_cluster_assignment)

Expression.matrix = t(ICR_subset_RNAseq_log2[sample_order,])

# z-score Expression.matrix
Expression.matrix.z = Expression.matrix
for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}


table_cluster_assignment$CMS = factor(table_cluster_assignment$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "Not predicted"))
table_cluster_assignment$MSI = factor(table_cluster_assignment$MSI, levels = c("MSI-H", "MSS", "Not determined"))


ha = HeatmapAnnotation(ICR = table_cluster_assignment$ICR_HML,
                       `CMS` = table_cluster_assignment$CMS,
                       MSI = table_cluster_assignment$MSI,
                       col = list(ICR = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  MSI = c("MSI-H" = "purple", "MSS" = "white", "Not determined" = "lightgrey"),
                                  `CMS` = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                                    "Not predicted" = "lightgrey")))

png(filename = paste0("./Figures/Trimmed_p/Heatmaps/ICR Heatmaps/v3_", ha_genes, "_ComplexHeatmap_", with_metastasis, ".png"), res = 600,
    width = 5.5, height = 3.5, units = "in")
Heatmap(Expression.matrix.z, cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha, name = "Expression\n z score"
)
dev.off()
