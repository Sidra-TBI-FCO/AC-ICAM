

## Create Heatmap of selected genes
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
ha_genes = "DUX4_targets" #"Yao_DUX4_genes"

# Load data
load("./Processed_Data/External/Gene_collections/DUX4_Yao_genes.Rdata")
load("./Processed_Data/External/Gene_collections/DUX4_56_Gene_Symbol.Rdata")

# Create folders and log file
dir.create("./Figures/Trimmed_p",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/09_Heatmap", showWarnings = FALSE)

load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")

ICR_genes = c("IFNG", "TBX21", "CD8A", "CD8B", "IL12B", "STAT1", "IRF1",
              "CXCL9", "CXCL10", "CCL5",
              "GNLY", "PRF1", "GZMA", "GZMB", "GZMH",
              "CD274", "CTLA4", "FOXP3", "IDO1", "PDCD1")

load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

if(ha_genes == "ICR genes"){
  subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% ICR_genes, ])
  subset_RNAseq_log2 = subset_RNAseq_log2[,ICR_genes]
}

if(ha_genes == "BRCAness genes"){
  BRCAness_genes = BRCAness_df$Gene.Symbol
  BRCAness_genes = BRCAness_genes[which(BRCAness_genes %in% rownames(RNASeq.QN.LOG2))]
  subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% BRCAness_genes, ])
  subset_RNAseq_log2 = subset_RNAseq_log2[,BRCAness_genes]
}

if(ha_genes == "Yao_DUX4_genes"){
  Yao_DUX4_genes = Yao_DUX4_genes[which(Yao_DUX4_genes %in% rownames(RNASeq.QN.LOG2))]
  subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% Yao_DUX4_genes, ])
  subset_RNAseq_log2 = subset_RNAseq_log2[, Yao_DUX4_genes]
}

if(ha_genes == "DUX4_targets"){
  DUX4_targets_vector = DUX4_targets_vector[which(DUX4_targets_vector %in% rownames(RNASeq.QN.LOG2))]
  subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% DUX4_targets_vector, ])
  subset_RNAseq_log2 = subset_RNAseq_log2[, DUX4_targets_vector]
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

Expression.matrix = t(subset_RNAseq_log2[sample_order,])

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

png(filename = paste0("./Figures/Trimmed_p/Gene_of_interest_plots/09_Heatmap/", ha_genes, "_ComplexHeatmap.png"), res = 600,
    width = 5.5, height = 3.5, units = "in")
Heatmap(Expression.matrix.z, cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha, name = "Expression\n z score"
)
dev.off()

