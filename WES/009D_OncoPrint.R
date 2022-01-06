
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ggplot2", "ComplexHeatmap")
ipak(required.packages)

# Set parameters
version = "D"
ha_genes = "DDR" # "DDR" or "BRCAness"
subset = "hypermutated" # "all_patients" "hypermutated"
order_by_cluster_heatmap = "" #"order_by_cluster_heatmap"

# Load data
load(paste0("./Analysis/WES/008_matrix_for_Oncoplot/008", version, "_matrix_for_oncoplot.Rdata"))
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
DDR = read.csv("./Processed_Data/External/Gene_collections/DDR_pathway_members.csv", stringsAsFactors = FALSE, check.names = FALSE)
BRCAness = read.csv("./Processed_Data/External/Gene_collections/Konstantinopoulus_BRCAness_Signature_Supplementary_Table_1.csv", stringsAsFactors = FALSE)
genes = DDR$Gene[which(DDR$NEW.ANNOTATION_CONSENSUS_Sept_11 %in% c("BRCA1", "BRCA2", "FA", "HR-NON BRCA"))]

DDR[DDR==""] = NA

if(ha_genes == "mat_genes"){
  genes = rownames(mat)
}
# Analysis

# Define functions for OncoPrint
col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "black")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#e6e6e6", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  }
)

#column_title = "Sidra-LUMC cohort \n BRCA1, BRCA2, FA, NON-BRCA HR somatic mutations"
column_title = "Sidra-LUMC cohort"
heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "MUT"), 
                            labels = c("Deep deletion", "Amplification", "Mutation"))

# Prepare for annotation
MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Rfcms = Rfcms[which(substring(rownames(Rfcms), 1, 4) %in% colnames(mat)),]
Rfcms$RF.predictedCMS = as.character(Rfcms$RF.predictedCMS)
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"

df = data.frame(Sample_ID = substring(rownames(Rfcms), 1, 4), CMS = Rfcms$RF.predictedCMS, 
                ICR = NA, MSI = NA, Mutational_load = NA)
rownames(df) = df$Sample_ID
df$MSI = MANTIS$MSI[match(df$Sample_ID, MANTIS$Sample_ID)]
df$Mutational_load = frequency_df$Non_silent_Mutation_frequency[match(substring(df$Sample_ID, 1, 3), frequency_df$Patient_ID)]
df$Mutational_load_per_MB = frequency_df$Nonsilent_mutational_burden_per_Mb[match(substring(df$Sample_ID, 1, 3), frequency_df$Patient_ID)]
df$ICR = table_cluster_assignment$ICR_HML[match(df$Sample_ID, substring(rownames(table_cluster_assignment), 1, 4))]
df$ICRscore = table_cluster_assignment$ICRscore[match(df$Sample_ID, substring(rownames(table_cluster_assignment), 1, 4))]

df$Mutation_cat = NA
df$Mutation_cat[which(df$Mutational_load_per_MB <= 12)] = "nonhypermutated"
df$Mutation_cat[which(df$Mutational_load_per_MB > 12)] = "hypermutated"

df$ICR = factor(df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))
df$CMS = factor(df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "mixed"))
df$MSI = factor(df$MSI, levels = c("MSI-H", "MSS"))
df$Mutation_cat = factor(df$Mutation_cat, levels = c("nonhypermutated", "hypermutated"))

# Define order based on mutational load
df = df[order(df$ICRscore),]
df = df[order(df$Mutational_load),]
df = df[order(df$CMS),]
df = df[order(df$ICR),]
df = df[order(df$Mutation_cat),]

if(subset == "all_patients"){}else{
  df = df[which(df$Mutation_cat == subset),]
}
sample_order = rownames(df)
mat = mat[,sample_order]

# Expression heatmap
colnames(RNASeq.QN.LOG2) = substring(colnames(RNASeq.QN.LOG2), 1, 4)
#Expression.matrix = RNASeq.QN.LOG2[rownames(mat),sample_order]

if(ha_genes == "DDR"){
  genes = unique(DDR$`Homology-dependent recombination (HDR)`)
  genes = genes[!is.na(genes)]
  genes = genes[which(genes %in% rownames(RNASeq.QN.LOG2))]
  # remove flat gene
  genes = genes[-which(genes == "SPO11")]
  #genes = DDR$Gene[which(DDR$NEW.ANNOTATION_CONSENSUS_Sept_11 %in% c("BRCA1", "BRCA2", "FA", "HR-NON BRCA"))]
  #genes = c(genes[c(1:2, 6)], "POLQ", genes[c(3:5, 7:20)])
}
if(ha_genes == "BRCAness"){
  genes = BRCAness$Gene.Symbol[which(BRCAness$Gene.Symbol %in% rownames(RNASeq.QN.LOG2))]
}

Expression.matrix = RNASeq.QN.LOG2[genes,sample_order]

# z-score Expression.matrix
Expression.matrix.z = Expression.matrix
for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

if(order_by_cluster_heatmap == "order_by_cluster_heatmap"){
  heatmap = Heatmap(Expression.matrix.z, cluster_rows = TRUE, cluster_columns = TRUE,
                    show_column_names = FALSE, name = "Expression\n z score")
  
  sample_order = colnames(Expression.matrix.z)[column_order(heatmap)] 
  
  mat = mat[,sample_order]
  Expression.matrix.z = Expression.matrix.z[,sample_order]
  df = df[sample_order,]
}

# Heatmap annotation
ha = HeatmapAnnotation(`Mutational load` = anno_barplot(df$Mutational_load, baseline = 0),
                       `Mutation status` = df$Mutation_cat,
                       ICR = df$ICR,
                       CMS = df$CMS,
                       MSI = df$MSI,
                       col = list(`Mutation status` = c("nonhypermutated" = "#A7EABD", "hypermutated" = "#EAAED0"),
                                  ICR = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  CMS = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                          "mixed" = "lightgrey"),
                                  MSI = c("MSI-H" = "purple", "MSS" = "white"))
)


# Plotting
dir.create("./Figures/WES/009_Oncoprint", showWarnings = FALSE)
png(paste0("./Figures/WES/009_Oncoprint/009", version, subset, "_", ha_genes, 
           order_by_cluster_heatmap ,"_all_Expression_Heatmap_OncoPrint_ICR_low.png"), res = 600, width = 11, height = 14, #7.5,
    units = "in")
ht_list= oncoPrint(mat,
                   alter_fun = alter_fun, col = col, 
                   column_title = column_title, heatmap_legend_param = heatmap_legend_param,
                   top_annotation = ha,
                   row_order = rownames(mat),
                   column_order = sample_order) %v%
  Heatmap(Expression.matrix.z, cluster_rows = FALSE, cluster_columns = FALSE,
          show_column_names = FALSE, name = "Expression\n z score")
  #Heatmap(matrix(rnorm(nrow(mat)*10), ncol = 10), name = "expr", width = unit(4, "cm"))
draw(ht_list)  
dev.off()

