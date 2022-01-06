
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("dplyr", "tidyverse")
ipak(required.packages)

# Set parameters
version = "I"

# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")
DDR = read.csv("./Processed_Data/External/Gene_collections/IND99_GENES_ANNOTATION_JR_V5_DB.csv", stringsAsFactors = FALSE)
all_patients = read.csv(paste0("./Analysis/WES/003_ICRscore_MUT_vs_WT/all_patients_cancer_driver_QGP_table_MUT_WT_ICRscore.csv"), stringsAsFactors = FALSE)
hypermutated = read.csv(paste0("./Analysis/WES/003_ICRscore_MUT_vs_WT/hypermutated_cancer_driver_QGP_table_MUT_WT_ICRscore.csv"), stringsAsFactors = FALSE)
nonhypermutated = read.csv(paste0("./Analysis/WES/003_ICRscore_MUT_vs_WT/nonhypermutated_cancer_driver_QGP_table_MUT_WT_ICRscore.csv"), stringsAsFactors = FALSE)

# Filter minimum number of mutations
if(version %in% c("B", "H", "I")){
  all_patients = all_patients[which(all_patients$delta < 0),]
  nonhypermutated = nonhypermutated[which(nonhypermutated$delta < 0),]
  hypermutated = hypermutated[which(hypermutated$delta < 0),]
}
if(version %in% c("F", "G")){
  all_patients = all_patients[which(all_patients$delta > 0),]
  nonhypermutated = nonhypermutated[which(nonhypermutated$delta > 0),]
  hypermutated = hypermutated[which(hypermutated$delta > 0),]
}

all_patients = all_patients[which(all_patients$n_MUT >4),]
all_patients = all_patients[which(all_patients$p_value < 0.05),]

nonhypermutated = nonhypermutated[which(nonhypermutated$n_MUT >4),]
nonhypermutated = nonhypermutated[which(nonhypermutated$p_value < 0.05),]

hypermutated = hypermutated[which(hypermutated$n_MUT >4),]
hypermutated = hypermutated[which(hypermutated$p_value < 0.05),]


if(version == "B"){
  genes = c(all_patients$Gene, nonhypermutated$Gene, hypermutated$Gene)
}
if(version == "C"){
  genes = DDR$Gene[which(DDR$NEW.ANNOTATION_CONSENSUS_Sept_11 %in% c("BRCA1", "BRCA2", "FA", "HR-NON BRCA"))]
}
if(version == "D"){
  genes = c("BRCA1", "BRCA2", "FANCA")
}
if(version == "E"){
  genes = c("BRCA1", "BRCA2", "FANCA", "POLQ")
}
if(version == "F"){
  genes = c(nonhypermutated$Gene, hypermutated$Gene)
}
if(version == "G"){
  genes = c(all_patients$Gene)
}
if(version == "H"){
  genes = c(all_patients$Gene)
}
if(version == "I"){
  genes = c(nonhypermutated$Gene, hypermutated$Gene)
}

# Gene count ranking
finalMafFiltered$Tumor_Sample_Barcode = as.character(finalMafFiltered$Tumor_Sample_Barcode)
finalMafFiltered$Hugo_Symbol = as.character(finalMafFiltered$Hugo_Symbol)
MAF_gene_Count= finalMafFiltered %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(count=n()) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%     
  summarise(count=n()) %>% group_by(Hugo_Symbol) %>% 
  summarise(count=n()) %>% arrange(desc(count))

matrix = matrix(nrow = length(unique(finalMafFiltered$Sample_ID)), ncol = nrow(MAF_gene_Count))
colnames(matrix) = MAF_gene_Count$Hugo_Symbol
rownames(matrix) = unique(finalMafFiltered$Sample_ID)

matrix[which(is.na(matrix))] = ""
matrix = matrix[,genes]
#matrix = matrix[,c(148, 1:100)] # 132 = POLE

i = 1
for(i in 1:ncol(matrix)){
  gene = colnames(matrix)[i]
  sub = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  samples = sub$Sample_ID
  matrix[,gene][which(rownames(matrix) %in% samples)] = "MUT"
}

mat = t(matrix)

dir.create("./Analysis/WES/008_matrix_for_Oncoplot", showWarnings = FALSE)
save(mat, file = paste0("./Analysis/WES/008_matrix_for_Oncoplot/008", 
                        version, "_matrix_for_oncoplot.Rdata"))

