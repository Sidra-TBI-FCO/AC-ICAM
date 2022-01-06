
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Load data
load("./Analysis/WES/GISTIC/001_Data_prep/del_amp_genes_vectors.Rdata")
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_thresholded_by_genes.Rdata")

Annotation = read.table("./Analysis/WES/GISTIC/001_Data_prep/Annotation_for_281_patients.txt", sep = "\t", header = TRUE)
load("../../TBI-LAB - General/Bioinformatics tools/GeneInfo/geneInfo.July2017.RData")

# Set parameters
q_value_filter = "" # "q_value_filter" or ""
shallow_included = "shallow_included" # "shallow_included" or "shallow_excluded"

# Data preparation

if(shallow_included == "shallow_included"){
  mat_amp = as.matrix(df_GISTIC)
  mat_amp[mat_amp == -1] = 0
  mat_amp[mat_amp == -2] = 0
  mat_amp[mat_amp == 1] = 1
  mat_amp[mat_amp == 2] = 1
  
  mat_del = as.matrix(df_GISTIC)
  mat_del[mat_del == 1] = 0
  mat_del[mat_del == 2] = 0
  mat_del[mat_del == -1] = 1
  mat_del[mat_del == -2] = 1
}

if(shallow_included == "shallow_excluded"){
  mat_amp = as.matrix(df_GISTIC)
  mat_amp[mat_amp == -1] = 0
  mat_amp[mat_amp == -2] = 0
  mat_amp[mat_amp == 1] = 0
  mat_amp[mat_amp == 2] = 1
  
  mat_del = as.matrix(df_GISTIC)
  mat_del[mat_del == 1] = 0
  mat_del[mat_del == 2] = 0
  mat_del[mat_del == -1] = 0
  mat_del[mat_del == -2] = 1
}


# Analysis for amplification

if(q_value_filter == "q_value_filter"){
  mat_amp = mat_amp[which(rownames(mat_amp) %in% amp_genes),]
}

#mat_amp = mat_amp[-which(rowSums(mat_amp) %in% c(0, 1)),]
dim(mat_amp)

Annotation$CNV = NA
Annotation = Annotation[-which(Annotation$ICR == "ICR Intermediate"),]
Annotation$ICR = factor(Annotation$ICR, levels = c("ICR Low", "ICR High"))

results = data.frame(Gene = rownames(mat_amp), p_val = NA, chi_sq_statistic = NA,
                     direction = NA)

i=2
for (i in 1:nrow(mat_amp)){
  gene = rownames(mat_amp)[i]
  Annotation$CNV = mat_amp[gene,][match(Annotation$Sample, colnames(mat_amp))]
  if(sum(Annotation$CNV) > 1){
    tbl = table(Annotation$ICR, Annotation$CNV)
    test = chisq.test(tbl)
    results$p_val[which(results$Gene == gene)] = test$p.value
    results$chi_sq_statistic[which(results$Gene == gene)] = test$statistic
    if(tbl["ICR High", "1"] > tbl["ICR Low", "1"]){
      results$direction[which(results$Gene == gene)] = "Amp_in_ICR_High"
    }
    if(tbl["ICR High", "1"] < tbl["ICR Low", "1"]){
      results$direction[which(results$Gene == gene)] = "Amp_in_ICR_Low"
    } 
  }else{
    results$p_val[which(results$Gene == gene)] = 1
  }
}
results$FDR = p.adjust(results$p_val, n = nrow(results), method = "fdr")
results = results[order(results$p_val),]

results$chr = geneInfo$chr[match(results$Gene, rownames(geneInfo))]
results$cytoband = geneInfo$band[match(results$Gene, rownames(geneInfo))]

results_sig = results[which(results$p_val < 0.05),]
table(results_sig$direction)

table(results_sig$chr)

dir.create("./Analysis/WES/GISTIC/004_Differentially_CNV_ICR_High_ICR_Low", showWarnings = FALSE)
write.csv(results, file = paste0("./Analysis/WES/GISTIC/004_Differentially_CNV_ICR_High_ICR_Low/Sign_amplified_genes_", q_value_filter, "_",
                                 shallow_included, ".csv"),
          row.names = FALSE)

results_amp = results

# Analysis for deletions

if(q_value_filter == "q_value_filter"){
  mat_del = mat_del[which(rownames(mat_del) %in% del_genes),]
}
#mat_del = mat_del[-which(rowSums(mat_del)  %in% c(0, 1)),]
dim(mat_del)

results = data.frame(Gene = rownames(mat_del), p_val = NA, chi_sq_statistic = NA,
                     direction = NA)

i=2
for (i in 1:nrow(mat_del)){
  gene = rownames(mat_del)[i]
  Annotation$CNV = mat_del[gene,][match(Annotation$Sample, colnames(mat_del))]
  if(sum(Annotation$CNV) > 1){
    tbl = table(Annotation$ICR, Annotation$CNV)
    test = chisq.test(tbl)
    results$p_val[which(results$Gene == gene)] = test$p.value
    results$chi_sq_statistic[which(results$Gene == gene)] = test$statistic
    if(tbl["ICR High", "1"] > tbl["ICR Low", "1"]){
      results$direction[which(results$Gene == gene)] = "Del_in_ICR_High"
    }
    if(tbl["ICR High", "1"] < tbl["ICR Low", "1"]){
      results$direction[which(results$Gene == gene)] = "Del_in_ICR_Low"
    } 
  }else{
    results$p_val[which(results$Gene == gene)] = 1
  }
}
results$FDR = p.adjust(results$p_val, n = nrow(results), method = "fdr")
results = results[order(results$p_val),]

results$chr = geneInfo$chr[match(results$Gene, rownames(geneInfo))]
results$cytoband = geneInfo$band[match(results$Gene, rownames(geneInfo))]

results_sig = results[which(results$p_val < 0.05),]
table(results_sig$direction)

results_del = results

dir.create("./Analysis/WES/GISTIC/004_Differentially_CNV_ICR_High_ICR_Low", showWarnings = FALSE)
write.csv(results, file = paste0("./Analysis/WES/GISTIC/004_Differentially_CNV_ICR_High_ICR_Low/Sign_deleted_genes_", q_value_filter, "_",
                                 shallow_included, ".csv"),
          row.names = FALSE)
