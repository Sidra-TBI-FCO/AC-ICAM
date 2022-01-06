

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_All_Immune_gene_signatures_table.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
labels = read.csv("./Processed_Data/External/Immune Traits Annotation/TableS2_Immune Traits_v2020921_REVISED.csv", stringsAsFactors = FALSE)

# Select right columns (signatures)
Thorsson = colnames(immune_sig_df)[grep("Thorsson", colnames(immune_sig_df))]
Benci = colnames(immune_sig_df)[grep("Benci", colnames(immune_sig_df))]
ConsensusTME = colnames(immune_sig_df)[grep("ConsensusTME", colnames(immune_sig_df))]
TIS = colnames(immune_sig_df)[grep("TIS", colnames(immune_sig_df))]
#ICR = colnames(immune_sig_df)[grep("ICR", colnames(immune_sig_df))][c(1,2)]
immune_sig_df = immune_sig_df[, c(Thorsson, Benci, ConsensusTME, TIS)] # 128 signatures remain

# Remove bindea signatures
Thorsson_bindea = colnames(immune_sig_df)[grep("Bindea", colnames(immune_sig_df))]
immune_sig_df = immune_sig_df[,-which(colnames(immune_sig_df) %in% Thorsson_bindea)] # 103 signatures remain

# Remove patients (from 366 -> 348 patients)
immune_sig_df = immune_sig_df[which(rownames(immune_sig_df) %in% colnames(RNASeq.QN.LOG2)),]

# Fix labels
labels = labels[which(labels$UseCase.RareVariants == 1),]
labels$Sayaman.InternalLabel = gsub("\\.", "\\_", labels$Sayaman.InternalLabel)
labels$Sayaman.InternalLabel = gsub("Sigs160_", "", labels$Sayaman.InternalLabel)

colnames(immune_sig_df) = gsub("Thorsson |  ", "", colnames(immune_sig_df))
colnames(immune_sig_df) = gsub("\\|", "", colnames(immune_sig_df))
colnames(immune_sig_df) = gsub("\\.", "\\_", colnames(immune_sig_df))
colnames(immune_sig_df) = gsub("ConsensusTME", "ConsensusTME -", colnames(immune_sig_df))
colnames(immune_sig_df) = gsub("Benci", "Expression Signature - Benci", colnames(immune_sig_df))
colnames(immune_sig_df) = gsub("ICR_ICR", "Expression Signature - ICR", colnames(immune_sig_df))
colnames(immune_sig_df) = gsub("TIS TIS_genes", "Expression Signature - TIS", colnames(immune_sig_df))

Friendly = labels$Sayaman.FriendlyLabel[match(colnames(immune_sig_df),labels$Sayaman.InternalLabel)]

colnames(immune_sig_df) = gsub("_", " ", colnames(immune_sig_df))
Friendly[which(is.na(Friendly))] = colnames(immune_sig_df)[which(is.na(Friendly))]
Friendly[which(Friendly == "Attractors G HLA-DPA1")] = "Attractor Metagene - G HLA DPA1"

colnames(immune_sig_df) = Friendly

save(immune_sig_df, file = "./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
