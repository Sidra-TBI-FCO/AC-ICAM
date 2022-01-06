
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("RColorBrewer", "forestplot", "stringr", "ggplot2")
ipak(required.packages)

# Set parameters
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = "" 
included_stages = "All"
Surv.cutoff.years = 20
Survival_outcome = "DFS_Def1" # "OS" "DFS_Def1"
exclude = c("Conpair_lower_90_percent", "non-epithelial")
Source = "v2_All"

# Manual adjustment
order = NA

# Load data
load(paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_", Source,"_", Survival_outcome,"_HR_table_exclude_",
            str_c(exclude, collapse = "_"), 
            "_stages_", str_c(included_stages, collapse = "_"),
            exclude_adjuvant, exclude_stage_IV, "_cutoff_", Surv.cutoff.years,".Rdata"))
load("./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/048_Correlation_plot_order_pathways.Rdata")
labels = read.csv("./Processed_Data/External/Immune Traits Annotation/TableS2_Immune Traits_v2020921_REVISED.csv", stringsAsFactors = FALSE)

# Cleanup
labels = labels[which(labels$UseCase.RareVariants == 1),]
labels$Sayaman.InternalLabel = gsub("\\.", "\\_", labels$Sayaman.InternalLabel)
labels$Sayaman.InternalLabel = gsub("Sigs160_", "", labels$Sayaman.InternalLabel)

HR_table$Signature = as.character(HR_table$Signature)

# Select right columns (signatures)
Thorsson = HR_table$Signature[grep("Thorsson", HR_table$Signature)]
Benci = HR_table$Signature[grep("Benci", HR_table$Signature)]
ConsensusTME = HR_table$Signature[grep("ConsensusTME", HR_table$Signature)]
TIS = HR_table$Signature[grep("TIS", HR_table$Signature)]
HR_table = HR_table[which(HR_table$Signature %in% c(Thorsson, Benci, ConsensusTME, TIS)),] # 128 signatures remain

# Remove bindea signatures
Thorsson_bindea = HR_table$Signature[grep("Bindea", HR_table$Signature)]
HR_table = HR_table[-which(HR_table$Signature %in% Thorsson_bindea),] # 103 signatures remain

# Friendly labels
HR_table$Signature = gsub("Thorsson |  ", "", HR_table$Signature)
HR_table$Signature = gsub("\\|", "", HR_table$Signature)
HR_table$Signature = gsub("\\.", "\\_", HR_table$Signature)
HR_table$Signature = gsub("ConsensusTME", "Leukocyte Subset ES (ConsensusTME) -", HR_table$Signature)
HR_table$Signature = gsub("Benci", "Expression Signature - Benci", HR_table$Signature)
HR_table$Signature = gsub("ICR_ICR", "Expression Signature - ICR", HR_table$Signature)
HR_table$Signature = gsub("TIS TIS_genes", "Expression Signature - TIS", HR_table$Signature)

Friendly = labels$Sayaman.FriendlyLabel[match(HR_table$Signature,labels$Sayaman.InternalLabel)]

HR_table$Signature = gsub("_", " ", HR_table$Signature)

Friendly[which(is.na(Friendly))] = HR_table$Signature[which(is.na(Friendly))]
Friendly[which(Friendly == "Attractors G HLA-DPA1")] = "Attractor Metagene - G HLA DPA1"
Friendly = gsub(".*\\- ", "", Friendly)
Friendly
HR_table$Signature = Friendly

# Re-order HR-table for heatmap
rownames(HR_table) = HR_table$Signature

HR_table = HR_table[order,]

dir.create("./Figures/Trimmed_p/014.2_HR_Heatmap", showWarnings = FALSE)
write.csv(HR_table, file = "./Figures/Trimmed_p/014.2_HR_Heatmap/HR_Heatmap_DFS_Def1.csv")
