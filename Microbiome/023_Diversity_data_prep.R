
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "stringr"))

# Load data
load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/clean_paired_rank_level_abundancies.Rdata")
alpha_diversity_ICR = read.csv("./Processed_Data/Microbiome/23_Alpha_diversity/From Arun/Alpha_ICR-HiLo.csv",
                               stringsAsFactors = FALSE)

alpha_diversity_AJCC = read.csv("./Processed_Data/Microbiome/23_Alpha_diversity/From Arun/Alpha_AJCC-1234.csv",
                                stringsAsFactors = FALSE)

alpha_diversity_tumor_normal = read.csv("./Processed_Data/Microbiome/23_Alpha_diversity/From Arun/Alpha_Tissue-NT-paired.csv",
                                        stringsAsFactors = FALSE)

# Extract relevant data
alpha_diversity_tumor_normal$SampleCode = str_pad(alpha_diversity_tumor_normal$SampleCode, pad = "0", 3)
alpha_diversity_tumor_normal$Sample_ID = paste(alpha_diversity_tumor_normal$SampleCode,
                                                   alpha_diversity_tumor_normal$SamplesDescription,
                                                   sep = "")
  
alpha_diversity_tumor_normal$Patient_ID = alpha_diversity_tumor_normal$SampleCode
alpha_diversity_tumor_normal$Tissue = alpha_diversity_tumor_normal$SamplesDescription

df = alpha_diversity_tumor_normal

df_alpha = df[, c("Sample_ID", "Patient_ID", "Tissue", "Observed", "Chao1",
                            "Shannon", "InvSimpson")]

which(df_alpha$Patient_ID %in% c("334", "393"))

included_samples = colnames(Genus_full_abundance)

df_alpha = df_alpha[which(df_alpha$Sample_ID %in% included_samples),]

df_alpha = df_alpha[order(df_alpha$Sample_ID),]
rownames(df_alpha) = 1:nrow(df_alpha)

dir.create("./Processed_Data/Microbiome/23_Alpha_diversity/Sample_based_filted", showWarnings = FALSE)
save(df_alpha, file = "./Processed_Data/Microbiome/23_Alpha_diversity/Sample_based_filted/alpha_diversity_TN_pairs.Rdata")
