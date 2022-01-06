
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "reshape2")
ipak(required.packages)

# Load data
aa_comparison = read.table("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/CDR3aaMiXCR.txt",
                           sep = "\t", header = TRUE)

# Quick check Gianni Monaco overview file
table(aa_comparison$SharedWithImmunoSeq)
aa_comparison_not_shared = aa_comparison[which(aa_comparison$SharedWithImmunoSeq == FALSE),]
aa_comparison_shared = aa_comparison[which(aa_comparison$SharedWithImmunoSeq == TRUE),]



# ImmunoSeq data prep
combined_rearrangements = read.table("./Processed_Data/TCR/CombinedRearrangements_v1.tsv", stringsAsFactors = FALSE,
                                     sep = "\t", header = TRUE)

df = combined_rearrangements
rownames(df) = df$Amino.Acid
df$Sum..Productive.Frequency. = NULL
df$Present.In = NULL
df$Amino.Acid = NULL
colnames(df) = gsub("\\X", "", colnames(df))
matrix = as.matrix(df)
matrix[matrix > 0] = 1  # make binary, detected yes or no (1/0)

ImmunoSeq_result = melt(matrix)

colnames(ImmunoSeq_result) = c("Amino acid", "Sample", "Detected")
ImmunoSeq_result = ImmunoSeq_result[-which(ImmunoSeq_result$Detected == 0),]
head(ImmunoSeq_result)

ImmunoSeq_result$Combination = paste(ImmunoSeq_result$`Amino acid`, ImmunoSeq_result$Sample, sep = "_")

ImmunoSeq_result$Tissue = substring(ImmunoSeq_result$Sample, 4, 6)
ImmunoSeq_result = ImmunoSeq_result[which(ImmunoSeq_result$Tissue  == "T"),] # only focus on the tumor

samples_ImmunoSeq = unique(as.character(ImmunoSeq_result$Sample))

# MiXCR data prep
MiXCR_aa = aa_comparison
rownames(MiXCR_aa) = MiXCR_aa$CD3aa
MiXCR_aa$SharedWithImmunoSeq = NULL
MiXCR_aa$SamplesInMiXCR = NULL
MiXCR_aa$SamplesInImmunoSeq = NULL
MiXCR_aa$PresentIn = NULL
MiXCR_aa$ProductiveFrequency = NULL
MiXCR_aa$CD3aa = NULL
colnames(MiXCR_aa) = gsub("\\X", "", colnames(MiXCR_aa))

matrix2 = as.matrix(MiXCR_aa)
matrix2[matrix2 > 0] = 1 # make binary, detected yes or no (1/0)

MiXCR_result = melt(matrix2)
head(MiXCR_result)

colnames(MiXCR_result) = c("Amino acid", "Sample", "Detected")
MiXCR_result = MiXCR_result[-which(MiXCR_result$Detected == 0),]
head(MiXCR_result)

MiXCR_result$Combination = paste(MiXCR_result$`Amino acid`, MiXCR_result$Sample, sep = "_")

MiXCR_result$Tissue = substring(MiXCR_result$Sample, 4, 6)
MiXCR_result = MiXCR_result[which(MiXCR_result$Tissue  == "T"),] # only focus on the tumor

samples_MiXCR = as.character(unique(MiXCR_result$Sample))

intersection = intersect(samples_ImmunoSeq, samples_MiXCR)

MiXCR_result = MiXCR_result[which(MiXCR_result$Sample %in% intersection),]
ImmunoSeq_result = ImmunoSeq_result[which(ImmunoSeq_result$Sample %in% intersection),]

# Shared or not
MiXCR_result$SharedWithImmunoSeq = "unique to MiXCR"
MiXCR_result$SharedWithImmunoSeq[which(MiXCR_result$Combination %in% ImmunoSeq_result$Combination)] = "shared"
table(MiXCR_result$SharedWithImmunoSeq, exclude = NULL)

ImmunoSeq_result$SharedWithMiXCR = "unique to ImmunoSeq"
ImmunoSeq_result$SharedWithMiXCR[which(ImmunoSeq_result$Combination %in% MiXCR_result$Combination)] = "shared"
table(ImmunoSeq_result$SharedWithMiXCR)

# 76.6% of CDR3 detected by MiXCR were also detected in the same samples by ImmunoSeq sequencing
length(which(MiXCR_result$SharedWithImmunoSeq == "shared")) / nrow(MiXCR_result)

# 0.22% of CDR3 detected by ImmunoSeq were also detected in the same samples by MiXCR sequencing
length(which(ImmunoSeq_result$SharedWithMiXCR == "shared")) / nrow(ImmunoSeq_result)
