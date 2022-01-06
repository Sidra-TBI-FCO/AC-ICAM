
# Generate tertiles of the cohort based on expression of gene of interest

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/heatmap.3.R"))

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper", "heatmap3", "plyr", "spatstat",
                      "heatmap3")
ipak(required.packages)

# Set parameters
gene_of_interest = "CDH1"
mean_of = "" #c("HLA-A", "HLA-B", "HLA-C")
within_ICR_group = "" # "" if in all, "ICR High" if only ICR High group etc
name_mean_of = "" #"HLA class I" # For name in output file
if(mean_of > 1){
  outputname = name_mean_of
}else{
  outputname = gene_of_interest
}

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

if(mean_of > 1){
  calc_matrix = RNASeq.QN.LOG2[mean_of,]
  calc_df = as.data.frame(t(calc_matrix))
  calc_df$Expression = rowMeans(calc_df)
  df = data.frame(Sample = rownames(calc_df), Expression = calc_df$Expression, ICR = NA)
}else{
  df = data.frame(Sample = colnames(RNASeq.QN.LOG2), Expression = RNASeq.QN.LOG2[gene_of_interest,], ICR = NA)
}

df$ICR = table_cluster_assignment$ICR_HML[match(df$Sample, rownames(table_cluster_assignment))]
df$ICR = factor(df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))

if(within_ICR_group == ""){}else{
  df = df[which(df$ICR == within_ICR_group),]
}

df = df[-which(substring(rownames(df), 1, 3) %in% excluded_df$Patient_ID),]

df$Expression_category = NA
first = unname(quantile(df$Expression, probs = seq(0, 1, length.out = 4))[2])
second = unname(quantile(df$Expression, probs = seq(0, 1, length.out = 4))[3])
df$Expression_category[which(df$Expression < first)] = "Low"
df$Expression_category[which(df$Expression >= first & df$Expression < second)] = "Medium"
df$Expression_category[which(df$Expression >= second)] = "High"

df$Expression_category = as.character(df$Expression_category)

dir.create("./Analysis/Trimmed_p/Gene_of_interest_tertiles", showWarnings = FALSE)
save(df, file = paste0("./Analysis/Trimmed_p/Gene_of_interest_tertiles/", outputname, "_tertiles_within", within_ICR_group,".Rdata"))
