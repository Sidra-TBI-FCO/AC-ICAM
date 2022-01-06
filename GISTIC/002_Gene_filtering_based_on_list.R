
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Load data
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_thresholded_by_genes.Rdata")
load("./Analysis/WES/GISTIC/001_Data_prep/del_amp_genes_vectors.Rdata")
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")
load("./Processed_Data/External/Gene_collections/501_genes_Michele.Rdata")
Immunegene_list = read.csv("./Processed_Data/External/Gene_collections/extensive Immunogene list_DB.csv", stringsAsFactors = FALSE)

immune_genes = Immunegene_list$Gene.Symbol

# Set parameters
select = "immune_genes" # "genes_Michele" "QGPC_genes" "immune_genes"

# Formatting
genes_of_interest = get(select)

# Filtering
select_amp_genes = amp_genes[which(amp_genes %in% genes_of_interest)]
select_del_genes = del_genes[which(del_genes %in% genes_of_interest)]

df_amp = df_GISTIC[select_amp_genes,]
df_del = df_GISTIC[select_del_genes,]

mat_amp = as.matrix(df_amp)
mat_amp[mat_amp == -1] = 0
mat_amp[mat_amp == 1] = 0
mat_amp[mat_amp == -2] = 0
mat_amp[mat_amp == 2] = 1

df = data.frame(rowSums(mat_amp))
df$percentage_amp = round(df$rowSums.mat_amp. / 281 *100, 2)
colnames(df)[1] = "Amplification in N samples"
df = df[order(df$`Amplification in N samples`, decreasing = TRUE),]

dir.create("./Analysis/WES/GISTIC/002_Gene_filtering_based_on_list", showWarnings = FALSE)
write.csv(df, file = paste0("./Analysis/WES/GISTIC/002_Gene_filtering_based_on_list/amp_genes_", select, ".csv"), 
          row.names = TRUE)

mat_del = as.matrix(df_del)
mat_del[mat_del == -1] = 0
mat_del[mat_del == 1] = 0
mat_del[mat_del == -2] = 1
mat_del[mat_del == 2] = 0

df = data.frame(rowSums(mat_del))
df$percentage_del = round(df$rowSums.mat_del. / 281 *100, 2)
colnames(df)[1] = "Deletion in N samples"
df = df[order(df$`Deletion in N samples`, decreasing = TRUE),]

dir.create("./Analysis/WES/GISTIC/002_Gene_filtering_based_on_list", showWarnings = FALSE)
write.csv(df, file = paste0("./Analysis/WES/GISTIC/002_Gene_filtering_based_on_list/del_genes_", select, ".csv"), 
          row.names = TRUE)

