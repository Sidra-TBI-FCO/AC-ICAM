
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set parameters
exclude = c("Conpair_lower_90_percent")
Geneset = "Immune_traits" # Selected.pathways" # Selected Pathways clean-up in script: 022_Dotted_Heatmap_CMS_ICR
# Immune traits, see script 048 in RNASeq pipeline, correlation matrix Sayaman immune traits

# Load data
if(Geneset == "Selected.pathways"){
  load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/Final_Selected.pathways.Sidra.LUMC_ES.Rdata"))
  ES = ES_selected_pathways
}

if(Geneset == "Immune_traits"){
  load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
  ES = t(as.matrix(immune_sig_df))
}

excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Prepare data
rownames(ES)
ES = data.frame(ES)
colnames(ES) = gsub("X", "", colnames(ES))
colnames(ES)

colnames(ES) = paste("SER-SILU-CC-P0", colnames(ES), sep = "")
colnames(ES) = gsub("T_P", "-PT-01", colnames(ES))
colnames(ES) = gsub("LM1", "-LM-01", colnames(ES))
colnames(ES) = gsub("LM2", "-LM-02", colnames(ES))
colnames(ES) = gsub("T_B", "-PT-02", colnames(ES))
colnames(ES) = gsub("T_C", "-PT-03", colnames(ES))

rownames(ES) = gsub("\\[", "", rownames(ES))
rownames(ES) = gsub("\\]", "", rownames(ES))
rownames(ES) = gsub("\\/", "-", rownames(ES))
rownames(ES) = gsub("\\ ", "_", rownames(ES))
rownames(ES) = gsub("\\_-_", "-", rownames(ES))
rownames(ES) = toupper(rownames(ES))

ES$geneset_id = rownames(ES)
ES = ES[,c(ncol(ES), 1:(ncol(ES) -1))]
rownames(ES) = NULL

write.table(ES,
            file = paste0("./Processed_Data/Shared_Data/For cbioportal/SER-SILU-GSVA_", Geneset,"_for_cbioportal.tsv"),
            sep = "\t", row.names = FALSE)

ES[1:103,2:349] = 0.01

write.table(ES,
            file = paste0("./Processed_Data/Shared_Data/For cbioportal/SER-SILU-GSVA_", Geneset,"_pvals_for_cbioportal.tsv"),
            sep = "\t", row.names = FALSE)
