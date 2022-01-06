
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Load data
load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Normalized_abundancies/s16sV1V3.2_sil_OTU_table.Rdata")
load("./Processed_Data/Microbiome/001_data_preparation/translation_table_Sample_Codes.Rdata")
load("./Overview Sample Stocks/Meta_data/Excluded_patients.Rdata")

# Subsetting data
translation_table$Tissue_type = substring(translation_table$SampleCode_new, 4, 6)
translation_table$Patient_ID = substring(translation_table$SampleCode_new, 1, 3)
translation_table = translation_table[which(translation_table$Tissue_type %in% c("T", "N")),] # Drop the liver meta samples, these will be used in a seperate analysis

profiled_samples = translation_table$SampleCode_new
length(unique(profiled_samples))
# 554 # without liver metastasis

patients_with_T = translation_table$Patient_ID[which(translation_table$Tissue_type == "T")]
patients_with_N = translation_table$Patient_ID[which(translation_table$Tissue_type == "N")]

patients_with_both = intersect(patients_with_T, patients_with_N)
# 250 patients with both

patients_with_both = patients_with_both[-which(patients_with_both %in% c("334", "393"))]

# 334","393
patients_with_both = patients_with_both[-which(patients_with_both %in% exclude$Patient_ID[which(exclude$Reason.excluded == "non-epithelial")])]


# Pairs
sample_pairs = translation_table$SampleCode_new[which(translation_table$Patient_ID %in% patients_with_both)]
OTU_table_dat = OTU_table_dat[, sample_pairs]
dim(OTU_table_dat)

save(OTU_table_dat, tax_table_dat, file = "./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Normalized_abundancies/clean_paired_V3.2_sil_OTU_table.Rdata")

load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/s16sV1V3.2_sil_OTU_table.Rdata")
dim(OTU_table_dat)
OTU_table_dat = OTU_table_dat[, sample_pairs]
dim(OTU_table_dat)

save(OTU_table_dat, tax_table_dat, file = "./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/clean_paired_V3.2_sil_OTU_table.Rdata")
