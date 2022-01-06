
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

dir.create("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies",
           showWarnings = FALSE)

# Load data
load("./Processed_Data/Microbiome/NormObjects.RData")
load("./Processed_Data/Microbiome/001_data_preparation/translation_table_Sample_Codes.Rdata")

# Calculate relative abundance
Microbiome_Relative = transform_sample_counts(s16sV1V3.2_sil, function(x) x / sum(x) )
OTU_table = otu_table(Microbiome_Relative)
OTU_table_dat = OTU_table@.Data
dim(OTU_table_dat)
# 16514   570
colSums(OTU_table_dat)

# filter such that only OTUs with a mean greater than 10^-5 are kept.
#Microbiome_Relative_Filtr = filter_taxa(Microbiome_Relative, function(x) mean(x) > 1e-5, TRUE)
#OTU_table = otu_table(Microbiome_Relative_Filtr)
#OTU_table_dat = OTU_table@.Data
#dim(OTU_table_dat)
# 3756   570

colnames(OTU_table_dat) = translation_table$SampleCode_new[match(colnames(OTU_table_dat), translation_table$Sample)]
colSums(OTU_table_dat)

OTU_table_dat = OTU_table_dat[, order(colnames(OTU_table_dat))]
save(OTU_table_dat, file = "./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/s16sV1V3.2_sil_OTU_table.Rdata")
