
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Load data
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_Seg_file.Rdata")
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_thresholded_by_genes.Rdata")

# 

colnames(df_GISTIC) = paste("SER-SILU-CC-P0", colnames(df_GISTIC), sep = "")
colnames(df_GISTIC) = gsub("T", "-PT-01-B-02", colnames(df_GISTIC))

df_GISTIC$Hugo_Symbol = rownames(df_GISTIC)
df_GISTIC = df_GISTIC[, c(ncol(df_GISTIC), 1:(ncol(df_GISTIC) - 1))]
rownames(df_GISTIC) = NULL

colnames(df_GISTIC) = substring(colnames(df_GISTIC), 1, 23)
write.table(df_GISTIC, row.names = FALSE, file = "./Processed_Data/Shared_Data/For cbioportal/data_CNA_SILU.txt", 
            sep= "\t", quote = FALSE)

seg_file$Sample = paste("SER-SILU-CC-P0", seg_file$Sample, sep = "")
seg_file$Sample = gsub("T", "-PT-01-B-02", seg_file$Sample)
seg_file$Segment_Mean = log(seg_file$Segment_Mean, 2)

seg_file$Sample = substring(seg_file$Sample, 1, 23)

colnames(seg_file) = c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

write.table(seg_file, row.names = FALSE, file = "./Processed_Data/Shared_Data/For cbioportal/data_cna_hg19_SILU.seg", 
            sep= "\t", quote = FALSE)
