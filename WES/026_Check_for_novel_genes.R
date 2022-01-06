
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters


# Load data
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")
Sidra_LUMC = read.csv("./Analysis/WES/011_Compare_with_TCGA/Oct_2021_Sidra_LUMC_frequency_compared_with_TCGA_NHS_PFHS.csv", stringsAsFactors = FALSE)

# Check
MUT_above_5_percent = Sidra_LUMC$Hugo.Symbol[which(Sidra_LUMC$Percentage.in.Sidra.LUMC.cohort > 5)]
overlap_with_QGPC = MUT_above_5_percent[which(MUT_above_5_percent %in% QGPC_genes)]
not_in_QGPC = MUT_above_5_percent[-which(MUT_above_5_percent %in% QGPC_genes)]

potentially_novel_df = Sidra_LUMC[which(Sidra_LUMC$Hugo.Symbol %in% not_in_QGPC),]

# remove all that a also above 5% in TCGA-COAD or NHS-HPFS
potentially_novel_df = potentially_novel_df[-which(potentially_novel_df$Percentage.in.TCGA.COAD > 5),]
potentially_novel_df = potentially_novel_df[-which(potentially_novel_df$Giannakis..Percentage.in.NHS.HPFS.Colon.Cancer > 5),]

potentially_novel_df = potentially_novel_df[, c("Hugo.Symbol",
                                                "Percentage.in.Sidra.LUMC.cohort",
                                                "Percentage.in.TCGA.COAD",
                                                "Giannakis..Percentage.in.NHS.HPFS.Colon.Cancer",
                                                "Gene.in.Colaprico.COAD",
                                                "Gene.in.Bailey.COADREAD")]
