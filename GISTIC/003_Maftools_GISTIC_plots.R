
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ComplexHeatmap", "reshape2")
ipak(required.packages)
ibiopak("maftools")


# Load data
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_Seg_file.Rdata")

sidra_lumc_gistic = readGistic(gisticAllLesionsFile = "./Processed_Data/WES/GISTIC/all_lesions.conf_95.txt",
                               gisticAmpGenesFile = "./Processed_Data/WES/GISTIC/amp_genes.conf_95.txt",
                               gisticDelGenesFile = "./Processed_Data/WES/GISTIC/del_genes.conf_95.txt", 
                               gisticScoresFile = "./Processed_Data/WES/GISTIC/scores.gistic", isTCGA = FALSE)

dir.create("./Figures/WES/GISTIC", showWarnings = FALSE)
dir.create("./Figures/WES/GISTIC/003_Maftools_gistic_plots", showWarnings = FALSE)

png("./Figures/WES/GISTIC/003_Maftools_gistic_plots/gisticChromPlot.png",
    res = 600, width = 9, height = 5, units = "in")
gisticChromPlot(gistic = sidra_lumc_gistic, markBands = "all")
dev.off()

png("./Figures/WES/GISTIC/003_Maftools_gistic_plots/gisticBubblePlot.png",
    res = 600, width = 6, height = 6, units = "in")
gisticBubblePlot(gistic = sidra_lumc_gistic)
dev.off()

png("./Figures/WES/GISTIC/003_Maftools_gistic_plots/segment_plot.png",
    res = 600, width = 10, height = 4, units = "in")
plotCBSsegments(cbsFile = "./Processed_Data/WES/GISTIC/298SampleMerged.seg")
dev.off()
