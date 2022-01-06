#################################################################
###
### This script creates a correlation matrix for the selected genes
###
### Data input:
### 
### Output :
### ICR_Correlation_plot.png"
#################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
colpattern = colorRampPalette(c("#FBFBFB", "#25AD22"))(100)
selected_genes = "ICR"                                                                                                   # Specify which genes will be correlated
genes_to_correlate = selected_genes
Log_file = paste0("./Logfiles/Correlation/ICR_correlation_plot_",gsub(":",".",gsub(" ","_",date())),".txt")
test = "pearson"
name_gene_collection = "ICR genes"

# Define parameters
if (selected_genes == "ICR") {
  load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))
  genes_to_correlate = ICR_genes
  name_gene_collection = "ICR genes"
}


# Create folders
dir.create("./Figures/Trimmed_p",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Correlation",showWarnings = FALSE)
dir.create("./Logfiles/Correlation",showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p/Correlation",showWarnings = FALSE)

cat("This is a log file for creating correlation plots",                                       # Set-up logfile
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("Genes to correlate = ", genes_to_correlate),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Scripts output :",
    "",
    "Creating heatmaps",
    file = Log_file,
    append = FALSE, sep= "\n")

# Make correlation plots
start.time <- Sys.time ()

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
unavailable.genes <- genes_to_correlate[-which(genes_to_correlate %in% row.names(RNASeq.QN.LOG2))]
subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% genes_to_correlate,])

# Correlation matrix
Selected_genes_cor <- cor (subset_RNAseq_log2,method=test)

# Correlation significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = test, conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
Selected_genes_cor_sign <- cor.mtest(subset_RNAseq_log2, 0.95)

png(paste0("./Figures/Trimmed_p/Correlation/", selected_genes, "_Correlation_plot.png"),res=600,height=11,width=11,unit="in")
dcex.before <- par("cex")
par(cex = 0.45)
lims=c(-1,1)
if (length(Selected_genes_cor[Selected_genes_cor<0]) == 0) {lims=c(0,1)}
annotation = data.frame (gene = rownames(Selected_genes_cor),color = c(rep("#01B050",20)),stringsAsFactors = FALSE)
annotation$color[annotation$gene %in% c("IDO1","CD274","CTLA4","FOXP3","PDCD1")] = "#CC0506"
annotation = annotation[corrMatOrder(Selected_genes_cor,order="FPC"),]

mean_correlation = round(mean(Selected_genes_cor),2)

corrplot.mixed (Selected_genes_cor,
                #type="lower",
                #p.mat = Selected_genes_cor_sign[[1]],                                                                      # add significance to correlations
                #col = colpattern,
                lower = "square",
                upper ="number",
                order="FPC",
                cl.lim=lims,                                                                                               # only positive correlations
                tl.pos ="lt",
                tl.col = as.character(annotation$color),
                insig= "pch",                                                                                              # remove insignificant correlations
                pch = "x",
                pch.cex= 1.5,
                tl.cex = 1.4/par("cex"),
                cl.cex = 1/par("cex"),
                cex.main = 1/par("cex"),
                mar=c(6,4.1,7,5))
title(main = list(paste0(name_gene_collection, " correlation plot"),
                  cex = 4), line = -2.5)
title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data is \n from JSREP project (LUMC and Sidra). \n",
                        "Mean correlation is: ", mean_correlation, " \n",
                        "Significance level of correlation is represented by the size of the squares."), cex = 3), line = 1.5)
#par(cex = cex.before)
dev.off()

cat(paste0("Mean correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")

save(Selected_genes_cor, Selected_genes_cor_sign, file = paste0("./Analysis/Trimmed_p/Correlation/", selected_genes, "_correlation_matrix_and_sig", ".Rdata"))
