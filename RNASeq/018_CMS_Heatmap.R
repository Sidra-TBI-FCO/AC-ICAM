
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper", "heatmap3", "plyr", "spatstat",
                      "heatmap3")
library(devtools)
install_github("Sage-Bionetworks/CMSclassifier")

source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")

library(synapseClient)
library(CMSclassifier)
CMS_genes_entrez = CMSclassifier::listModelGenes("RF")
load("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/GeneInfo/geneInfo.July2017.RData")
CMS_gene_Info = geneInfo[which(geneInfo$EntrezID %in% CMS_genes_entrez), ]
CMS_genes = row.names(CMS_gene_Info)

ipak(required.packages)

# Set Parameters

# Load data
load("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/ICR genes/ICR_genes.RData")

# Create folders and log file
dir.create("./Figures",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/018_CMS_Heatmap", showWarnings = FALSE)

load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

CMS_subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% CMS_genes, ])

Rfcms[is.na(Rfcms$RF.predictedCMS),6] = "Not_predicted"
Rfcms = Rfcms[order(Rfcms$RF.predictedCMS),]
Rfcms$ICR = table_cluster_assignment$ICR_HML[match(row.names(Rfcms), row.names(table_cluster_assignment))]
Rfcms$Sample_ID = row.names(Rfcms)

annotation = Rfcms[,c("ICR", "RF.predictedCMS")]

annotation$RF.predictedCMS = as.factor(annotation$RF.predictedCMS)

annotation.blot = as.matrix(annotation)
annotation.blot[annotation.blot=="CMS1"]="#FF9F21"
annotation.blot[annotation.blot=="CMS2"]="#0074AF"
annotation.blot[annotation.blot=="CMS3"]="#E97AA8"
annotation.blot[annotation.blot=="CMS4"]="#009E74"
annotation.blot[annotation.blot=="Not_predicted"]="lightgrey"

annotation.blot[annotation.blot=="ICR High"]="red"
annotation.blot[annotation.blot=="ICR Medium"]="green"
annotation.blot[annotation.blot=="ICR Low"]="blue"

my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
expression.matrix = t(CMS_subset_RNAseq_log2)[,rownames(annotation.blot)]

png(paste0("./Figures/Trimmed_p/018_CMS_Heatmap/CMS_Heatmap1.png"),res=600,height=12,width=10,unit="in")
heatmap.3 (expression.matrix,
           col=my.palette,
           Colv=NA,
           scale = "row",
           side.height.fraction = 0.3,
           ColSideColors= annotation.blot,
           ColSideLabs = c( paste0("ICR Low/Medium/High: ", gsub(",","/",toString(count(annotation, "ICR")$freq))),
                            paste0("CMS1/ 2/ 3/ 4/ NP: ", gsub(",","/",toString(count(annotation, "RF.predictedCMS")$freq)))),
           font_size_col_Labs = 1.3,
           cexCol = 0.9,
           labCol=FALSE,
           cexRow = 0.36,
           margins=c(15,15))

title(main = list(paste0("CMS Heatmap Sidra-LUMC cohort"),cex = 1.5),
      outer = FALSE, line = -4 , adj = 0.55)
legend("topright", legend = c("CMS1", "CMS2", "CMS3", "CMS4", "Not Predicted", "ICR Low", "ICR Medium", "ICR High"),
       col = c("#FF9F21", "#0074AF","#E97AA8", "#009E74", "lightblue", "blue", "green", "red"), lty= 0.5,lwd = 0.5, cex = 0.9, pch= 15, pt.cex = 1)
dev.off()
