
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)

# Load data
load("./Data/RNASeq/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Data/RNASeq_TCGA_COAD/COAD_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
filtered.norm.RNAseqData = log(filtered.norm.RNAseqData +1, 2) 

genes_in_both = intersect(rownames(filtered.norm.RNAseqData), rownames(RNASeq.QN.LOG2))
filtered.norm.RNAseqData = filtered.norm.RNAseqData[genes_in_both,]
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[genes_in_both,]


## JSREP
write.table(RNASeq.QN.LOG2, sep = "\t", file = "./Data/RNASeq/Normalized/intersect_JSREP.clean.dataset.EDAseq.QN.HPC.txt", quote = FALSE)
filterCommonGenes(input.f=paste0("./Data/RNASeq/Normalized/intersect_JSREP.clean.dataset.EDAseq.QN.HPC.txt"),
                  output.f=paste0("./Analysis/ESTIMATE/intersect_JSREP.clean.ESTIMATE.input.gct"),
                  id=c("GeneSymbol","EntrezID"))

estimateScore(input.ds ="./Analysis/ESTIMATE/intersect_JSREP.clean.ESTIMATE.input.gct",
              output.ds = "./Analysis/ESTIMATE/intersect_JSREP.clean.ESTIMATE.score.gct",
              platform= "illumina")

estimate.gct<-read.table("./Analysis/ESTIMATE/intersect_JSREP.clean.ESTIMATE.score.gct", skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate.gct) = estimate.gct$NAME
estimate.gct$NAME = NULL
estimate.gct$Description = NULL

ESTIMATE = t(estimate.gct)
rownames(ESTIMATE) = gsub("X", "", rownames(ESTIMATE))
save(ESTIMATE, file = "./Analysis/ESTIMATE/intersect_JSREP_clean_ESTIMATE_scores.Rdata")


## TCGA

write.table(filtered.norm.RNAseqData, sep = "\t", file = "./Data/RNASeq_TCGA_COAD/intersect_COAD_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.txt", quote = FALSE)

filterCommonGenes(input.f=paste0("./Data/RNASeq_TCGA_COAD/intersect_COAD_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.txt"),
                  output.f=paste0("./Analysis/ESTIMATE/intersect_TCGA.COAD.whitelist.filtered.ESTIMATE.input.gct"),
                  id=c("GeneSymbol","EntrezID"))

estimateScore(input.ds ="./Analysis/ESTIMATE/intersect_TCGA.COAD.whitelist.filtered.ESTIMATE.input.gct",
              output.ds = "./Analysis/ESTIMATE/intersect_TCGA.COAD.whitelist.filtered.ESTIMATE.score.gct",
              platform= "illumina")

estimate_tcga.gct<-read.table("./Analysis/ESTIMATE/intersect_TCGA.COAD.whitelist.filtered.ESTIMATE.score.gct", skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate_tcga.gct) = estimate_tcga.gct$NAME
estimate_tcga.gct$NAME = NULL
estimate_tcga.gct$Description = NULL

ESTIMATE_tcga = t(estimate_tcga.gct)
save(ESTIMATE_tcga, file = "./Analysis/ESTIMATE/intersect_TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")


