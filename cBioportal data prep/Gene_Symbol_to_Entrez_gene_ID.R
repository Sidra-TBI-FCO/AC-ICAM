
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak("mygene")

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.Complete.dataset.EDAseq.QN.HPC.Rdata")

# Prepare dataframe
geneInfo$hgnc_symbol = rownames(geneInfo)
gene.query <- queryMany(geneInfo$hgnc_symbol, scopes = "symbol", fields=c("ensembl.gene", "entrezgene"), species = "human")
gene.table <- as.data.frame(gene.query)

results = data.frame(Gene_symbol = rownames(RNASeq.QN.LOG2), Entrez_ID = NA)
results$Entrez_ID = gene.table$entrezgene[match(results$Gene_symbol, gene.table$query)]

write.csv(results, file = "./Processed_Data/Shared_Data/For cbioportal/SER-SILU-RNASeq-GeneSymbol_to_EntrezID.csv", row.names = FALSE)

