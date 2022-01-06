
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/As.data.frame.list.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "mygene", "reshape2", "org.Hs.eg.db")
ipak(required.packages)

# Load data
load("~/Dropbox (SMCP)/DB-LAB/bioinformatics tools/GSEA list/OLD/immune.gene.lists.v3.Rdata")
load("~/Dropbox (SMCP)/DB-LAB/bioinformatics tools/Selected Pathways/Selected.pathways.3.4.RData")

# Limma package try alias2Symbol.

# Set parameters
Gene.set_name = "Bindea_ORIG"

# Translate to entrez id
Gene.set = get(Gene.set_name)
DF = melt(Gene.set)
colnames(DF) = c("HugoSymbol", "Pathway")

gene_query = queryMany(DF$HugoSymbol, scopes = "symbol", fields=c("ensembl.gene", "entrezgene"), species = "human")
gene_table = as.data.frame(gene_query)
gene_table = gene_table[-which(is.na(gene_table$entrezgene)),]
DF$Entrez_ID = gene_table$entrezgene[match(DF$HugoSymbol, gene_table$query)]
DF = DF[-which(is.na(DF$Entrez_ID)),]

# Covert to gmt format
pathways <- unique(DF$Pathway)
Full_list <-list()
i=1
for (i in 1:length(pathways)){
  pathway = pathways[i]
  list.item = list(DF$Entrez_ID[DF$Pathway==pathway])
  names(list.item) <- pathway
  Full_list <- c(Full_list, list.item)
}

DF_new = as.data.frame.list(Full_list)
DF_new$Name = rownames(DF_new)
if(Gene.set_name == "Bindea_ORIG"){
  DF_new$Link = "https://doi.org/10.1016/j.immuni.2013.10.003"
}
if(Gene.set_name == "Selected.pathways"){
  DF_new$Link[grep("IPA", DF_new$Name)] = "IPA"
  DF_new$Link = "to_be_filled"
}

rownames(DF_new) = NULL
DF_new = DF_new[,c((ncol(DF_new) -1), ncol(DF_new) , 1:(ncol(DF_new) - 2))]

DF_new$Name = gsub("\\[", "", DF_new$Name)
DF_new$Name = gsub("\\]", "", DF_new$Name)
DF_new$Name = gsub("\\/", "-", DF_new$Name)
DF_new$Name = gsub("\\ ", "_", DF_new$Name)
DF_new$Name[which(DF_new$Name == "TPW_Immunogenic_Cell_Death_(ICD)")] = "TPW_Immunogenic_Cell_Death"
DF_new$Name = toupper(DF_new$Name)

matrix = as.matrix(DF_new)
matrix[is.na(matrix)] = ""

write.table(matrix, file = paste0("./Processed_Data/Shared_Data/For cbioportal/", Gene.set_name, ".gmt"), sep = "\t", 
            col.names = FALSE, row.names = FALSE)


