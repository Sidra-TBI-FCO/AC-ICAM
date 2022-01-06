
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full"    # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
min_non_zero_samples = 0.10
subset = ""
signatures_inclusion = "" # "restricted_subset_signatures"

# Load data
annotated = read.csv(paste0("./Analysis/Microbiome/017.2_Annotate_results/April_2021_", subset, min_non_zero_samples, "_v5_Annotated_Spearman_results.csv"),
                     stringsAsFactors = FALSE)
annotated = annotated[-which(is.na(annotated$rho)),]

# Genus to include based on min non zero samples
genus_to_include = annotated$Genus

load(paste0("./Analysis/Microbiome/017_Spearman_continous/", Type, "/", Rank, "/v5_April_2021_",
            "_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank,".Rdata"))
complete_results = results_all

complete_results = complete_results[which(complete_results$Name %in% genus_to_include),]
rownames(complete_results) = paste(complete_results$Name, complete_results$signature, sep = "_")

load(paste0("./Analysis/Microbiome/017_Spearman_continous/", Type, "/", Rank, "/v5_April_2021_",
            "nonhypermutated_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank,".Rdata"))
nonhypermutated = results_all
nonhypermutated = nonhypermutated[which(nonhypermutated$Name %in% genus_to_include),]
rownames(nonhypermutated) = paste(nonhypermutated$Name, nonhypermutated$signature, sep = "_")

load(paste0("./Analysis/Microbiome/017_Spearman_continous/", Type, "/", Rank, "/v5_April_2021_",
            "hypermutated_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank,".Rdata"))
hypermutated = results_all
hypermutated = hypermutated[which(hypermutated$Name %in% genus_to_include),]
rownames(hypermutated) = paste(hypermutated$Name, hypermutated$signature, sep = "_")


load(paste0("./Analysis/Microbiome/017_Spearman_continous/", Type, "/", Rank, "/v5_April_2021_",
            "Left sided_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank,".Rdata"))
left = results_all
left = left[which(left$Name %in% genus_to_include),]
rownames(left) = paste(left$Name, left$signature, sep = "_")

load(paste0("./Analysis/Microbiome/017_Spearman_continous/", Type, "/", Rank, "/v5_April_2021_",
            "Right sided_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank,".Rdata"))
right = results_all
right = right[which(right$Name %in% genus_to_include),]
rownames(right) = paste(right$Name, right$signature, sep = "_")

complete_results = merge(complete_results, nonhypermutated, by = "row.names")
rownames(complete_results) = complete_results$Row.names
complete_results$Row.names = NULL

complete_results = merge(complete_results, hypermutated, by = "row.names")
rownames(complete_results) = complete_results$Row.names
complete_results$Row.names = NULL

complete_results$Name.y = NULL
complete_results$signature.y = NULL
complete_results$Name = NULL
complete_results$signature = NULL

colnames(complete_results)
colnames(complete_results) = c(Rank, "Signature", "p value overall", "rho overall", 
                               "p value nonhypermutated", "rho nonhypermutated",
                               "p value hypermutated", "rho hypermutated")

complete_results = merge(complete_results, right, by = "row.names")
rownames(complete_results) = complete_results$Row.names
complete_results$Row.names = NULL

complete_results = merge(complete_results, left, by = "row.names")
rownames(complete_results) = complete_results$Row.names
complete_results$Row.names = NULL

complete_results$Name.y = NULL
complete_results$signature.y = NULL
complete_results$Name.x = NULL
complete_results$signature.x = NULL

colnames(complete_results)
colnames(complete_results) = c(Rank, "Signature", "p value overall", "rho overall", 
                               "p value nonhypermutated", "rho nonhypermutated",
                               "p value hypermutated", "rho hypermutated",
                               "p value right sided", "rho right sided",
                               "p value left sided", "rho left sided")

#complete_results = complete_results[-which(is.na(complete_results$`p value overall`)),] # Remove p-value NA (due to complete absence of the bacterial rank)
complete_results = complete_results[order(complete_results$`p value overall`),]

if(signatures_inclusion == "restricted_subset_signatures"){
  complete_results = complete_results[c(grep("ConsensusTME", complete_results$Signature), grep("ICR_SCORE", complete_results$Signature),
                                       grep("ImmunoSeq", complete_results$Signature), grep("MiXCR", complete_results$Signature),
                                       grep("Mutation", complete_results$Signature), grep("PDL1_data", complete_results$Signature),
                                       grep("Wolf.MHC.I_19272155", complete_results$Signature), grep("Wolf.MHC.II_19272155", complete_results$Signature),
                                       grep("CSF1_response", complete_results$Signature), grep("CHANG_CORE_SERUM_RESPONSE_UP", complete_results$Signature),
                                       grep("TGFB_score_21050467", complete_results$Signature), grep("Module3_IFN_score", complete_results$Signature),
                                       grep("LIexpression_score", complete_results$Signature)),]
}

dir.create(paste0("./Analysis/Microbiome/018_Spearman_combine"), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/018_Spearman_combine/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/018_Spearman_combine/",  Type, "/", Rank), showWarnings = FALSE)

save(complete_results, file = paste0("./Analysis/Microbiome/018_Spearman_combine/",  Type, "/", Rank,
                                "/April_2021_", signatures_inclusion, min_non_zero_samples, "Combined_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank, ".Rdata"))
write.csv(complete_results, file = paste0("./Analysis/Microbiome/018_Spearman_combine/",  Type, "/", Rank,
                                     "/April_2021_", signatures_inclusion, min_non_zero_samples, "Combined_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank, ".csv"),
          row.names = FALSE)

