
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "dplyr", "plyr", "ggpubr", "plotly"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full"    # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
FDR_cutoff = 0.001 # step one to get the relevant interactions included in the Sankey
FDR_cutoff_step2 = 0.05 # step two for plotting between included genus and signatures
include_ICR = ""
signatures_inclusion = "" # "restricted_subset_signatures"
min_non_zero_samples = 0.10
colors_by_rho = ""


colors = c("Expression Signature - GRANS PCA 16704732" = "rgba(31, 119, 180,0.25)", 
           "Expression Signature - PDL1 data" = "rgba(255, 127, 14,0.25)", 
           "Expression Signature - IL8 21978456" = "rgba(44, 160, 44,0.25)", 
           "Attractor Metagene - G SIGLEC9" = "rgba(214, 39, 40,0.25)", 
           "Expression Signature - TREM1 data" = "rgba(148, 103, 189,0.25)", 
           "Expression Signature - Chemokine12 score"  = "rgba(140, 86, 75, 0.25)", 
           "Leukocyte Subset ES - Macrophages M1" = "rgba(227, 119, 194, 0.25)", 
           "Leukocyte Subset ES - Immune Score" = "rgba(252, 186, 3,0.25)", 
           "Expression Signature - IFNG score 21050467" = "rgba(188, 189, 34,0.25)", 
           "Leukocyte Subset ES - Monocytes" = "rgba(23, 190, 207,0.25)",
           "Leukocyte Subset ES - Macrophages" = "rgba(77, 5, 232, 0.25)",
           "Leukocyte Subset ES - Mast cells" = "rgba(127, 127, 127,0.25)", 
           "Leukocyte Subset ES - Dendritic cells" = "rgba(106, 114, 184,0.25)",
           "Attractor Metagene - G LILRB4" = "rgba(184, 112, 135, 0.25)"
           )


# Load data
load(paste0("./Analysis/Microbiome/018_Spearman_combine/",  Type, "/", Rank,
            "/April_2021_", signatures_inclusion, min_non_zero_samples, "Combined_Spearman_by_all_Thorsson_ConsensusTME_Benci_", Type, "_", Rank, ".Rdata"))
labels = read.csv("./Processed_Data/External/Immune Traits Annotation/TableS2_Immune Traits_v2020921_REVISED.csv", stringsAsFactors = FALSE)
labels = labels[which(labels$UseCase.RareVariants == 1),]
labels$Sayaman.InternalLabel = gsub("\\.", "\\_", labels$Sayaman.InternalLabel)
labels$Sayaman.InternalLabel = gsub("Sigs160_", "", labels$Sayaman.InternalLabel)

complete_results$FDR = p.adjust(complete_results$`p value overall`, method = "fdr", n = nrow(complete_results))

if(include_ICR == "include_ICR"){
  extra = complete_results[which(complete_results$Signature == "Thorsson |  ICR.ICR_SCORE"),] 
  extra = extra[1:2,] # select first 2 rows
}

write.csv(complete_results, file = "./Analysis/Microbiome/019_Sankey_table_in_correct_order/Supplementary_Table_S5.csv", row.names = FALSE)

included_signatures = unique(complete_results$Signature[which(complete_results$FDR < FDR_cutoff)])
included_genus = unique(complete_results$Genus_full[which(complete_results$FDR < FDR_cutoff)])

complete_results = complete_results[which(complete_results$Genus_full %in% included_genus),]
complete_results = complete_results[which(complete_results$Signature %in% included_signatures),]

complete_results = complete_results[which(complete_results$FDR < FDR_cutoff_step2),]

if(include_ICR == "include_ICR"){
  complete_results = rbind(complete_results, extra)
}

complete_results$Phylum = gsub(".*\\D_1__", "", complete_results[, Rank])
complete_results$Phylum = gsub("\\ D_2_.*", "", complete_results$Phylum)
complete_results$Genus_full = gsub(".*\\D_5__", "", complete_results$Genus_full)
complete_results$Genus_full = gsub("\\[", "", complete_results$Genus_full)
complete_results$Genus_full = gsub("\\]", "", complete_results$Genus_full)
complete_results$Genus_full = paste0(complete_results$Genus_full, " (", complete_results$Phylum, ")")

complete_results$Signature = gsub("\\|", "", complete_results$Signature)
complete_results$Signature = gsub("\\.", "\\_", complete_results$Signature)
complete_results$Signature = gsub("ConsensusTME", "Leukocyte Subset ES", complete_results$Signature)
complete_results$Signature = gsub("Benci", "Expression Signature - Benci", complete_results$Signature)
complete_results$Signature = gsub("_", " ", complete_results$Signature)

df = complete_results

plot_df = ddply(df, .(df$Signature, df$Genus_full, df$`p value overall`, df$`rho overall`), nrow)
all_nodes = unique(c(plot_df$`df$Signature`,plot_df$`df$Genus_full`))
sankey_df = plot_df

colnames(sankey_df) = c("Source", "Target", "P-Val", "Rho", "Value")
sankey_df$Color = NA
names(all_nodes) = 0:(length(all_nodes)-1)

for (i in 1:nrow(sankey_df)){
  value = sankey_df$Source[i]
  new_value = colors[which(names(colors) == value)]
  sankey_df$Color[i] = new_value
}

if(colors_by_rho == "colors_by_rho"){
  sankey_df$Color[which(sankey_df$Rho < 0)] = "rgba(65, 131, 215, 0.3)"
  sankey_df$Color[which(sankey_df$Rho > 0)] = "rgba(217, 30, 24, 0.3)"
}

for (i in 1:nrow(sankey_df)){
  value = sankey_df$Source[i]
  new_value = names(all_nodes)[which(all_nodes == value)]
  sankey_df$Source[i] = new_value
}
for (i in 1:nrow(sankey_df)){
  value = sankey_df$Target[i]
  new_value = names(all_nodes)[which(all_nodes == value)]
  sankey_df$Target[i] = new_value
}

sankey_df$`P-Val` = paste0("p value = ", signif(sankey_df$`P-Val`, 3))

p <- plot_ly(
  type = "sankey",
  orientation = "h",
  arrangement = "snap",
  
  node = list(
    #label = all_nodes,
    label = c(""),
    #color = c(unique(sankey_df$Color), rep("white", length(unique(sankey_df$Target)))),
    color = "white",
    pad = 15,
    thickness = 6,
    line = list(
      color = "black",
      width = 1
    )
  ),
  
  link = list(
    source = sankey_df$Source,
    target = sankey_df$Target,
    value = sankey_df$Value,
    color = sankey_df$Color,
    label = sankey_df$`P-Val`
  )
) %>% 
  layout(
    title = "",
    font = list(
      size = 15,
      family = "arial",
      color = "black"
    )
  ) 

p


order_sign = c("Expression Signature - IL8 21978456",
               "Expression Signature - GRANS PCA 16704732",
               "Expression Signature - PDL1 data",
               "Expression Signature - TREM1 data",
               "Expression Signature - IFNG score 21050467",
               "Attractor Metagene - G LILRB4",
               "Attractor Metagene - G SIGLEC9",
               "Leukocyte Subset ES - Macrophages",
               "Leukocyte Subset ES - Macrophages M1",
               "Leukocyte Subset ES - Monocytes",
               "Leukocyte Subset ES - Dendritic cells",
               "Leukocyte Subset ES - Immune Score",
               "Leukocyte Subset ES - Mast cells",
               "Expression Signature - Chemokine12 score")

order_genus = c("Selenomonas (Firmicutes)", "Selenomonas 3 (Firmicutes)",
                "Blautia (Firmicutes)", "Eubacterium hallii group (Firmicutes)",
                "Ruminococcaceae UCG-013 (Firmicutes)", "Faecalibacterium (Firmicutes)",
                "Fusicatenibacter (Firmicutes)", "Subdoligranulum (Firmicutes)")


df$Signature = gsub("\\(ConsensusTME) ", "", df$Signature)
df$Signature = factor(df$Signature, levels = order_sign)
df$Genus_full = factor(df$Genus_full, levels = order_genus)

df = df[order(df$Signature),]
df = df[order(df$Genus_full),]

df$`p value overall` = -log10(df$`p value overall`)
df$`p value nonhypermutated` = -log10(df$`p value nonhypermutated`)
df$`p value hypermutated` = -log10(df$`p value hypermutated`)
df$`p value right sided` = -log10(df$`p value right sided`)
df$`p value left sided` = -log10(df$`p value left sided`)

write.csv(df, file = "./Analysis/Microbiome/019_Sankey_table_in_correct_order/April_2021_v5_019_Sankey_table_in_correct_order.csv",
          row.names = FALSE)
dir.create("./Analysis/Microbiome/019_Sankey_table_in_correct_order", showWarnings = FALSE)
save(df, file = "./Analysis/Microbiome/019_Sankey_table_in_correct_order/April_2021_v5_019_Sankey_table_in_correct_order.Rdata")

#########
# Not for main sankey plot, supplementary analysis
# Mirrored for normal tissue
masterfile = read.csv("./Analysis/Microbiome/017.2_Annotate_results/April_2021_0.1_v5_Annotated_Spearman_results.csv",
                      stringsAsFactors = FALSE)
colnames(masterfile)[2] = "Signature"
colnames(masterfile)[1] = "Genus_full"

masterfile$Phylum = gsub(".*\\D_1__", "", masterfile[, Rank])
masterfile$Phylum = gsub("\\ D_2_.*", "", masterfile$Phylum)
masterfile$Genus_full = gsub(".*\\D_5__", "", masterfile$Genus_full)
masterfile$Genus_full = gsub("\\[", "", masterfile$Genus_full)
masterfile$Genus_full = gsub("\\]", "", masterfile$Genus_full)
masterfile$Genus_full = paste0(masterfile$Genus_full, " (", masterfile$Phylum, ")")

masterfile$Signature = gsub("\\|", "", masterfile$Signature)
masterfile$Signature = gsub("\\.", "\\_", masterfile$Signature)
masterfile$Signature = gsub("ConsensusTME", "Leukocyte Subset ES", masterfile$Signature)
masterfile$Signature = gsub("Benci", "Expression Signature - Benci", masterfile$Signature)
masterfile$Signature = gsub("_", " ", masterfile$Signature)

df$paste = paste(df$Genus_full, df$Signature)

df2 = masterfile
df2$paste = paste(df2$Genus_full, df2$Signature)
df2 = df2[which(df2$paste %in% df$paste),]
  
#df2$Signature = gsub("\\(ConsensusTME) ", "", df2$Signature)
df2$Signature = factor(df2$Signature, levels = order_sign)
df2$Genus_full = factor(df2$Genus_full, levels = order_genus)

df2 = df2[order(df2$Signature),]
df2 = df2[order(df2$Genus_full),]

df2$p_val = -log10(df2$p_val)
df2$N_p_val = -log10(df2$N_p_val)
df2$delta_p_val = -log10(df2$delta_p_val)

write.csv(df2, file = "./Analysis/Microbiome/019_Sankey_table_in_correct_order/April_2021_v5_Masterfile_subset_in_correct_Sankey_order.csv",
          row.names = FALSE)

