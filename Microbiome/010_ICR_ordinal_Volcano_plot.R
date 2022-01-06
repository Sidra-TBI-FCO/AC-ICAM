

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggrepel"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Species_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)

# Load data
load(paste0("./Analysis/Microbiome/008_ICR_cluster_Spearman/",  Type, "/", Rank,
            "/Spearman_by_ICR_cluster_", Type, "_", Rank, ".Rdata"))

# Make Volcano plot
results$Name = as.character(results$Name)
results$p_value_log10 = -log10(results$p_val)
results$color = "non-significant"
results$color[results$p_val < 0.05 & results$rho > 0] = "p value < 0.05 correlate with ICR High"
results$color[results$p_val < 0.05 & results$rho < 0] = "p value < 0.05 correlate with ICR Low"

#results = results[order(results$ratio, decreasing = TRUE),]
results$label = gsub(".*D_6__", "", results$Name)
results$Phylum = gsub(".*D_1__", "", results$Name)
results$Phylum = gsub("\\ D_2__.*", "", results$Phylum)
results$label = paste(results$label, " (", results$Phylum, ")", sep = "")
results$Plot_label = results$label
results$Plot_label[-which(results$p_val < 0.01)] = NA
results$Plot_label[which(results$p_val < 0.05 & results$rho > 0)] = results$label[which(results$p_val < 0.05 & results$rho > 0)]
  
results = results[-which(results$Name == "Unknown"),]

plot = ggplot(results, aes(x= rho, y= p_value_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "p value < 0.05 correlate with ICR High" = "red", "p value < 0.05 correlate with ICR Low" = "blue")) +
  theme_bw() +
  xlim(-0.3, 0.3) +
  ylab("-log10 pvalue") +
  xlab("Spearman correlation coefficient") +
  #scale_x_log10() +
  geom_text_repel(aes(x = rho, p_value_log10, 
                      label = Plot_label), size = 1.5,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

dir.create("./Figures/Microbiome/010_ICR_Spearman_ordinal_Volcano", showWarnings = FALSE)
png("./Figures/Microbiome/010_ICR_Spearman_ordinal_Volcano/010_ICR_Spearman_ordinal_Volcano.png", width = 5, height = 4, units = "in", res = 600)
plot(plot)
dev.off()
