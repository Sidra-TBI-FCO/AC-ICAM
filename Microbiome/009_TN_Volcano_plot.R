

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggrepel"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)

# Load data
load(paste0("./Analysis/Microbiome/007_Paired_wilcoxon_test_TN/",  Type, "/", Rank,
            "/Wilcoxon_test_", Type, "_", Rank, ".Rdata"))

# Make Volcano plot
results$Name = as.character(results$Name)
results$p_value_log10 = -log10(results$p_val)
results$color = "non-significant"
results$color[results$FDR < 0.1 & results$Direction == "Enriched in tumor"] = "FDR < 0.1 enriched in Tumor"
results$color[results$FDR < 0.1 & results$Direction == "Enriched in normal"] = "FDR < 0.1 enriched in Normal"

#results = results[order(results$ratio, decreasing = TRUE),]
results$label = gsub(".*D_5__", "", results$Name)
results$Phylum = gsub(".*D_1__", "", results$Name)
results$Phylum = gsub("\\ D_2__.*", "", results$Phylum)
results$label = paste(results$label, " (", results$Phylum, ")", sep = "")
results$Plot_label = results$label
results$Plot_label[-which(results$FDR < 0.0001)] = NA
results$Plot_label[which(results$log_ratio_paired < 0.5)] = NA
#results$Plot_label[which(results$p_value_log10 < 10)] = NA
#results$Plot_label[which(results$log_ratio > 0.01)] = results$label[which(results$log_ratio > 0.01)]


results = results[-which(results$Name == "Unknown"),]
#results = results[-which(results$mean_N < 0.000001 & results$mean_T < 0.00001),] # delete low abundance rank level microbiota


plotA = ggplot(results, aes(x= mean_ratio_paired, y= p_value_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "FDR < 0.1 enriched in Tumor" = "#1E90FF", "FDR < 0.1 enriched in Normal" = "#FF00FF")) +
  theme_bw() +
  #xlim(-2, 2) +
  ylab("-log10 pvalue") +
  xlab("Ratio (log10)") +
  geom_line(y = 1.3, linetype = "dashed", color = "darkgrey") +
  #geom_line(y = 2.1, linetype = "dashed", color = "darkgrey") + # change this for final version (FDR = 0.05)
  geom_text_repel(aes(x = mean_ratio_paired, p_value_log10, 
                      label = Plot_label), size = 2.5,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

dir.create("./Figures/Microbiome/009_TN_Volcano", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/009_TN_Volcano/009_TN_Volcano_mean_ratio_paired_", Type, "_abundance_", Rank, ".png"), width = 5, height = 4, units = "in", res = 600)
plot(plotA)
dev.off()

results$Ratio_mean_log10 = log(results$Ratio_mean,10)
results$Plot_label[which(results$Ratio_mean_log10 < 0.5 & results$p_value_log10 < 10)] = NA

plotB = ggplot(results, aes(x= Ratio_mean_log10, y= p_value_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "FDR < 0.1 enriched in Tumor" = "#1E90FF", "FDR < 0.1 enriched in Normal" = "#FF00FF")) +
  theme_bw() +
  xlim(-2, 2) +
  ylab("-log10 pvalue") +
  xlab("Log10 ratio \n (mean tumor / mean normal)") +
  geom_line(y = 1.3, linetype = "dashed", color = "darkgrey") +
  geom_line(y = 2.1, linetype = "dashed", color = "darkgrey") + # change this for final version (FDR = 0.05)
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "darkgrey") +
  #geom_line(x = 0.5, linetype = "dashed", color = "darkgrey") +
  #geom_line(x = -0.5, linetype = "dashed", color = "darkgrey") +
  geom_text_repel(aes(x = Ratio_mean_log10, p_value_log10, 
                      label = Plot_label), size = 0,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

dir.create("./Figures/Microbiome/009_TN_Volcano", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/009_TN_Volcano/009_TN_Volcano_Ratio_means_log10_", Type, "_abundance_", Rank, ".png"), width = 5, height = 4, units = "in", res = 600)
plot(plotB)
dev.off()


plot2 = ggplot(results, aes(x= mean_delta_paired, y= p_value_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "FDR < 0.1 enriched in Tumor" = "#1E90FF", "FDR < 0.1 enriched in Normal" = "#FF00FF")) +
  theme_bw() +
  xlim(-0.045, 0.045) +
  ylab("-log10 pvalue") +
  xlab("Difference (delta) between proportions") +
  geom_text_repel(aes(x = mean_delta_paired, p_value_log10, 
                      label = Plot_label), size = 2.5,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

dir.create("./Figures/Microbiome/009_TN_Volcano", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/009_TN_Volcano/009_TN_Volcano_mean_delta_paired_", Type, "_abundance_", Rank, ".png"), width = 5, height = 4, units = "in", res = 600)
plot(plot2)
dev.off()

results$Proportion_difference_log = log(results$Proportion_difference, 10)
# Negative proportion difference cannot be log transformed!

plot3 = ggplot(results, aes(x= Proportion_difference_log, y= p_value_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "FDR < 0.1 enriched in Tumor" = "#1E90FF", "FDR < 0.1 enriched in Normal" = "#FF00FF")) +
  theme_bw() +
  xlim(-2, 2) +
  ylab("-log10 pvalue") +
  xlab("Log difference (delta) between proportions") +
  geom_text_repel(aes(x = Proportion_difference_log, p_value_log10, 
                      label = Plot_label), size = 2.5,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

dir.create("./Figures/Microbiome/009_TN_Volcano", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/009_TN_Volcano/009_TN_Volcano_Delta_log_", Type, "_abundance_", Rank, ".png"), width = 5, height = 4, units = "in", res = 600)
plot(plot3)
dev.off()

