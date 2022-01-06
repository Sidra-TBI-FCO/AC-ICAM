
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

# Load data
load("./Analysis/TCR/7_Overlapping_sequences_analysis/7_Overlapping_sequences_TN_pairs.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

# Set parameter
variable = "Fraction_Tumor"

# Analysis
results$Total_sequences = results$Number_Tumor_Restricted + results$Number_Normal_Restricted + results$Number_Overlap
results$Fraction_Overlap = results$Number_Overlap / results$Total_sequences
results$Fraction_Tumor = results$Number_Tumor_Restricted / results$Total_sequences

TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]

# Plot
plot_df = data.frame(Patient = results$Patient, Variable = results[, variable], ICRscore = NA, ICR_cluster = NA)
plot_df$Patient = as.character(plot_df$Patient)
plot_df$ICRscore = table_cluster_assignment$ICRscore[match(plot_df$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
plot_df$ICR_cluster = table_cluster_assignment$ICR_HML[match(plot_df$Patient, substring(rownames(table_cluster_assignment), 1, 3))]

plot = ggplot(plot_df, aes(x = plot_df$ICRscore, y = plot_df$Variable, label = plot_df$Patient)) +
  geom_point(aes(color = plot_df$ICR_cluster)) +
  scale_color_manual(values = c("red", "green", "blue")) +
  geom_abline() +
  ylab(gsub("_", " ", variable)) +
  xlab("ICR score") +
  theme_bw() +
  stat_cor(method = "pearson") +
  geom_smooth(method="lm") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black")) +
  geom_text(hjust = 0, nudge_x = 0.05) +
  xlim(3, 8.5)

dev.new()
plot(plot)
