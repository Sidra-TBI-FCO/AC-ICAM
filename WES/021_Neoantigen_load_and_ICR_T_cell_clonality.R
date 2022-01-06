
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Load data
load("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")

# plot
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
plot_df = frequency_df

plot_df$ICRscore = table_cluster_assignment$ICRscore[match(plot_df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
plot_df$ICR_cluster = table_cluster_assignment$ICR_HML[match(plot_df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
plot_df$ICR_cluster = factor(plot_df$ICR_cluster, levels = c("ICR Low","ICR Medium",  "ICR High"))
plot_df$TCR_clonality = TCR_Overview$productive_clonality[match(plot_df$Patient_ID, TCR_Overview$Patient_ID)]

dir.create("./Figures/WES/021_Neoantigen_load_and_immune_signature_correlation", showWarnings = FALSE)

plot = ggplot(data = plot_df, aes(x = ICRscore, y = Neoantigen_count)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  scale_y_log10() +
  theme_bw() +
  xlab("ICR score") +
  ylab("Neoantigen count \n Median WT score > 500 nM") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

png("./Figures/WES/021_Neoantigen_load_and_immune_signature_correlation/Neoantigen_load_ICRscore.png",
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()

plot_df = plot_df[-which(plot_df$Neoantigen_count == 0),] # remove the 0 sample to show statistics

plot = ggplot(plot_df, aes(x = ICR_cluster, y = Neoantigen_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.4, aes(color = ICR_cluster)) +
  scale_color_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
  #scale_color_manual(values = MANTIS_colors) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  stat_compare_means(method = "t.test", comparisons = list(c("ICR High", "ICR Medium"),
                                                           c("ICR High", "ICR Low")),
                     label = "p.signif") +
  ylab("Neoantigen count") +
  xlab("")

png("./Figures/WES/021_Neoantigen_load_and_immune_signature_correlation/Neoantigen_load_ICR_cluster.png",
    res = 600, width = 4, height = 5, units = "in")
plot(plot)
dev.off()

plot_df = plot_df[-which(is.na(plot_df$TCR_clonality)),]

plot = ggplot(data = plot_df, aes(x = TCR_clonality, y = Neoantigen_count)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  scale_y_log10() +
  theme_bw() +
  xlab("TCR productive clonality") +
  ylab("Neoantigen count \n Median WT score > 500 nM") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

png("./Figures/WES/021_Neoantigen_load_and_immune_signature_correlation/TCR_clonality_Neoantigen_count.png",
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()
