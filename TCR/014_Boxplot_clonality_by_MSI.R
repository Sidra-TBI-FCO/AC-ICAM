
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable = "productive_clonality"

# Load data
load("./Analysis/TCR/2_Corrections_for_DNA_input/TCR_Overview_DNA_input_corrected.Rdata")
load("./Analysis/Trimmed_p/Deconvolution_and_GSEA/Bindea_ORIG_ES.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

TCR_Overview$HLM_cluster = factor(TCR_Overview$HLM_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
levels(TCR_Overview$HLM_cluster) = c("Low", "Medium", "High")
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
TCR_Overview$CMS = Rfcms$RF.predictedCMS[match(TCR_Overview$Patient_ID, substring(rownames(Rfcms), 1, 3))]
TCR_Overview$CMS = factor(TCR_Overview$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
TCR_Overview$MSI = MANTIS$MSI[match(TCR_Overview$Patient_ID, MANTIS$Patient_ID)]
TCR_Overview$MSI = factor(TCR_Overview$MSI, levels = c("MSS", "MSI-H"))
TCR_Overview = TCR_Overview[-which(is.na(TCR_Overview$MSI)),]

# Plotting
df_plot = TCR_Overview[, c("ICRscore", "MSI", variable)]

plot = ggplot(df_plot, aes(x = MSI, y = productive_clonality)) +
  #geom_boxplot(outlier.shape = NA, aes(fill = df_plot$MSI)) +
  scale_fill_manual(values = c("MSI-H" = "purple",  "MSS" = "white")) +
  #geom_jitter(width = 0.2, size = 1.1) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.position = "none") +
  stat_compare_means(method = "t.test", comparisons = list(c("MSS", "MSI-H")),
  label = "p.signif") +
  ylab(gsub("_", " ", variable)) +
  xlab("") +
  ylim(0, 0.4)

violin_plot = plot + geom_violin(aes(fill = df_plot$MSI)) +
  geom_boxplot(width=.1, outlier.shape = NA, aes(fill = df_plot$MSI)) +
  ylab("") +
  theme(axis.text.x = element_blank())

dir.create("./Figures/TCR/014_Boxplot_clonality_MSI", showWarnings = FALSE)
png(filename = "./Figures/TCR/014_Boxplot_clonality_MSI/Boxplot_clonality_MSI_status.png",
    width = 3, height = 3.8, units = "in", res = 600)
plot(plot)
dev.off()

png(filename = "./Figures/TCR/014_Boxplot_clonality_MSI/Violin_plot_clonality_MSI_status_no_stats.png",
    width = 3, height = 3.8, units = "in", res = 600)
plot(violin_plot)
dev.off()
