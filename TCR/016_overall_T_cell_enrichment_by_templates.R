
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable = "T_cells"

# Load data
load("./Analysis/TCR/2_Corrections_for_DNA_input/TCR_Overview_DNA_input_corrected.Rdata")
load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/ConsensusTME_COAD_ES.Rdata"))

immune_sig_df = as.data.frame(t(ES))

i=1
for (i in 1:ncol(immune_sig_df)){
  col = colnames(immune_sig_df)[i]
  immune_sig_df[, col] = (immune_sig_df[, col] - min(immune_sig_df[, col]))/(max(immune_sig_df[,col])-min(immune_sig_df[,col]))
}

immune_sig_df$T_cells = NA
immune_sig_df$T_cells = (immune_sig_df$T_cells_CD4 + 
                           immune_sig_df$T_cells_CD8) / 2

TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
TCR_Overview$T_cells = immune_sig_df$T_cells[match(TCR_Overview$sample_name,
                                                   substring(rownames(immune_sig_df), 1, 4))]

plot2 = ggplot(TCR_Overview, aes(x = TCR_Overview[,variable], y = productive_templates)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  xlab(gsub("\\_", " ", variable)) +
  ylab("productive templates") +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson", size = 5) +
  geom_smooth(method="lm")


correlation = cor.test(x = TCR_Overview[,variable], y = TCR_Overview$productive_templates, method = "pearson")
n = 114
t = 0.6 / (sqrt((1-(0.6)^2)/(114-2)))
p = 2*pt(-abs(t),df=n-2)
-log10(p)

dir.create("./Figures/TCR/016_overall_T_cells_enrichment", showWarnings = FALSE)
png(paste0("./Figures/TCR/016_overall_T_cells_enrichment/016_Scatterplot_", variable, "_TCR_productive_templates.png"), res = 600,
    width = 3.2, height = 3, units = "in")
plot(plot2)
dev.off()
