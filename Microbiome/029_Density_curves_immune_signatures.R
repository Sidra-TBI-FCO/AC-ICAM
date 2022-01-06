
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Set parameters

# Load data
load(paste0("./Analysis/Microbiome/017_Spearman_continous/Redundant_immune_sig_df_used_as_input.Rdata"))

signatures = colnames(immune_sig_df)

dir.create("./Figures/Microbiome/029_Density_curves_immune_signatures", showWarnings = FALSE)
i=1
for (i in 1:length(signatures)){
  signature = signatures[i]
  df = immune_sig_df
  df[,signature] = as.numeric(df[,signature])
  plot = ggplot(df, aes(x = get(signature))) +
    geom_density() +
    xlab(signature) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.x = element_text(color = "black", size = 14),
          axis.title.y = element_text(color = "black", size = 14))
  
  png(paste0("./Figures/Microbiome/029_Density_curves_immune_signatures/", signature, "_density.png"),
      res=600, units = "in", width = 4, height = 4)
  plot(plot)
  dev.off()
}