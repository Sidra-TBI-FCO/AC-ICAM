
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable_of_interest = "ImmuneScore" #"StromalScore" #"ImmuneScore" #"ESTIMATEScore"
intersect = ""
download.method = "Biolinks"

# Create directories
dir.create("./Figures/Trimmed_p/ESTIMATE", showWarnings = FALSE)

# Load data
if(download.method == "Biolinks"){
  load("../NGS_Data_TCGA_COAD_Jessica/Analysis/008_ESTIMATE/Biolinks_TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")
}else{
  load(paste0("./Analysis/Trimmed_p/ESTIMATE/", intersect, "TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata"))
}
ESTIMATE_tcga = ESTIMATE

load(paste0("./Analysis/Trimmed_p/ESTIMATE/", intersect, "JSREP_clean_ESTIMATE_scores.Rdata"))

ESTIMATE = as.data.frame(ESTIMATE)
ESTIMATE_tcga = as.data.frame(ESTIMATE_tcga)

ESTIMATE$cohort = "Sidra-LUMC cohort"
ESTIMATE_tcga$cohort = "TCGA"
df_plot = rbind(ESTIMATE, ESTIMATE_tcga)

plot = ggplot(df_plot, aes(x = df_plot[,variable_of_interest], color = cohort)) +
  geom_density(aes(fill = cohort), alpha = .4) +
  #xlab(paste0(variable_of_interest)) +
  scale_color_manual(values = c("#ff8c00", "#006DFF")) +
  scale_fill_manual(values = c("#ff8c00", "#006DFF")) +
  theme_bw() +
  xlab(variable_of_interest) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15))

png(paste0("./Figures/Trimmed_p/ESTIMATE/", download.method, "_", intersect, "Density_plot_", variable_of_interest, "_JSREP_vs_TCGA.png"), res = 600, height = 4, width = 5, units = "in")
plot(plot)
dev.off()
