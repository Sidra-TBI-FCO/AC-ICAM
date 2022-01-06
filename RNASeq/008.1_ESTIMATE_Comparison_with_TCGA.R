
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable_of_interest = "StromalScore" #"ImmuneScore" #"StromalScore" #"ESTIMATEScore"
download.method = "Biolinks"
intersect = ""

# Create directories
dir.create("./Figures/Trimmed_p/ESTIMATE", showWarnings = FALSE)

# Load data
# Load data
if(download.method == "Biolinks"){
  load("../NGS_Data_TCGA_COAD_Jessica/Analysis/008_ESTIMATE/Biolinks_TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")
}else{
  load(paste0("./Analysis/Trimmed_p/ESTIMATE/", intersect, "TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata"))
}
ESTIMATE_tcga = ESTIMATE

load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")

ESTIMATE = as.data.frame(ESTIMATE)
ESTIMATE_tcga = as.data.frame(ESTIMATE_tcga)

ESTIMATE$cohort = "JSREP"
ESTIMATE_tcga$cohort = "TCGA"
df_plot = rbind(ESTIMATE, ESTIMATE_tcga)

my_comparisons = list(c("JSREP", "TCGA"))
plot = ggboxplot(df_plot, x = "cohort", y = variable_of_interest,
                 color = "cohort", palette = c("JSREP" = "#FFC000", "TCGA" = "#006DFF"),
                 add = "jitter", outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

png(paste0("./Figures/Trimmed_p/ESTIMATE/", download.method, "_", variable_of_interest, "_JSREP_vs_TCGA.png"), res = 600, height = 6, width = 6, units = "in")
plot(plot)
dev.off()
