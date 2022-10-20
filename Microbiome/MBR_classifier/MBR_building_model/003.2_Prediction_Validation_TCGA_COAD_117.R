
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "survival", "survivalAnalysis", "data.table", 
       "corrplot", "creditmodel", "dplyr", "glmnet", "doParallel", "extrafont",
       "extrafontdb", "stringr", "stringi", "factoextra", "ggfortify", "survminer",
       "gbm", "survcomp"))

source("./R code/Microbiome/Risk_model_glmnet/get_functions.R")

# Load data
load("./Analysis/Microbiome/glmnet_model_microbiome/best_cv_fit.Rdata")
load("./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_TCGA_COAD_117_for_glmnet.Rdata")

#Build the survival model data frame
tcga_out <- get_predictions(rev_survival_df = tcga_survival_df, colnames_of_interest = colnames_of_interest, y_cols = y_cols, best_cvfit = best_cvfit)
tcga_y_pred <- as.numeric(tcga_out[[1]])

tcga_output_df <- build_risk_score_df(rev_survival_df = tcga_survival_df, tcga_y_pred, median_y_test)

#Build model for 5 years
g_tcga <- make_km_plot(tcga_output_df, years = 5)
print(g_tcga)

ggsave(filename = "./Figures/Microbiome/glmnet_model_microbiome/002_Predictions/KM_Plot_TCGA_Full_5years.pdf", plot = g_tcga$plot, device = pdf(), units="in", height=6, width=8, dpi=300)
dev.off()

save(tcga_output_df, file = "./Analysis/Microbiome/glmnet_model_microbiome/002_Predictions/002_Predictions_TCGA_COAD_117_Predictions.Rdata")
