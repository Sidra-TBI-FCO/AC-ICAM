
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
load("./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_AC_ICAM_42_for_glmnet.Rdata")

#Build the survival model data frame
subset_42_out <- get_predictions(rev_survival_df = rev_subset_42_survival_df, colnames_of_interest = colnames_of_interest, 
                                 y_cols = y_cols, best_cvfit = best_cvfit)
subset_42_y_pred <- as.numeric(subset_42_out[[1]])

#Builds the risk score dataframe
subset_42_output_df <- build_risk_score_df(rev_survival_df = rev_subset_42_survival_df, y_pred_test = subset_42_y_pred, median_y_test)

#See the model fit
subset_42_output_df$status <- as.numeric(as.numeric(subset_42_output_df$time<5) & subset_42_output_df$status)
g_subset_42 <- make_km_plot(subset_42_output_df, years=5)
ggsave(filename = "./Figures/Microbiome/glmnet_model_microbiome/002_Predictions/KM_Plot_Validation_42_Full_5years.pdf", plot = g_subset_42$plot, device = pdf(), units="in", height=6, width=8, dpi=300)
dev.off()

save(subset_42_output_df, file = "./Analysis/Microbiome/glmnet_model_microbiome/002_Predictions/AC_ICAM_42_Predictions.Rdata")
