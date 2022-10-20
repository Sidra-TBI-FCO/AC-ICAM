

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# BiocManager::install("survcomp")
ipak(c("ggplot2", "survival", "survivalAnalysis", "data.table", 
       "corrplot", "creditmodel", "dplyr", "glmnet", "doParallel", "extrafont",
       "extrafontdb", "stringr", "stringi", "factoextra", "ggfortify", "survminer",
       "gbm", "survcomp"))

loadfonts()
registerDoParallel(cores=10)
source("./R code/Microbiome/Risk_model_glmnet/get_functions.R")

# Load data
load("./Analysis/Microbiome/glmnet_model_microbiome/001_Prepared_data_AC_ICAM_246_for_glmnet.Rdata")

#Perform CV
# cvfit_list <- perform_CV(rev_survival_df = rev_survival_df, train_data = train_data, colnames_of_interest = colnames_of_interest, y_cols = y_cols)

#Get the best CV CI #Here CV stands for cross-validation
load("./Processed_Data/Microbiome/All_Input_Data_for_Risk_Model/CV_Performance.Rdata")
cv_CI <- NULL
for (i in 1:10)
{
  best_CV_CI <- max(cvfit_list[[i]]$cvm)
  cv_CI <- c(cv_CI,best_CV_CI)
}

#Make the boxplot for nested CV 
boxplot(cv_CI,names=c("CV CI"))

#Get the best cv model closest to median
best_cv_index <- which.min(abs(cv_CI-as.numeric(summary(cv_CI)[5])))
best_cvfit <- cvfit_list[[best_cv_index]]

#Build the survival model data frame
train_out <- get_predictions(rev_survival_df, colnames_of_interest, y_cols, best_cvfit)
y_pred_test <- as.numeric(train_out[[1]])

#Build the model for risk score
median_y_test <- mean(y_pred_test)
output_df <- build_risk_score_df(rev_survival_df, y_pred_test, median_y_test)

dir.create("./Figures/Microbiome/glmnet_model_microbiome", showWarnings = FALSE)
dir.create("./Figures/Microbiome/glmnet_model_microbiome/002_Predictions", showWarnings = FALSE)

save(best_cvfit, median_y_test, file = "./Analysis/Microbiome/glmnet_model_microbiome/best_cv_fit.Rdata")

#Make the plot of risk score
g_output <- ggplot(data=output_df,aes(x=samples, y=prediction, color=risk_cat)) + 
  scale_color_manual(labels = c("Low Risk","High Risk"),values=c("blue","red"))+
  geom_point(size=2) + theme_classic() + 
  theme(axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 11),
        legend.position = "bottom", legend.text = element_text(size=11))
ggsave(filename = "./Figures/Microbiome/glmnet_model_microbiome/002_Predictions/Training_Based_Risk_Categories.pdf", 
       plot=g_output, device = pdf(), height = 4, width = 5, units = "in", dpi=300)
dev.off()

#PCA 
res.pca <- prcomp(rev_survival_df[,c(4:ncol(rev_survival_df))], scale = F, center = F)
groups <- as.factor(as.numeric(as.vector(y_pred_test<median_y_test)))
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("blue",  "red"),
             addEllipses = F, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Risk Category",
             repel = TRUE)

#Make the KM plot
g <- make_km_plot(output_df, years=10)
ggsave(filename="./Figures/Microbiome/glmnet_model_microbiome/002_Predictions/KM_Plot_Training_Full.pdf",plot=g$plot, device = pdf(), height = 6, width = 8, units="in", dpi=300)
dev.off()

dir.create("./Analysis/Microbiome/glmnet_model_microbiome/002_Predictions", showWarnings = FALSE)
save(output_df, file = "./Analysis/Microbiome/glmnet_model_microbiome/002_Predictions/AC_ICAM_246_Predictions.Rdata")

