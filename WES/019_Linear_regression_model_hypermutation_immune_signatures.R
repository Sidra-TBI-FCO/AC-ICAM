

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c()
ipak(required.packages)

# Set parameters
Gene.set = "ConsensusTME_COAD_WJedit_source"
co_variate = "ESTIMATE_Immune" # "ESTIMATE_Immune" or "ICRscore"

# Load data
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/", Gene.set, "_ES.Rdata"))

# Data preparation
ESTIMATE = data.frame(ESTIMATE)
frequency_df$Patient_ID = as.character(frequency_df$Patient_ID)

plot_df = frequency_df
plot_df$ESTIMATE_Immune = ESTIMATE$ImmuneScore[match(plot_df$Patient_ID, 
                                                     substring(rownames(ESTIMATE), 1, 3))]

plot_df$ICRscore = table_cluster_assignment$ICRscore[match(plot_df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]

plot_df$Mutation_cat = NA
plot_df$Mutation_cat[which(plot_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
plot_df$Mutation_cat[which(plot_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"
plot_df$Mutation_cat = factor(plot_df$Mutation_cat, levels = c("nonhypermutated", "hypermutated"))

# Analysis
variables = rownames(ES)

# Logistic Regression
# where F is a binary factor and
# x1-x3 are continuous predictors
# fit <- glm(F~x1+x2+x3,data=mydata,family=binomial())

results = data.frame(Signature = variables, p_value = NA, estimate = NA)

i=22
for (i in 1:length(variables)){
  var = variables[i]
  plot_df$ES = ES[var,][match(plot_df$Patient_ID, substring(colnames(ES), 1, 3))]
  fit = glm(Mutation_cat ~ ES + get(co_variate), data = plot_df, family = "binomial")
  table = summary(fit)
  test_table = data.frame(table$coefficients)
  results$p_value[which(results$Signature == var)] = test_table["ES","Pr...z.."]
  results$estimate[which(results$Signature == var)] = test_table["ES","Estimate"]
}

dir.create("./Analysis/WES/019_Linear_regression_model_hypermutation_status", showWarnings = FALSE)
write.csv(results, file = "./Analysis/WES/019_Linear_regression_model_hypermutation_status/results.csv",
          row.names = FALSE)
