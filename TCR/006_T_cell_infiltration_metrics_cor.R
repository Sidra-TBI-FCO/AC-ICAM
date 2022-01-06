
### T cell metrics

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters and load data
Source = "ConsensusTME_COAD"

load("./Analysis/TCR/2_Corrections_for_DNA_input/TCR_Overview_DNA_input_corrected.Rdata")
load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/", Source, "_ES.Rdata"))
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

ES = as.data.frame(ES)
ES["ICRscore",] = table_cluster_assignment$ICRscore[match(colnames(ES), rownames(table_cluster_assignment))]

ES = as.matrix(ES)
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]

variables = rownames(ES)
N.variables = length(variables)

results_ICR = data.frame(Variable = variables, cor_coef = NA)
results_productive_templates = data.frame(Variable = variables, cor_coef = NA)
results_productive_templates_cor = data.frame(Variable = variables, cor_coef = NA)
results_clonality = data.frame(Variable = variables, cor_coef = NA)

i = 20
for (i in 1:N.variables){
  variable = variables[i]
  # Append ICR Cluster data to the overview file
  TCR_Overview[,variable] = ES[variable, ][match(TCR_Overview$sample_name, substring(colnames(ES), 1, 4))]
  TCR_Overview[,variable] = as.numeric(TCR_Overview[,variable])
  
  plot = ggplot(TCR_Overview, aes(x = ICRscore, y = TCR_Overview[,variable])) +
    geom_point(aes(color = HLM_cluster), size = 0.8) +
    ylab(gsub("\\_", " ", variable)) +
    scale_color_manual(values = c("red", "green", "blue")) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15),
          legend.position = "none",
          aspect.ratio = 1/1) +
    labs(color = "ICR cluster") +
    stat_cor(method = "pearson") +
    geom_smooth(method="lm")
  
  correlation = cor(x = TCR_Overview$ICRscore, y = TCR_Overview[,variable], method = "pearson")
  results_ICR$cor_coef[which(results_ICR$Variable == variable)] = correlation
  
  dir.create("./Figures/TCR/", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_Scatterplot_ICR_", variable, ".png"), res = 600,
      width = 3.2, height = 3, units = "in")
  plot(plot)
  dev.off()
  
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
  
  correlation = cor(x = TCR_Overview[,variable], y = TCR_Overview$productive_templates, method = "pearson")
  results_productive_templates$cor_coef[which(results_ICR$Variable == variable)] = correlation
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_i_Scatterplot_", variable, "_TCR_productive_templates.png"), res = 600,
      width = 3.2, height = 3, units = "in")
  plot(plot2)
  dev.off()
 
  plot3 = ggplot(TCR_Overview, aes(x = TCR_Overview[,variable], y = productive_templates_cor)) +
    geom_point(aes(color = HLM_cluster), size = 0.8) +
    xlab(gsub("\\_", " ", variable)) +
    scale_color_manual(values = c("red", "green", "blue")) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15),
          legend.position = "none",
          aspect.ratio = 1/1) +
    labs(color = "ICR cluster") +
    stat_cor(method = "pearson") +
    geom_smooth(method="lm")
  
  correlation = cor(x = TCR_Overview[,variable], y = TCR_Overview$productive_templates_cor, method = "pearson")
  results_productive_templates_cor$cor_coef[which(results_ICR$Variable == variable)] = correlation
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_Scatterplot_", variable, "_productive_templates_cor.png"), res = 600,
      width = 3.2, height = 3, units = "in")
  plot(plot3)
  dev.off()
  
  plot4 = ggplot(TCR_Overview, aes(x = TCR_Overview[,variable], y = productive_clonality)) +
    geom_point(aes(color = HLM_cluster), size = 0.8) +
    xlab(gsub("\\_", " ", variable)) +
    ylab("productive clonality") +
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
    #xlab("") +
    #ylab("")
  
  correlation = cor(x = TCR_Overview[,variable], y = TCR_Overview$productive_clonality, method = "pearson")
  results_clonality$cor_coef[which(results_ICR$Variable == variable)] = correlation
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_i_new_Scatterplot_", variable, "_TCR_cell_clonality.png"), res = 600,
      width = 3.2, height = 3, units = "in")
  plot(plot4)
  dev.off()
}

dir.create("./Analysis/TCR/006_T_cell_infiltration_cor", showWarnings = FALSE)
write.csv(results_ICR, file = paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_ICR.csv"), row.names = FALSE)
write.csv(results_productive_templates, file = paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_productive_templates.csv"), row.names = FALSE)
write.csv(results_productive_templates_cor, file = paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_productive_templates_cor.csv"), row.names = FALSE)
write.csv(results_clonality, file = paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_clonality.csv"), row.names = FALSE)

ESTIMATE = data.frame(ESTIMATE)
rownames(ESTIMATE) = gsub("_P", "", rownames(ESTIMATE))
TCR_Overview$`ESTIMATE Immune` = ESTIMATE$ImmuneScore[match(TCR_Overview$sample_name, rownames(ESTIMATE))]
TCR_Overview$`ESTIMATE Score` = ESTIMATE$ESTIMATEScore[match(TCR_Overview$sample_name, rownames(ESTIMATE))]
TCR_Overview$`ESTIMATE Stromal` = ESTIMATE$StromalScore[match(TCR_Overview$sample_name, rownames(ESTIMATE))]

variable = "ESTIMATE Immune"

plot4 = ggplot(TCR_Overview, aes(x = get(variable), y = productive_templates)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  ylab("productive templates") +
  xlab(variable) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson") +
  geom_smooth(method="lm")


dir.create("./Figures", showWarnings = FALSE)
dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_Scatterplot_", variable, "_productive_templates.png"), res = 600,
    width = 3.2, height = 3, units = "in")
plot(plot4)
dev.off()

plot6 = ggplot(TCR_Overview, aes(x = TCR_Overview[,variable], y = productive_clonality)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  xlab(variable) +
  ylab("productive clonality") +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        aspect.ratio = 1/1,
        legend.position = "none") +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson", size = 6) +
  geom_smooth(method="lm") +
  xlab("") +
  ylab("")


dir.create("./Figures", showWarnings = FALSE)
png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_Scatterplot_", variable, "_Productive_clonality.png"), res = 600,
    width = 3.2, height = 3, units = "in")
plot(plot6)
dev.off()
