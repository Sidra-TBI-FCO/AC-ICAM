
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(required.packages = c("ggplot2", "dplyr"))

# Set parameters
exclude_medium = "exclude_medium" # "include_medium" or "exclude_medium" or "include_medium_cat"
stage_filter = ""

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/frequency_df.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

Merged_dataset$ICRscore = table_cluster_assignment$ICRscore[match(Merged_dataset$Patient_ID,
                                                                  substring(rownames(table_cluster_assignment), 1, 3))]
Merged_dataset$ICR_cluster = table_cluster_assignment$ICR_HML[match(Merged_dataset$Patient_ID,
                                                                    substring(rownames(table_cluster_assignment), 1, 3))]

Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

dim(RNASeq.QN.LOG2)
table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(RNASeq.QN.LOG2)),]
median_ICRscore = median(table_cluster_assignment$ICRscore)

if(exclude_medium == "exclude_medium"){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$ICR_cluster == "ICR Medium"),]
  color_vals = c("ICR High" = "red", "ICR Low" = "blue")
}

if(exclude_medium == "include_medium"){
  Merged_dataset$ICR_by_median = NA
  Merged_dataset$ICR_by_median[which(Merged_dataset$ICRscore < median_ICRscore)] = "ICR Low"
  Merged_dataset$ICR_by_median[which(Merged_dataset$ICRscore >= median_ICRscore)] = "ICR High"
  table(Merged_dataset$ICR_by_median, Merged_dataset$Immunoedited)
  Merged_dataset$ICR_cluster = Merged_dataset$ICR_by_median
  color_vals = c("ICR High" = "red", "ICR Low" = "blue")
}

if(exclude_medium == "include_medium_cat"){
  # No need for changing ICR subgroups
  color_vals = c("ICR High" = "red", "ICR Medium" = "lightgrey", "ICR Low" = "blue")
}

table(Merged_dataset$ajcc_pathologic_tumor_stage, exclude = NULL)
if(stage_filter == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage == stage_filter),]
}

TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
Merged_dataset$TCR_productive_clonality = TCR_Overview$productive_clonality[match(Merged_dataset$Patient_ID,
                                                                                  TCR_Overview$Patient_ID)]

  
# Scatterplot ICR
scatterplot = ggplot(data = Merged_dataset, aes(x = Immunoediting_score, y = ICRscore)) +
  geom_point(aes(color = ICR_cluster)) +
  scale_color_manual(values = color_vals) +
  #stat_cor(method = "pearson", size = 6) +
  #scale_y_log10() +
  theme_bw() +
  ylab("ICR score") +
  xlab("Immunoediting score") +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) #+
  #geom_smooth(method="lm")


png(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/022c_Scatterplot_ICRscore_by_immunoediting_status_",
           exclude_medium, stage_filter, ".png"),
    res = 600, height = 5, width = 5, units = "in")
plot(scatterplot)
dev.off()

Merged_dataset$Immunoedited_ICR = paste(Merged_dataset$Immunoedited, Merged_dataset$ICR_cluster)
# Scatterplot alternative
scatterplot = ggplot(data = Merged_dataset, aes(x = Immunoediting_score, y = ICRscore)) +
  geom_point(aes(color = Immunoedited_ICR)) +
  scale_color_manual(values = c("immunoedited ICR High" = "#FF3806",
                                "less immunoedited ICR High" = "#FF9F00",
                                "immunoedited ICR Low" = "#009CFC",
                                "less immunoedited ICR Low" = "#5233FC")) +
  #stat_cor(method = "pearson", size = 6) +
  #scale_y_log10() +
  theme_bw() +
  ylab("ICR score") +
  xlab("Immunoediting score") +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        aspect.ratio = 1/1) #+

png(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/022c_Alternative_Scatterplot_ICRscore_by_immunoediting_status_",
           exclude_medium, stage_filter, ".png"),
    res = 600, height = 5, width = 5, units = "in")
plot(scatterplot)
dev.off()


# Boxplot ICR
plot = ggplot(Merged_dataset, aes(x = Immunoedited, y = ICRscore)) +
  geom_boxplot(outlier.shape = NA, aes(color = Immunoedited)) +
  scale_color_manual(values = c("immunoedited" = "#FF5800", "less immunoedited" = "#8085E9")) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  xlab("") +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(c("immunoedited", "less immunoedited")))
  
png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/022c_boxplot_ICRscore_by_immunoediting_status.png",
    res = 600, width = 3, height = 4.3, units = "in")
plot(plot)  
dev.off()

TCR_subset = Merged_dataset[-which(is.na(Merged_dataset$TCR_productive_clonality)),]
# O/E Immunoediting score
scatterplot = ggplot(data = TCR_subset, aes(x = Immunoediting_score, y = TCR_productive_clonality)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  #scale_y_log10() +
  theme_bw() +
  ylab("TCR clonality") +
  xlab("Immunoediting score") +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")
  

png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/022c_Scatterplot_TCR_productive_clonality_by_immunoediting_status.png",
    res = 600, height = 5, width = 5, units = "in")
plot(scatterplot)
dev.off()

boxplot = ggplot(TCR_subset, aes(x = Immunoedited, y = TCR_productive_clonality)) +
  geom_boxplot(outlier.shape = NA, aes(color = Immunoedited)) +
  scale_color_manual(values = c("immunoedited" = "#FF5800", "less immunoedited" = "#8085E9")) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  xlab("") +
  ylab("TCR productive clonality") +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(c("immunoedited", "less immunoedited")))

png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/022c_Boxplot_TCR_productive_clonality_by_immunoediting_status.png",
    res = 600, height = 4.3, width = 3, units = "in")
plot(boxplot)
dev.off()

### Stacked barchart CMS
df_plot = Merged_dataset[-which(is.na(Merged_dataset$CMS)),]
#df_plot = df_plot[which(df_plot$ajcc_pathologic_tumor_stage %in% c("3", "4")),]
DF = df_plot %>%
  group_by(Immunoedited, CMS) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

colors = c("CMS1" = "#FF9F21",  "CMS2" = "#0074AF", 
           "CMS3" = "#E97AA8", "CMS4" = "#009E74")

stacked_barplot = ggplot(DF, aes(x = Immunoedited, y =perc, fill = CMS)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "CMS", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= colors)

png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/022c_Stacked_barplot_immunoediting_by_CMS_stage.png",
    res = 600, height = 4.3, width = 3, units = "in")
plot(stacked_barplot)
dev.off()

save(Merged_dataset, file = paste0("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset_", exclude_medium, ".Rdata"))

