
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak("ggpubr")

# Load data
load("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")

# Set parameters
subset = ""  #"nonhypermutated" "hypermutated"

# Hypermutation status
frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

# Ratio calculation
frequency_df$Ratio = frequency_df$Neoantigen_count / frequency_df$Non_silent_Mutation_frequency

frequency_df = frequency_df
frequency_df$ICRscore = table_cluster_assignment$ICRscore[match(frequency_df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
frequency_df$ICR_cluster = table_cluster_assignment$ICR_HML[match(frequency_df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
frequency_df$ICR_cluster = factor(frequency_df$ICR_cluster, levels = c("ICR Low","ICR Medium",  "ICR High"))
frequency_df$Pathologic_stage = clinical_data$ajcc_pathologic_tumor_stage[match(frequency_df$Patient_ID, clinical_data$Patient_ID)]

if(subset == ""){}else{
  frequency_df = frequency_df[which(frequency_df$Mutation_cat == subset),]
}

dir.create("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration", showWarnings = FALSE)

plot = ggplot(data = frequency_df, aes(x = ICRscore, y = Ratio)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  #scale_y_log10() +
  theme_bw() +
  xlab("ICR score") +
  ylab("Ratio \n Neoantigen count / Mutational load") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm") +
  ylim(0, 0.34)

png(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/", subset, "_Ratio_Neoantigens_Mutational_load_by_ICRscore.png"),
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()

if(subset == "hypermutated"){}else{
  frequency_df = frequency_df[-which(frequency_df$Neoantigen_count == 0),] # remove the 0 sample to show statistics
}

plot = ggplot(data = frequency_df, aes(x = Pathologic_stage, y = Ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.4, aes(color = ICR_cluster)) +
  scale_color_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
  #scale_color_manual(values = MANTIS_colors) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  ylab("Ratio \n Neoantigen count / Mutational load") +
  xlab("Pathologic stage") 

png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/Ratio_Neoantigens_Mutational_load_by_Stage.png",
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()

plot = ggplot(data = frequency_df, aes(x = Pathologic_stage, y = Neoantigen_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.4, aes(color = ICR_cluster)) +
  scale_color_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
  #scale_color_manual(values = MANTIS_colors) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  ylab("Neoantigen count") +
  xlab("Pathologic stage")

png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/Neoantigens_by_Stage.png",
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()

plot = ggplot(data = frequency_df, aes(x = Pathologic_stage, y = Non_silent_Mutation_frequency)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.4, aes(color = ICR_cluster)) +
  scale_color_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
  #scale_color_manual(values = MANTIS_colors) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  ylab("Nonsynonymous mutational load") +
  xlab("Pathologic stage")

png("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/Nonsynonymous_mut_load_by_Stage.png",
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()

#### Alternative #####

fit = lm(Neoantigen_count ~ Non_silent_Mutation_frequency, data = frequency_df)
summary(fit)
print(fit)
# Formula: Neoantigen_count = -2.38770 + (0.09171  * Non_silent_Mutation_frequency)
# For nonhypermutated formula: Neoantigen_count = -0.61144 + (0.08724  * Non_silent_Mutation_frequency)
# For hypermutated formula: Neoantigen_count = -0.61144 + (0.08724  * Non_silent_Mutation_frequency)
if(subset == ""){
  frequency_df$Expected_neoantigen_count = -2.38770 + (0.09171 * frequency_df$Non_silent_Mutation_frequency)}
if(subset == "nonhypermutated"){
  frequency_df$Expected_neoantigen_count = -0.61144 + (0.08724 * frequency_df$Non_silent_Mutation_frequency)}
if(subset == "hypermutated"){
  frequency_df$Expected_neoantigen_count = -11.88801 + (0.09503 * frequency_df$Non_silent_Mutation_frequency)}

frequency_df$Immunoediting_score = (frequency_df$Neoantigen_count / frequency_df$Expected_neoantigen_count) #/ frequency_df$Non_silent_Mutation_frequency

plot = ggplot(data = frequency_df, aes(x = ICRscore, y = Immunoediting_score)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  #scale_y_log10() +
  theme_bw() +
  xlab("ICR score") +
  ylab("Immunoediting score") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

plot(plot)
###########################

