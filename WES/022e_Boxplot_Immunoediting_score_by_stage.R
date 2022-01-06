
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

# Load data
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")
#load("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata")


# Boxplot
plot = ggplot(Merged_dataset, aes(x = ajcc_pathologic_tumor_stage, y = Immunoediting_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Immunoedited), width = 0.1) +
  scale_color_manual(values = c("immunoedited" = "#FF5800", "less immunoedited" = "#8085E9")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  xlab("") +
  ylab("Immunoediting score") +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(c("1", "2"),
                                        c("1", "3"),
                                        c("1", "4"),
                                        c("2", "3"),
                                        c("2", "4")))

dev.new()
plot(plot)
