
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "dplyr")
ipak(required.packages)

# Set parameters
gene_of_interest = "MUC2"

# Load data
load(paste0("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata"))
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")

# Data preparation
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"
plot_df = frequency_df
plot_df$CMS = Rfcms$RF.predictedCMS[match(plot_df$Patient_ID, substring(rownames(Rfcms), 1, 3))]

Maf_gene = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene_of_interest),]
Mut_patients = unique(Maf_gene$Patient_ID)

plot_df$Mutation = "WT"
plot_df$Mutation[which(plot_df$Patient_ID %in% Mut_patients)] = "MUT"
table(plot_df$Mutation)
plot_df$Mutation = factor(plot_df$Mutation, levels = c("WT", "MUT"))

# Stacked barplot
DF1 <- plot_df %>%
  group_by(CMS, Mutation) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot = ggplot(DF1, aes(x = CMS, y =perc*100, fill = Mutation)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = paste0(gene_of_interest, " mutation"), face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_text(size = 19, colour = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values=  c("WT" = "lightblue",  "MUT" = "orange"))

dir.create("./Figures/WES/007_Gene_of_interest_by_CMS", showWarnings = FALSE)
png(paste0("./Figures/WES/007_Gene_of_interest_by_CMS/", gene_of_interest, "_by_CMS_stacked_barplot.png"),
    width = 5.5, height = 5, units = "in", res = 600)
plot(plot)
dev.off()
