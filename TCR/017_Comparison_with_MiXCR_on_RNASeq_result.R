

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Load data
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
MiXCR = read.csv("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Diversity.csv", 
                     stringsAsFactors = FALSE)

# Set parameters
Type = "TRB" # ALL IGH IGK IGL TRA TRB TRD TRG 

# Filter type
MiXCR$Type = substring(MiXCR$X, 8, 10)
MiXCR = MiXCR[which(MiXCR$Type == Type),]

# Prepare data
colnames(MiXCR)[1] = "Sample_ID"
MiXCR$Patient_ID = substring(MiXCR$Sample_ID, 1, 3)
MiXCR$Tissue = substring(MiXCR$Sample_ID, 4, 6)

# Filter only tumor samples
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
MiXCR = MiXCR[which(MiXCR$Tissue == "T-P"),]

# Append to overview
TCR_Overview$MiXCR_Pielou = MiXCR$Pielou[match(TCR_Overview$Patient_ID, MiXCR$Patient_ID)]
TCR_Overview$MiXCR_Clonality = 1-TCR_Overview$MiXCR_Pielou
TCR_Overview$MiXCR_Shannon = MiXCR$Shannon[match(TCR_Overview$Patient_ID, MiXCR$Patient_ID)]

plot = ggplot(TCR_Overview, aes(ICRscore, MiXCR_Clonality)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  xlab("ICR score") +
  ylab("TCR Clonality") +
  #xlab("ICR score (from RNASeq)") +
  #ylab(paste0("1- Pielou's Evenness", Type, "\n(from MiXCR on RNASeq)")) +
  #ylim(0, 0.4) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(color = "black", size = 20),
        axis.title.y = element_text(color = "black", size = 20),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson", size = 6,
           label.y = 0.4) +
  geom_smooth(method="lm")

png(paste0("./Figures/TCR/017_Comparison_with_MiXCR_on_RNASeq/1-Pielous_Eveness_", Type, "_and_ICRscore.png"),
    #res = 600, width = 3.2, height = 3, units = "in")
    res = 600, width = 4, height = 4, units = "in")
plot(plot)
dev.off()

plot = ggplot(TCR_Overview, aes(ICRscore, MiXCR_Pielou)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  xlab("ICR score") +
  ylab("Pielou's Evenness") +
  #xlab("ICR score (from RNASeq)") +
  #ylab(paste0("1- Pielou's Evenness", Type, "\n(from MiXCR on RNASeq)")) +
  #ylim(0, 0.4) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(color = "black", size = 20),
        axis.title.y = element_text(color = "black", size = 20),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson", size = 6,
           label.y = 0.6) +
  geom_smooth(method="lm")

png(paste0("./Figures/TCR/017_Comparison_with_MiXCR_on_RNASeq/Pielous_Eveness_", Type, "_and_ICRscore.png"),
    #res = 600, width = 3.2, height = 3, units = "in")
    res = 600, width = 4, height = 4, units = "in")
plot(plot)
dev.off()

plot = ggplot(TCR_Overview, aes(ICRscore, MiXCR_Shannon)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  xlab("ICR score") +
  ylab("TCR Shannon") +
  #xlab("ICR score (from RNASeq)") +
  #ylab(paste0("1- Pielou's Evenness", Type, "\n(from MiXCR on RNASeq)")) +
  #ylim(0, 0.4) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(color = "black", size = 20),
        axis.title.y = element_text(color = "black", size = 20),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson", size = 6,
           label.y = 0.6) +
  geom_smooth(method="lm")

png(paste0("./Figures/TCR/017_Comparison_with_MiXCR_on_RNASeq/TCR_Shannon_", Type, "_and_ICRscore.png"),
    #res = 600, width = 3.2, height = 3, units = "in")
    res = 600, width = 4, height = 4, units = "in")
plot(plot)
dev.off()


plot = ggplot(TCR_Overview, aes(productive_clonality, MiXCR_Pielou)) +
  geom_point(aes(color = HLM_cluster), size = 0.8) +
  xlab("Productive clonality \n(from Adaptive ImmunoSeq Assay \n using DNA)") +
  ylab(paste0("Pielou's Evenness ", Type, "\n(from MiXCR on RNASeq)")) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  stat_cor(method = "pearson", size = 5,
           label.y = 0.7) +
  geom_smooth(method="lm")

dir.create("./Figures/TCR/017_Comparison_with_MiXCR_on_RNASeq", showWarnings = FALSE)
png(paste0("./Figures/TCR/017_Comparison_with_MiXCR_on_RNASeq/Pielous_Eveness_", Type, "_and_Productive_Clonality.png"),
    #res = 600, width = 3.2, height = 3, units = "in")
    res = 600, width = 4, height = 4, units = "in")
plot(plot)
dev.off()

plot2 = ggplot(TCR_Overview, aes(productive_clonality, MiXCR_Clonality)) +
  geom_point(size = 0.8) +
  #xlab("Productive clonality \n(from Adaptive ImmunoSeq Assay \n using DNA)") +
  xlab("Productive clonality") +
  ylab(paste0("MiXCR clonality")) +
  #ylab(paste0("1- Pielou's Evenness ", Type, "\n(from MiXCR on RNASeq)")) +
  #scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.position = "none",
        aspect.ratio = 1/1) +
  labs(color = "ICR cluster") +
  #stat_cor(method = "pearson", size = 4.5,
   #        label.y = 0.4) +
  geom_smooth(method="lm")

png(paste0("./Figures/TCR/017_Comparison_with_MiXCR_on_RNASeq/1-Pielous_Eveness_", Type, "_and_Productive_Clonality.png"),
    #res = 600, width = 3.2, height = 3, units = "in")
    res = 600, width = 3, height = 3, units = "in")
plot(plot2)
dev.off()

