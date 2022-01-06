
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
required.packages = c("ggplot2", "ggpubr", "graphics")
ipak(required.packages)

# Set parameters
Gene1 = "ITGAE"
Gene2 = "TGFB1" #"IFNG"
#Gene3 = "LILRB2" #"KIR2DL4"
#Gene4 = "KIR2DL4" #"LILRB1"
#Gene5 = "HLA-A"

CMS_subset = "All" # "All" or "CMS1"
MSI_subset = "All" # "All" or "MSI-H"

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Plot
plot_df = data.frame(Sample = colnames(RNASeq.QN.LOG2), Gene1 = RNASeq.QN.LOG2[Gene1,], Gene2= RNASeq.QN.LOG2[Gene2,],
                     CMS = NA)
                     #Gene3 = RNASeq.QN.LOG2[Gene3,], Gene4= RNASeq.QN.LOG2[Gene4,])

Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"
plot_df$CMS = Rfcms$RF.predictedCMS[match(plot_df$Sample, rownames(Rfcms))]

plot_df$MSI = MANTIS$MSI[match(substring(plot_df$Sample, 1, 4), MANTIS$Sample_ID)]

if(CMS_subset == "All"){}else{
  plot_df = plot_df[which(plot_df$CMS == CMS_subset),]
}

if(MSI_subset == "All"){}else{
  plot_df = plot_df[which(plot_df$MSI == MSI_subset),]
}

plot = ggplot(plot_df, aes(x = Gene1, y = Gene2)) +
  geom_point(size = 0.4) +
  #geom_point(size = 0.4, aes(color = CMS)) +
  #scale_color_manual(values = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
   #                             "mixed" = "grey")) +
  xlab(Gene1) +
  ylab(Gene2) +
  theme_bw() +
  theme(aspect.ratio = 1/1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", size = 3)
  

dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/03_Scatterplot_multiple_genes", showWarnings = FALSE)
png(paste0("./Figures/Trimmed_p/Gene_of_interest_plots/03_Scatterplot_multiple_genes/GOI_03_", Gene1, "_", Gene2,
           "_in_specific_", CMS_subset, "_", MSI_subset, ".png"),
    res = 600, width = 4, height = 3, units = "in")
plot(plot)
dev.off()


## Correlation matrix
reg <- function(x, y, col) abline(lm(y~x), col=col)

upper.panel<-function(x, y, col.smooth = "red"){
  #points(x,y, pch=19, cex = 0.1, col=c("orange", "grey", "purple")[ICR_cluster_assignment_allcancers$ED])
  points(x, y, pch=19, cex = 0.1, col = "black")
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

#####
reg <- function(x, y, col) abline(lm(y~x), col=col) 

panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                      cex = 0.1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.1, font = 4)
}
colnames(plot_df) = c("Sample", Gene1, Gene2, Gene3, Gene4)

dev.new()
pairs(plot_df[, 2:5], 
      upper.panel = panel.lm,
      lower.panel = panel.cor)
