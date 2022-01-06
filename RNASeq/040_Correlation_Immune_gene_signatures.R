
# Calculate pearson correlation between all immune gene signature scores in Sidra-LUMC dataset

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("RColorBrewer", "forestplot", "corrplot")
ipak(required.packages)

# Set parameters
test = "pearson"

# Load data
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")

# Correlation
ssGSEA_cor <- cor(immune_sig_df,method=test)

# Correlation significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

ssGSEA_cor_sign <- cor.mtest(ssGSEA_cor, 0.95)

cex.before <- par("cex")
par(cex = 1)
lims=c(-1,1)
if (length(ssGSEA_cor[ssGSEA_cor<0]) == 0) {lims=c(0,1)}

#annotation = data.frame (gene = rownames(Hallmark_GSEA_cor),color = c(rep("#CC0506",54)),stringsAsFactors = FALSE)
#annotation$color[annotation$gene]
#annotation = annotation[corrMatOrder(Hallmark_GSEA_cor,order="FPC"),]

mean_correlation = round(mean(ssGSEA_cor),2)

dir.create("./Figures/Trimmed_p/040_Correlation_Immune_signatures", showWarnings = FALSE)
png("./Figures/Trimmed_p/040_Correlation_Immune_signatures/April_2021_Corrpot_hclust.png", res = 600,
    units = "in", width = 40, height = 40)
corrplot.mixed (ssGSEA_cor,
                #type="lower",
                #p.mat = ICR_cor_sign[[1]],                                                                      # add significance to correlations
                #col = colpattern,
                lower = "square",
                upper ="square",
                order="hclust",
                cl.lim=lims,                                                                                               # only positive correlations
                tl.pos ="lt",
                tl.col = "black" , #as.character(annotation$color),
                insig= "pch",                                                                                              # remove insignificant correlations
                pch = "x",
                pch.cex= 4,
                tl.cex = 1/par("cex"),
                cl.cex = 1/par("cex"),
                cex.main = 1/par("cex"),
                mar=c(6,4.1,7,5))
title(main = list(paste0( " correlation between immune gene signatures. \n ","Mean: ", mean_correlation,"."),
                  cex = 3), line = -6)
par(cex = cex.before)
dev.off()

dir.create("./Analysis/Trimmed_p/040_Correlation_Immune_signatures", showWarnings = FALSE)
save(ssGSEA_cor, ssGSEA_cor_sign, file = paste0("./Analysis/Trimmed_p/040_Correlation_Immune_signatures/Immune_signatures_correlation_", test, ".Rdata"))

png(paste0("./Figures/Trimmed_p/048_Immune_traits_Pearson_correlation/048_2021_April_Pearson_correlation_heatmap_",
           Immune_only, ".png"), res = 600, 
    width = 13, height = 13, units = "in")
hm = heatmap.2(mat_cor,
               # dendrogram control
               Rowv = TRUE,
               distfun = dist,
               hclustfun = hclust,
               dendrogram = c("both"),
               col = "bluered",
               trace = "none",
               margins=c(17,17),
               key = "T", density = "none")
dev.off()