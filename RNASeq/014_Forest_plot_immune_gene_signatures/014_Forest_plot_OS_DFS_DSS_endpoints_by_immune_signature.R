
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = ""
included_stages = c(3)
Surv.cutoff.years = 10
exclude = c("Conpair_lower_90_percent", "non-epithelial")
immune_signature = "ICR |  ICRscore"

# Load data
load(paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_OS",
            "_HR_table_exclude_", str_c(exclude, collapse = "_"),
            "_included_stages_", str_c(included_stages, collapse = "_"),
            exclude_adjuvant, exclude_stage_IV, "_cutoff_", Surv.cutoff.years, ".Rdata"))

HR_table_OS = HR_table

load(paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_DSS",
            "_HR_table_exclude_", str_c(exclude, collapse = "_"),
            "_included_stages_", str_c(included_stages, collapse = "_"),
            exclude_adjuvant, exclude_stage_IV, "_cutoff_", Surv.cutoff.years, ".Rdata"))

HR_table_DSS = HR_table

load(paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_DFS_Def1",
            "_HR_table_exclude_", str_c(exclude, collapse = "_"), 
            "_included_stages_", str_c(included_stages, collapse = "_"),
            exclude_adjuvant, exclude_stage_IV, "_cutoff_", Surv.cutoff.years, ".Rdata"))

HR_table_DFS = HR_table

HR_table = data.frame(Outcome = c("OS", "DSS", "DFS"), 
                      p_value = NA,
                      HR = NA,
                      CI_lower = NA,
                      CI_upper = NA)
HR_table[which(HR_table$Outcome == "OS"), 2:5] = HR_table_OS[which(HR_table_OS$Signature == immune_signature),2:5]
HR_table[which(HR_table$Outcome == "DSS"), 2:5] = HR_table_DSS[which(HR_table_OS$Signature == immune_signature),2:5]
HR_table[which(HR_table$Outcome == "DFS"), 2:5] = HR_table_DFS[which(HR_table_OS$Signature == immune_signature),2:5]


## Forest plot seperate script
n_signatures = nrow(HR_table)
x = n_signatures + 2

HR.matrix = as.matrix(HR_table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

#HR_table = HR_table[order(HR_table$HR),]

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR_table$HR[1:n_signatures]), NA),
    lower = c(NA,HR_table$CI_lower[c(1:n_signatures)], NA),
    upper = c(NA,HR_table$CI_upper[c(1:n_signatures)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")


HR_table$p_value = signif(HR_table$p_value, 3)
HR_table$HR = signif(HR_table$HR, 3)
tabletext<-cbind(
  c("Outcome", as.character(HR_table$Outcome)[c(1:n_signatures)]),
  c("p-value", HR_table$p_value[c(1:n_signatures)]),
  c("HR",      HR_table$HR[c(1:n_signatures)]))


dir.create("./Figures/Trimmed_p/014_Forest_plot_immune_signatures", showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/014_Forest_plot_immune_signatures/014_OS_DSS_DFS", showWarnings = FALSE)
pdf(file = paste0("./Figures/Trimmed_p/014_Forest_plot_immune_signatures/014_OS_DSS_DFS/014_Forest_plot_OS_DSS_DFS_",
                  gsub(" |  ", "_", immune_signature), "_exclude_", 
                  str_c(exclude, collapse = "_"),
                  "_included_stages_", str_c(included_stages, collapse = "_"),
                  exclude_adjuvant, exclude_stage_IV,
                  "_cutoff_", Surv.cutoff.years, ".pdf"),
    height = 1.5, width = 4)

forestplot(mean = HR.matrix[,"HR"],
           lower = HR.matrix[,"CI_lower"],
           upper = HR.matrix[,"CI_upper"],
           labeltext = tabletext[-1,],
           new_page = FALSE,
           zero = 1,
           #is.summary=c(TRUE,rep(FALSE,n.cells),TRUE,rep(FALSE,n.cells),TRUE,FALSE),
           clip=c(0.001,55),
           xlog=TRUE,
           #xlim = c(0, 4),
           boxsize = .25,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 7), xlab = gpar(fontsize = 7),
                            ticks = gpar(fontsize = 10))
)
dev.off()

