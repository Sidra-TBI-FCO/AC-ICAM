
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "survminer", "survival", "coin"))

# Set parameters
Surv.cutoff.years = 20  
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Phylum" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full" 
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
Tissue = "T"
micr = "D_1__Fusobacteria" # name of the species/phylum of interest

# Load data
load(file= paste0("./Analysis/Microbiome/020_Categories_by_median/020_Categorized_", Type, "_abundance_by_median_", Tissue, "_", Rank,".Rdata"))
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

clinical_data = clinical_data[which(clinical_data$Patient_ID %in% rownames(abundance_cat)),]

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset$ICR_HML = factor(Merged_dataset$ICR_HML, levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]

# add microbiome variable
Merged_dataset$Microbiome_var = abundance_cat[,micr][match(Merged_dataset$Patient_ID, rownames(abundance_cat))]

# add combined microbiome and ICR variable
Merged_dataset$ICR_micr = paste(Merged_dataset$ICR_HML, Merged_dataset$Microbiome_var, sep = "+")
Merged_dataset = Merged_dataset[-which(Merged_dataset$ICR_HML == "ICR Medium"),] # exclude ICR medium
table(Merged_dataset$ICR_micr)

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$vital_status == "Alive", c("vital_status", "last_contact_days_to", "ICR_HML", 
                                                                    "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var",
                                                                    "ICR_micr")]
colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var", "ICR_micr")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$vital_status == "Dead", c("vital_status", "death_days_to", "ICR_HML",
                                                                  "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var",
                                                                  "ICR_micr")]
colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var", "ICR_micr")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"

TS.Surv$ICR_micr = factor(TS.Surv$ICR_micr, levels = c("ICR Low+Low", "ICR Low+High",
                                                             "ICR High+Low", "ICR High+High"))

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv$ICR_micr,conf.type = "log-log")

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
mHR = coxph(formula = msurv ~ Cal.Surv$ICR_micr,data = Cal.Surv, subset = Cal.Surv$ICR_micr %in% c("ICR High+Low", "ICR High+High"))

logrank_test(Surv(Time, Status) ~ ICR_micr, data = TS.Surv)
summary(mHR)

# plots

ggsurv_plot = ggsurvplot(mfit,
                         xlim = c(0, 90),
                         break.time.by = 30,
                         data = TS.Surv,
                         censor = TRUE,
                         risk.table = TRUE,
                         tables.y.text.col = FALSE,
                         tables.y.text = FALSE,
                         tables.height = 0.3,
                         tables.theme = theme_cleantable(),
                         tables.col = "strata",
                         risk.table.pos = "out",
                         legend = "none",
                         ylab = "",
                         xlab = "Time in months",
                         fontsize = 4.5,
                         font.x = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #pval = TRUE,
                         palette = c("darkred", "red", "darkblue", "blue")
)

dir.create("./Figures/Microbiome/022_KM_relative_abundance_and_ICR_groups", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/022_KM_relative_abundance_and_ICR_groups/KM_ICR_cluster_and_", Type, "_abundance_in_", Tissue, "_tissue_", Rank,
           "_", micr, ".png"),
    res=600,height=3.8,width=4.2,unit="in")  # set filename
print(ggsurv_plot)
dev.off()
