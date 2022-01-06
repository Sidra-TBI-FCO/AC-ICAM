
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "survminer", "survival", "coin", "survival"))

# Set parameters
Surv.cutoff.years = 20  
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full" 
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
Tissue = "T"


# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))
rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)
abundance_N = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "N")]
colnames(abundance_N) = substring(colnames(abundance_N), 1, 3)
abundance_N = t(abundance_N)

abundance = get(paste("abundance_", Tissue, sep = ""))

load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

clinical_data = clinical_data[which(clinical_data$Patient_ID %in% rownames(abundance)),]

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset$ICR_HML = factor(Merged_dataset$ICR_HML, levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]

ncol(abundance)

i = 5
for (i in 1:ncol(abundance)){
  col = colnames(abundance)[i]
  abundance[, col] = (abundance[, col] - min(abundance[, col]))/(max(abundance[,col])-min(abundance[,col]))
}

#constant_microbiome = which(is.na(colMeans(abundance)))
#abundance = abundance[, -constant_microbiome]
#ncol(abundance)

results = data.frame(Name = colnames(abundance), p_val_logrank = NA, p_val_cox = NA,
                     HR = NA)

i = 1
for(i in 1:ncol(abundance)){
  micr = colnames(abundance)[i]
  # add microbiome variable
  Merged_dataset$Microbiome_var = abundance[,micr][match(Merged_dataset$Patient_ID, rownames(abundance))]
  
  # time / event object creation
  Y = Surv.cutoff.years * 365
  TS.Alive = Merged_dataset[Merged_dataset$DFS.Status == "Disease Free", c("DFS.Status", "DFS.Time", "ICR_HML", 
                                                                      "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var")]
  colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = Merged_dataset[Merged_dataset$DFS.Status == "Event", c("DFS.Status", "DFS.Time", "ICR_HML",
                                                                    "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var")]
  colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Event"
  
  if(length(unique(TS.Surv$Microbiome_var)) == 1){next}
  
  # survival curve
  msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
  mfit = survfit(msurv~TS.Surv$Microbiome_var,conf.type = "log-log")
  
  # Calculations (Needs manual adaptation!)
  mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
  pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
  pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
  
  Cal.Surv = TS.Surv
  mHR = coxph(formula = msurv ~ Cal.Surv$Microbiome_var,data = Cal.Surv)
  
  #logrank_test(Surv(Time, Status) ~ Microbiome_var, data = TS.Surv)
  test = summary(mHR)
  logrank_pval = test$logtest[3]
  cox_pval = test$coefficients[5]
  HR = test$coefficients[2]
  
  results$p_val_logrank[which(results$Name == micr)] = logrank_pval
  results$p_val_cox[which(results$Name == micr)] = cox_pval
  
  results$HR[which(results$Name == micr)] = HR
  
}

results = results[order(results$p_val_cox),]

dir.create("./Analysis/Microbiome/025_Coxph_continuous_Genus_level", showWarnings = FALSE)
write.csv(results, file = "./Analysis/Microbiome/025_Coxph_continuous_Genus_level/025_Coxph_stats_Genus_continuous_PFS.csv",
          row.names = FALSE)