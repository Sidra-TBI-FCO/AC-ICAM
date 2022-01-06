
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c()
ipak(required.packages)

# Set parameters
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = ""
included_stages = c(1, 2, 3)
Surv.cutoff.years = 10
exclude = c("Conpair_lower_90_percent", "non-epithelial") #"Conpair_lower_90_percent" 
# c("Conpair_lower_90_percent", "non-epithelial")


# Load data
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_Immune_gene_signatures_table.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

## Survival Analysis with all continuous variable in seperate script:
# Exclude adjuvant treated
if(exclude_adjuvant == "adjuvant_treated_excluded"){
  clinical_data = clinical_data[which(clinical_data$Adjuvant_treatment == "No"),]
}

# Exclude Stage IV
if(exclude_stage_IV == "exclude_stage_IV" ){
  clinical_data = clinical_data[-which(clinical_data$ajcc_pathologic_tumor_stage == 4),]
}

if(included_stages == "All"){}else{
  clinical_data = clinical_data[which(clinical_data$ajcc_pathologic_tumor_stage %in% included_stages),]
}

if(exclude == ""){
}else{
  clinical_data = clinical_data[-which(clinical_data$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
}

clinical_data = clinical_data[-which(is.na(clinical_data$DSS.Status)),]

rownames(clinical_data) = clinical_data$Patient_ID
rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1, 3)

clinical_data = merge(clinical_data, immune_sig_df, by = "row.names")

for (i in 1:ncol(immune_sig_df)){
  col = colnames(immune_sig_df)[i]
  immune_sig_df[, col] = (immune_sig_df[, col] - min(immune_sig_df[, col]))/(max(immune_sig_df[,col])-min(immune_sig_df[,col]))
}


HR_table = data.frame(Signature = colnames(immune_sig_df), p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)
i=1
for (i in 1:ncol(immune_sig_df)){
  Group.of.interest = colnames(immune_sig_df)[i]
  Y = Surv.cutoff.years * 365
  # time / event object creation
  TS.Alive = clinical_data[clinical_data$DSS.Status == "Alive", c("DSS.Status", "DSS.Time", Group.of.interest)]
  colnames(TS.Alive) = c("Status","Time", Group.of.interest)
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$DSS.Status == "Dead", c("DSS.Status", "DSS.Time", Group.of.interest)]
  colnames(TS.Dead) = c("Status","Time", Group.of.interest)
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  HR_table$p_value[which(HR_table$Signature == Group.of.interest)] = p_value
  HR_table$CI_lower[which(HR_table$Signature == Group.of.interest)] = CI_lower
  HR_table$CI_upper[which(HR_table$Signature == Group.of.interest)] = CI_upper
  HR_table$HR[which(HR_table$Signature == Group.of.interest)] = HR
}

dir.create("./Analysis/Trimmed_p/013_Immune_signatures_survival", showWarnings = FALSE)
save(HR_table, file = paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_DSS_HR_table_exclude_",
                             str_c(exclude, collapse = "_"),
                             "_included_stages_", str_c(included_stages, collapse = "_"),
                             exclude_adjuvant, exclude_stage_IV, "_cutoff_", Surv.cutoff.years,".Rdata"))



