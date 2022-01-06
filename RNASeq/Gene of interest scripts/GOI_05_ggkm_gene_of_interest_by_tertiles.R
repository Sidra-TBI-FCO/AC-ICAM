
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 10                                        # SET cut-off
Group.of.interest = "CDH1"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
within_ICR_group = ""
exclude_adjuvant = ""
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Create folders and log file
dir.create("./Figures/Trimmed_p/",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Gene_of_interest_plots/05_Kaplan_Meier_plots", showWarnings = FALSE)

# Load data
load(paste0("./Analysis/Trimmed_p/Gene_of_interest_tertiles/", Group.of.interest, "_tertiles_within", within_ICR_group,".Rdata"))
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

table_cluster_assignment$Patient_ID = rownames(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")

if(exclude == ""){
}else{
  Merged_dataset = Merged_dataset[-which(Merged_dataset$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
}

Merged_dataset$Expression_category = df$Expression_category[match(Merged_dataset$Patient_ID, substring(df$Sample, 1, 3))]
Merged_dataset$Expression_category = factor(Merged_dataset$Expression_category, levels = c("High", "Medium", "Low"))

# Exclude adjuvant treated
if(exclude_adjuvant == "exclude_adjuvant"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$Adjuvant_treatment == "No"),]
}

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$vital_status == "Alive", c("vital_status", "last_contact_days_to", "Expression_category")]
colnames(TS.Alive) = c("Status","Time", "Expression_category")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$vital_status == "Dead", c("vital_status", "death_days_to", "Expression_category")]
colnames(TS.Dead) = c("Status","Time", "Expression_category")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv[,"Expression_category"],conf.type = "log-log")

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
#Cal.Surv[,"Expression_category"] = as.factor(Cal.Surv[,"Expression_category"])
#Cal.Surv[,"Expression_category"] = relevel(Cal.Surv[,"Expression_category"], "ICR3")
Cal.Surv[,"Expression_category"] = as.factor(Cal.Surv[,"Expression_category"])
mHR = coxph(formula = msurv ~ Cal.Surv[,"Expression_category"],data = Cal.Surv, subset = Cal.Surv[, "Expression_category"] %in% c("High", "Low"))
mHR.extract = extract.coxph(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)


# plots
png(paste0("./Figures/Trimmed_p/Gene_of_interest_plots/05_Kaplan_Meier_plots/", exclude_adjuvant,"Overall_Survival_", Group.of.interest, "_within_", 
           within_ICR_group, "_Surv_cutoff_years_",Surv.cutoff.years,".png"),res=600,height=6,width=8,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Expression_category"]),
     ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     xlabs = "Time in months",
     palette = c("red", "green", "blue"),
     PLOT_P = signif(p[2],3),
     PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3),
     PLOT_CI1 = CI[2,1],
     PLOT_CI2 = CI[2,2])
dev.off()

# Cox regression (continuous)
uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
summary(uni_variate_ICRscore)

# Multivariate
multivariate_ICR_and_stage = coxph(formula = Surv(Time, Status) ~ ICRscore + pathologic_stage, data = TS.Surv)
summary(multivariate_ICR_and_stage)

